// parallel_for over a persistent thread pool. The caller runs as worker 0, items are claimed
// from one atomic counter, and the first exception is rethrown on the caller. A nested or
// concurrent call runs on a temporary team of its own.
#ifndef FA_THREADING_PARALLEL_FOR_H
#define FA_THREADING_PARALLEL_FOR_H

#include "threading/thread_name.h"

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <fstream>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace fa { namespace cpu { namespace threading {

// Linux cgroup (v2 or v1) CPU quota rounded up to whole cores, or -1 when there is none.
inline int cgroup_cpu_quota_cores() {
#if defined(__linux__)
    {
        std::ifstream v2("/sys/fs/cgroup/cpu.max");
        std::string quota_tok;
        long long period = 0;
        if (v2 >> quota_tok >> period && quota_tok != "max" && period > 0) {
            try {
                const long long quota = std::stoll(quota_tok);
                if (quota > 0) return static_cast<int>((quota + period - 1) / period);
            } catch (...) {}
        }
    }
    {
        std::ifstream fq("/sys/fs/cgroup/cpu/cpu.cfs_quota_us");
        std::ifstream fp("/sys/fs/cgroup/cpu/cpu.cfs_period_us");
        long long quota = -1, period = 0;
        if (fq >> quota && fp >> period && quota > 0 && period > 0)
            return static_cast<int>((quota + period - 1) / period);
    }
#endif
    return -1;
}

// Hardware concurrency, capped by the cgroup CPU quota so a container is not oversubscribed.
inline int default_thread_count() {
    const unsigned hc = std::thread::hardware_concurrency();
    int n = hc == 0 ? 1 : static_cast<int>(hc);
    const int quota = cgroup_cpu_quota_cores();
    if (quota > 0 && quota < n) n = quota;
    return n < 1 ? 1 : n;
}

namespace detail {

// Claims and runs items for one worker until none remain. The first exception is kept and
// stops further claims.
template <class Body>
inline void drain_counter(std::atomic<int64_t>& counter, int64_t n, int tid,
                          Body&& body, std::mutex& exc_mutex,
                          std::exception_ptr& first_exc) {
    try {
        for (;;) {
            const int64_t i = counter.fetch_add(1, std::memory_order_relaxed);
            if (i >= n) break;
            body(i, tid);
        }
    } catch (...) {
        std::lock_guard<std::mutex> lock(exc_mutex);
        if (!first_exc) first_exc = std::current_exception();
        counter.store(n, std::memory_order_relaxed);
    }
}

// Process-wide worker team, started lazily and joined at exit. It grows to the largest
// thread count requested and never shrinks.
class PersistentPool {
public:
    static PersistentPool& instance() {
        static PersistentPool pool;
        return pool;
    }

    PersistentPool(const PersistentPool&) = delete;
    PersistentPool& operator=(const PersistentPool&) = delete;

    // Runs the job on the team and returns true, rethrowing the first body exception once all
    // workers are idle. Returns false without running anything when the pool is busy.
    bool try_run(int n_threads, int64_t n,
                 const std::function<void(int64_t, int)>& fn) {
        std::unique_lock<std::mutex> submit_lk(submit_mutex_, std::try_to_lock);
        if (!submit_lk.owns_lock()) return false;

        ensure_workers(n_threads - 1);  // tids 1 .. n_threads-1

        {
            std::lock_guard<std::mutex> lk(m_);
            job_fn_ = &fn;
            job_n_ = n;
            job_threads_ = n_threads;
            job_counter_.store(0, std::memory_order_relaxed);
            job_first_exc_ = nullptr;
            remaining_ = n_threads - 1;  // active workers that must report done
            ++generation_;               // publish the job to the workers
        }
        cv_work_.notify_all();

        // The calling thread participates as tid 0.
        drain_counter(job_counter_, n, 0,
                      [this](int64_t i, int tid) { (*job_fn_)(i, tid); },
                      exc_mutex_, job_first_exc_);

        std::exception_ptr e;
        {
            std::unique_lock<std::mutex> lk(m_);
            cv_done_.wait(lk, [&] { return remaining_ == 0; });
            e = job_first_exc_;
            job_fn_ = nullptr;  // no worker dereferences it past this point
        }
        if (e) std::rethrow_exception(e);
        return true;
    }

    // Joins the workers; the next try_run spawns a new team. job_threads_ = 0 matters: a
    // respawned worker starts with seen = 0, wakes on the last job's generation and must skip
    // it rather than call the cleared job_fn_.
    void join_workers() {
        std::unique_lock<std::mutex> submit_lk(submit_mutex_);
        std::vector<std::thread> leaving;
        {
            std::lock_guard<std::mutex> lk(m_);
            shutdown_ = true;
            leaving.swap(workers_);
        }
        cv_work_.notify_all();
        for (auto& th : leaving) th.join();
        std::lock_guard<std::mutex> lk(m_);
        shutdown_ = false;
        job_threads_ = 0;
    }

    ~PersistentPool() {
        {
            std::lock_guard<std::mutex> lk(m_);
            shutdown_ = true;
        }
        cv_work_.notify_all();
        for (auto& th : workers_) th.join();
    }

private:
    PersistentPool() = default;

    // Grows the team to at least `want` workers.
    void ensure_workers(int want) {
        if (want < 0) want = 0;
        std::lock_guard<std::mutex> lk(m_);
        while (static_cast<int>(workers_.size()) < want) {
            const int tid = static_cast<int>(workers_.size()) + 1;  // 1-based
            workers_.emplace_back([this, tid] { worker_loop(tid); });
        }
    }

    void worker_loop(int tid) {
        name_current_thread("fa-pool", tid);
        uint64_t seen = 0;
        for (;;) {
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_work_.wait(lk, [&] {
                    return shutdown_ || generation_ != seen;
                });
                if (shutdown_) return;
                seen = generation_;
                // Workers beyond this job's thread count sit it out.
                if (tid >= job_threads_) continue;
            }
            drain_counter(job_counter_, job_n_, tid,
                          [this](int64_t i, int t) { (*job_fn_)(i, t); },
                          exc_mutex_, job_first_exc_);
            {
                std::lock_guard<std::mutex> lk(m_);
                --remaining_;
            }
            cv_done_.notify_one();
        }
    }

    std::mutex submit_mutex_;  // one job at a time
    std::mutex m_;             // guards the job slot + team state
    std::mutex exc_mutex_;     // guards job_first_exc_ from the body catch path
    std::condition_variable cv_work_;
    std::condition_variable cv_done_;

    std::vector<std::thread> workers_;
    bool shutdown_ = false;
    uint64_t generation_ = 0;

    // Current job slot (valid for generation_ while remaining_ > 0).
    const std::function<void(int64_t, int)>* job_fn_ = nullptr;
    int64_t job_n_ = 0;
    int job_threads_ = 0;
    std::atomic<int64_t> job_counter_{0};
    int remaining_ = 0;
    std::exception_ptr job_first_exc_;
};

// Spawn-and-join team for a nested or concurrent parallel_for.
template <class Func>
inline void legacy_spawn_join(int n_threads, int64_t n, Func& func) {
    std::atomic<int64_t> next{0};
    std::mutex exc_mutex;
    std::exception_ptr first_exc;

    auto body = [&](int64_t i, int tid) { func(i, tid); };
    auto worker = [&](int tid) {
        drain_counter(next, n, tid, body, exc_mutex, first_exc);
    };

    std::vector<std::thread> pool;
    pool.reserve(static_cast<size_t>(n_threads - 1));
    for (int t = 1; t < n_threads; ++t) pool.emplace_back(worker, t);
    worker(0);  // the calling thread is tid 0
    for (auto& th : pool) th.join();

    if (first_exc) std::rethrow_exception(first_exc);
}

}  // namespace detail

// Joins the shared team and releases its threads, for a process done with parallel_for (for
// example after loading the index). A later parallel_for starts a new team. Waits for a
// running job, but a parallel_for started on another thread afterwards respawns the team.
inline void join_persistent_pool() {
    detail::PersistentPool::instance().join_workers();
}

// Runs func(i, tid) for every i in [0, n), with tid in [0, n_threads) for per-thread storage.
// Blocks until all items are done.
template <class Func>
void parallel_for(int n_threads, int64_t n, Func func) {
    if (n <= 0) return;
    if (n_threads < 1) n_threads = 1;
    if (n_threads == 1 || n == 1) {
        for (int64_t i = 0; i < n; ++i) func(i, 0);
        return;
    }
    if (static_cast<int64_t>(n_threads) > n) n_threads = static_cast<int>(n);

    std::function<void(int64_t, int)> fn(
        [&func](int64_t i, int tid) { func(i, tid); });
    if (detail::PersistentPool::instance().try_run(n_threads, n, fn)) return;

    detail::legacy_spawn_join(n_threads, n, func);
}

}}}  // namespace fa::cpu::threading

#endif  // FA_THREADING_PARALLEL_FOR_H
