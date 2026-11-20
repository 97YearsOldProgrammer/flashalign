// The persistent thread pool, parallel_for, the serial pipeline, the chunk stream and thread names.
#pragma once

#include <pthread.h>

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <deque>
#include <exception>
#include <fstream>
#include <functional>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <utility>
#include <vector>

/**** thread_name ****/
// Thread names for `top -H` and profilers. Linux keeps at most 15 characters; errors are
// ignored.

namespace fa { namespace cpu { namespace threading {

inline void name_current_thread(const char* name) {
#if defined(__APPLE__)
    (void)pthread_setname_np(name);
#elif defined(__linux__)
    (void)pthread_setname_np(pthread_self(), name);
#else
    (void)name;
#endif
}

// "prefix-N", e.g. "fa-map-17".
inline void name_current_thread(const char* prefix, int index) {
    char buf[16];
    std::snprintf(buf, sizeof buf, "%s-%d", prefix, index);
    name_current_thread(buf);
}

}}}  // namespace fa::cpu::threading

/**** parallel_for ****/
// parallel_for over a persistent thread pool. The caller runs as worker 0, items are claimed
// from one atomic counter, and the first exception is rethrown on the caller. A nested or
// concurrent call runs on a temporary team of its own.

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

/**** pipeline ****/
// The three stages on the calling thread (`-t 1`): produce(), process() and consume() in turn,
// in input order.

namespace fa { namespace cpu { namespace threading {

template <class Item, class Result, class Produce, class Process, class Consume>
void run_pipeline_serial(Produce produce, Process process, Consume consume) {
    for (;;) {
        std::optional<Item> item = produce();
        if (!item) break;
        consume(process(std::move(*item)));
    }
}

}}}  // namespace fa::cpu::threading

/**** chunk_stream ****/
// An ordered stream of chunks under one bound on the weight alive. A reader thread waits for
// room, then submits a chunk; a writer thread collects the chunks oldest first and releases
// each one's weight once it is written. A chunk's weight is alive from submit() to release(),
// and the reader waits before it reads, so the weight alive stays under the bound plus one
// chunk. A worker that finds the oldest chunk's items all claimed moves on to the next one,
// so one slow item delays only its own chunk.
//
// The stream owns `n_threads` workers with stable tids in [0, n_threads). Other threads may
// lend time through run_one(), one item per call.
//
// A chunk is quiescent once every item is claimed and no worker is inside it. A worker
// leaves a chunk by decrementing its count under the mutex, so a collector that sees zero
// also sees all of the workers' writes. The first exception in a chunk stops that chunk's
// remaining items and is rethrown by its collect(); later chunks still run.

namespace fa { namespace cpu { namespace threading {

template <class Payload>
class ChunkStream {
public:
    // body(item, tid): item in [0, n_items), tid the worker's slot.
    using Body = std::function<void(std::int64_t, int)>;

    // `max_alive` bounds the weight alive. Both arguments are clamped to at least 1.
    ChunkStream(int n_threads, std::uint64_t max_alive)
        : n_threads_(n_threads < 1 ? 1 : n_threads),
          max_alive_(max_alive < 1 ? 1 : max_alive) {
        workers_.reserve(static_cast<std::size_t>(n_threads_));
        for (int tid = 0; tid < n_threads_; ++tid)
            workers_.emplace_back([this, tid] { worker_loop(tid); });
    }

    ChunkStream(const ChunkStream&) = delete;
    ChunkStream& operator=(const ChunkStream&) = delete;

    // Stops new claims, waits for running items (owned workers are joined, run_one callers
    // waited out) and drops uncollected payloads. No thread may still be in wait_room() or
    // collect(): the caller joins its reader and ends its writer first.
    ~ChunkStream() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            shutdown_ = true;
            for (auto& chunk : chunks_)
                chunk->claim.store(chunk->n, std::memory_order_relaxed);
        }
        work_.notify_all();
        for (auto& worker : workers_)
            worker.join();
        std::unique_lock<std::mutex> lock(mutex_);
        done_.wait(lock, [this] { return helpers_in_flight_ == 0; });
    }

    // Reader thread. Blocks while the weight alive is at or over the bound; false after
    // abort().
    bool wait_room() {
        std::unique_lock<std::mutex> lock(mutex_);
        room_.wait(lock, [this] { return aborted_ || alive_ < max_alive_; });
        return !aborted_;
    }

    // Reader thread, before a wait_room(): whether the stream is short of input. It turns
    // short under three quarters of the bound and stays short until the bound is reached.
    bool short_of_input() {
        std::lock_guard<std::mutex> lock(mutex_);
        if (alive_ >= max_alive_)
            short_ = false;
        else if (alive_ < max_alive_ - max_alive_ / 4)
            short_ = true;
        return short_;
    }

    // Reader thread. Queues a chunk of `n_items` weighing `weight`. Items are claimed in
    // index order, so items sorted longest first give LPT scheduling. Dropped after abort().
    void submit(std::int64_t n_items, std::uint64_t weight, Body body,
                Payload payload) {
        std::unique_ptr<Chunk> chunk(new Chunk(
            n_items < 0 ? 0 : n_items, std::move(body), std::move(payload)));
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (aborted_)
                return;
            alive_ += weight;
            chunks_.push_back(std::move(chunk));
        }
        work_.notify_all();
        // An empty chunk is collectable at once, and the writer may wait on an empty stream.
        done_.notify_all();
    }

    // Reader thread: the end of the input.
    void close() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            closed_ = true;
        }
        done_.notify_all();
    }

    // Writer thread. Waits until the oldest chunk is quiescent, removes it and returns its
    // payload, or rethrows its first item exception; its weight stays alive until release().
    // std::nullopt once closed and empty, or after abort().
    std::optional<Payload> collect() {
        std::unique_ptr<Chunk> chunk;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            done_.wait(lock, [this] {
                return aborted_ || (chunks_.empty() ? closed_
                                                    : quiescent(*chunks_.front()));
            });
            if (aborted_ || chunks_.empty())
                return std::nullopt;
            chunk = std::move(chunks_.front());
            chunks_.pop_front();
        }
        if (chunk->first_exception)
            std::rethrow_exception(chunk->first_exception);
        return std::move(chunk->payload);
    }

    // Writer thread: a collected chunk of `weight` is written. Wakes the reader.
    void release(std::uint64_t weight) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (weight > alive_) {
                throw std::logic_error(
                    "ChunkStream::release of more than is alive");
            }
            alive_ -= weight;
        }
        room_.notify_all();
    }

    // Any thread, on an error: wait_room() and collect() return at once, and unclaimed
    // items are dropped.
    void abort() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            aborted_ = true;
            for (auto& chunk : chunks_)
                chunk->claim.store(chunk->n, std::memory_order_relaxed);
        }
        room_.notify_all();
        done_.notify_all();
    }

    // Runs at most one item on a thread the stream does not own. Returns false when no chunk
    // has an unclaimed item; true otherwise, even if another thread took the last item first.
    // `tid` is only passed to the body and need not lie in [0, n_threads). Callers must
    // stop calling before the stream is destroyed.
    bool run_one(int tid) {
        Chunk* chunk = nullptr;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (shutdown_)
                return false;
            chunk = oldest_claimable();
            if (chunk == nullptr)
                return false;
            ++chunk->workers;      // pins the chunk against collection
            ++helpers_in_flight_;  // pins the stream against destruction
        }
        const std::int64_t item =
            chunk->claim.fetch_add(1, std::memory_order_relaxed);
        if (item < chunk->n) {
            try {
                chunk->body(item, tid);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_);
                if (!chunk->first_exception)
                    chunk->first_exception = std::current_exception();
                chunk->claim.store(chunk->n, std::memory_order_relaxed);
            }
        }
        bool became_quiescent = false;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            --chunk->workers;
            became_quiescent = quiescent(*chunk);
            if (--helpers_in_flight_ == 0 && shutdown_) {
                // Notify under the lock: the waiting destructor may destroy done_ as soon as
                // the mutex is released.
                done_.notify_all();
                return true;
            }
        }
        // `chunk` may already be collected and destroyed.
        if (became_quiescent)
            done_.notify_all();
        return true;
    }

private:
    struct Chunk {
        Chunk(std::int64_t items, Body item_body, Payload item_payload)
            : n(items), body(std::move(item_body)),
              payload(std::move(item_payload)) {}

        std::int64_t n = 0;
        // Next item to hand out; set to `n` to stop claims, and may exceed n.
        std::atomic<std::int64_t> claim{0};
        int workers = 0;                        // guarded by mutex_
        std::exception_ptr first_exception;     // guarded by mutex_
        Body body;
        Payload payload;
    };

    // Called with mutex_ held. `claim` may be read relaxed: claimers join `workers` under the
    // mutex before touching `claim` and leave after their last write.
    static bool quiescent(const Chunk& chunk) {
        return chunk.claim.load(std::memory_order_relaxed) >= chunk.n &&
               chunk.workers == 0;
    }

    // Called with mutex_ held: the oldest chunk with an unclaimed item, or nullptr.
    Chunk* oldest_claimable() {
        for (auto& chunk : chunks_) {
            if (chunk->claim.load(std::memory_order_relaxed) < chunk->n)
                return chunk.get();
        }
        return nullptr;
    }

    // Claims and runs items until the chunk has none left.
    void drain(Chunk& chunk, int tid) {
        for (;;) {
            const std::int64_t item =
                chunk.claim.fetch_add(1, std::memory_order_relaxed);
            if (item >= chunk.n)
                return;
            try {
                chunk.body(item, tid);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_);
                if (!chunk.first_exception)
                    chunk.first_exception = std::current_exception();
                chunk.claim.store(chunk.n, std::memory_order_relaxed);
                return;
            }
        }
    }

    void worker_loop(int tid) {
        name_current_thread("fa-map", tid);
        for (;;) {
            Chunk* chunk = nullptr;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                for (;;) {
                    if (shutdown_)
                        return;
                    chunk = oldest_claimable();
                    if (chunk != nullptr) {
                        ++chunk->workers;  // pins it against collection
                        break;
                    }
                    work_.wait(lock);
                }
            }
            drain(*chunk, tid);
            bool became_quiescent = false;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                --chunk->workers;
                became_quiescent = quiescent(*chunk);
            }
            // `chunk` may already be collected and destroyed.
            if (became_quiescent)
                done_.notify_all();
        }
    }

    std::mutex mutex_;              // guards every member below but the threads and constants
    std::condition_variable work_;  // a chunk can be claimed / teardown began
    std::condition_variable done_;  // a chunk became collectable / the stream ended
    std::condition_variable room_;  // weight was released / abort
    std::deque<std::unique_ptr<Chunk>> chunks_;  // oldest at the front
    std::vector<std::thread> workers_;
    int n_threads_;
    std::uint64_t max_alive_;
    std::uint64_t alive_ = 0;
    bool short_ = true;  // short_of_input()'s last answer
    bool closed_ = false;
    bool aborted_ = false;
    bool shutdown_ = false;
    // run_one callers inside an item. The destructor waits for zero.
    int helpers_in_flight_ = 0;
};

}}}  // namespace fa::cpu::threading
