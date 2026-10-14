// The persistent thread pool, parallel_for, the ordered three-stage pipeline, the batch window and thread names.
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
#include <queue>
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
// Ordered three-stage pipeline: produce() on a reader thread, process() on the caller and
// consume() on a writer thread, joined by bounded FIFO queues.
//
// process() returns std::nullopt while it holds a batch back (a stage with a window of
// batches in flight). At end of input, flush() is called until it returns std::nullopt and
// each value it returns goes to the writer in order.

namespace fa { namespace cpu { namespace threading {

// Single-producer, single-consumer blocking queue. close() ends the stream after the backlog
// drains; abort() makes push and pop fail at once so no stage stays blocked on an error.
template <class T>
class BoundedQueue {
public:
    explicit BoundedQueue(std::size_t capacity)
        : capacity_(capacity < 1 ? 1 : capacity) {}

    // False when the queue is closed or aborted; the producer should stop.
    bool push(T value) {
        std::unique_lock<std::mutex> lock(mutex_);
        not_full_.wait(lock, [&] {
            return aborted_ || closed_ || queue_.size() < capacity_;
        });
        if (aborted_ || closed_) return false;
        queue_.push(std::move(value));
        lock.unlock();
        not_empty_.notify_one();
        return true;
    }

    // std::nullopt once closed and drained, or at once on abort.
    std::optional<T> pop() {
        std::unique_lock<std::mutex> lock(mutex_);
        not_empty_.wait(lock, [&] {
            return aborted_ || !queue_.empty() || closed_;
        });
        if (aborted_ || queue_.empty()) return std::nullopt;
        T value = std::move(queue_.front());
        queue_.pop();
        lock.unlock();
        not_full_.notify_one();
        return value;
    }

    void close() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            closed_ = true;
        }
        not_empty_.notify_all();
        not_full_.notify_all();
    }

    void abort() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            aborted_ = true;
        }
        not_empty_.notify_all();
        not_full_.notify_all();
    }

private:
    std::mutex mutex_;
    std::condition_variable not_full_;
    std::condition_variable not_empty_;
    std::queue<T> queue_;
    std::size_t capacity_;
    bool closed_ = false;
    bool aborted_ = false;
};

namespace detail {

// Default flush for a process stage that holds nothing back.
struct NoFlush {
    std::nullopt_t operator()() const { return std::nullopt; }
};

}  // namespace detail

// Runs the pipeline to completion. `depth` is each queue's capacity in batches. An exception
// in any stage aborts both queues and is rethrown on the caller after the threads join.
template <class Item, class Result, class Produce, class Process, class Consume,
          class Flush = detail::NoFlush>
void run_pipeline(Produce produce, Process process, Consume consume,
                  std::size_t depth = 2, Flush flush = Flush()) {
    BoundedQueue<Item> to_process(depth);
    BoundedQueue<Result> to_consume(depth);

    std::exception_ptr reader_exc;
    std::exception_ptr writer_exc;
    std::exception_ptr worker_exc;

    std::thread reader([&] {
        name_current_thread("fa-reader");
        try {
            for (;;) {
                std::optional<Item> item = produce();
                if (!item) break;
                if (!to_process.push(std::move(*item))) break;  // downstream aborted
            }
            to_process.close();
        } catch (...) {
            reader_exc = std::current_exception();
            to_process.abort();
            to_consume.abort();
        }
    });

    std::thread writer([&] {
        name_current_thread("fa-writer");
        try {
            for (;;) {
                std::optional<Result> result = to_consume.pop();
                if (!result) break;
                consume(std::move(*result));
            }
        } catch (...) {
            writer_exc = std::current_exception();
            to_process.abort();
            to_consume.abort();
        }
    });

    // The calling thread runs process(), so its own worker pool gets the compute threads.
    try {
        // False once the writer has aborted: then nothing is flushed.
        bool downstream_open = true;
        for (;;) {
            std::optional<Item> item = to_process.pop();
            if (!item) break;
            std::optional<Result> result = process(std::move(*item));
            if (!result) continue;  // held back by the process stage
            if (!to_consume.push(std::move(*result))) {
                downstream_open = false;  // downstream aborted
                break;
            }
        }
        // End of input: emit what the process stage still holds, in order.
        while (downstream_open) {
            std::optional<Result> result = flush();
            if (!result) break;
            if (!to_consume.push(std::move(*result))) {
                downstream_open = false;
                break;
            }
        }
        to_consume.close();
    } catch (...) {
        worker_exc = std::current_exception();
        to_process.abort();
        to_consume.abort();
    }

    reader.join();
    writer.join();

    if (worker_exc) std::rethrow_exception(worker_exc);
    if (reader_exc) std::rethrow_exception(reader_exc);
    if (writer_exc) std::rethrow_exception(writer_exc);
}

// All three stages on the calling thread (`-t 1`), in the same output order as run_pipeline.
template <class Item, class Result, class Produce, class Process, class Consume,
          class Flush = detail::NoFlush>
void run_pipeline_serial(Produce produce, Process process, Consume consume,
                         Flush flush = Flush()) {
    for (;;) {
        std::optional<Item> item = produce();
        if (!item) break;
        std::optional<Result> result = process(std::move(*item));
        if (!result) continue;  // held back by the process stage
        consume(std::move(*result));
    }
    // End of input: emit what the process stage still holds, in order.
    for (;;) {
        std::optional<Result> result = flush();
        if (!result) break;
        consume(std::move(*result));
    }
}

}}}  // namespace fa::cpu::threading

/**** batch_window ****/
// A bounded window of ordered batches. Up to `window` batches are in flight; a worker that
// finds the oldest batch exhausted moves on to the next one instead of waiting at a barrier,
// so one slow item delays only its own batch. collect() returns batches oldest first.
//
// The window owns `n_threads` workers with stable tids in [0, n_threads); the calling thread
// only submits and collects. submit(), collect() and the destructor belong to that one
// thread. Other threads may lend time through run_one(), one item per call.
//
// A batch is quiescent once every item is claimed and no worker is inside it. A worker
// leaves a batch by decrementing its count under the mutex, so a collector that sees zero
// also sees all of the workers' writes. The first exception in a batch stops that batch's
// remaining items and is rethrown by its collect(); later batches still run.

namespace fa { namespace cpu { namespace threading {

template <class Payload>
class BatchWindow {
public:
    // body(item, tid): item in [0, n_items), tid the worker's slot.
    using Body = std::function<void(std::int64_t, int)>;

    // Both arguments are clamped to at least 1.
    BatchWindow(int n_threads, std::size_t window)
        : n_threads_(n_threads < 1 ? 1 : n_threads),
          window_(window == 0 ? 1 : window) {
        workers_.reserve(static_cast<std::size_t>(n_threads_));
        for (int tid = 0; tid < n_threads_; ++tid)
            workers_.emplace_back([this, tid] { worker_loop(tid); });
    }

    BatchWindow(const BatchWindow&) = delete;
    BatchWindow& operator=(const BatchWindow&) = delete;

    // Stops new claims, waits for running items (owned workers are joined, run_one callers
    // waited out) and drops uncollected payloads. A body that waits on the destroying
    // thread deadlocks here.
    ~BatchWindow() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            shutdown_ = true;
            for (auto& batch : batches_)
                batch->claim.store(batch->n, std::memory_order_relaxed);
        }
        work_.notify_all();
        for (auto& worker : workers_)
            worker.join();
        std::unique_lock<std::mutex> lock(mutex_);
        done_.wait(lock, [this] { return helpers_in_flight_ == 0; });
    }

    // Queues a batch of `n_items`. Items are claimed in index order, so items sorted longest
    // first give LPT scheduling. Throws std::logic_error when the window is full, since only
    // the calling thread could free a slot.
    void submit(std::int64_t n_items, Body body, Payload payload) {
        std::unique_ptr<Batch> batch(new Batch(
            n_items < 0 ? 0 : n_items, std::move(body), std::move(payload)));
        const bool nothing_to_do = batch->n == 0;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (batches_.size() >= window_) {
                throw std::logic_error(
                    "BatchWindow::submit with the window already full");
            }
            batches_.push_back(std::move(batch));
        }
        work_.notify_all();
        // No worker will touch an empty batch, so wake a waiting collector here.
        if (nothing_to_do) done_.notify_all();
    }

    // Waits for the oldest batch to finish, removes it and returns its payload, or rethrows
    // its first item exception. Precondition: in_flight() > 0.
    Payload collect() {
        std::unique_ptr<Batch> batch;
        {
            std::unique_lock<std::mutex> lock(mutex_);
            if (batches_.empty()) {
                throw std::logic_error(
                    "BatchWindow::collect with no batch in flight");
            }
            done_.wait(lock, [this] { return quiescent(*batches_.front()); });
            batch = std::move(batches_.front());
            batches_.pop_front();
        }
        if (batch->first_exception)
            std::rethrow_exception(batch->first_exception);
        return std::move(batch->payload);
    }

    // Batches submitted but not yet collected, complete or not.
    std::size_t in_flight() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return batches_.size();
    }

    // Runs at most one item on a thread the window does not own. Returns false when no batch
    // has an unclaimed item; true otherwise, even if another thread took the last item first.
    // `tid` is only passed to the body and need not lie in [0, n_threads). Callers must
    // stop calling before the window is destroyed.
    bool run_one(int tid) {
        Batch* batch = nullptr;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (shutdown_)
                return false;
            batch = oldest_claimable();
            if (batch == nullptr)
                return false;
            ++batch->workers;      // pins the batch against collection
            ++helpers_in_flight_;  // pins the window against destruction
        }
        const std::int64_t item =
            batch->claim.fetch_add(1, std::memory_order_relaxed);
        if (item < batch->n) {
            try {
                batch->body(item, tid);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_);
                if (!batch->first_exception)
                    batch->first_exception = std::current_exception();
                batch->claim.store(batch->n, std::memory_order_relaxed);
            }
        }
        bool became_quiescent = false;
        {
            std::lock_guard<std::mutex> lock(mutex_);
            --batch->workers;
            became_quiescent = quiescent(*batch);
            if (--helpers_in_flight_ == 0 && shutdown_) {
                // Notify under the lock: the waiting destructor may destroy done_ as soon as
                // the mutex is released.
                done_.notify_all();
                return true;
            }
        }
        // `batch` may already be collected and destroyed.
        if (became_quiescent)
            done_.notify_all();
        return true;
    }

    std::size_t window() const { return window_; }

private:
    struct Batch {
        Batch(std::int64_t items, Body item_body, Payload item_payload)
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
    static bool quiescent(const Batch& batch) {
        return batch.claim.load(std::memory_order_relaxed) >= batch.n &&
               batch.workers == 0;
    }

    // Called with mutex_ held: the oldest batch with an unclaimed item, or nullptr.
    Batch* oldest_claimable() {
        for (auto& batch : batches_) {
            if (batch->claim.load(std::memory_order_relaxed) < batch->n)
                return batch.get();
        }
        return nullptr;
    }

    // Claims and runs items until the batch has none left.
    void drain(Batch& batch, int tid) {
        for (;;) {
            const std::int64_t item =
                batch.claim.fetch_add(1, std::memory_order_relaxed);
            if (item >= batch.n)
                return;
            try {
                batch.body(item, tid);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex_);
                if (!batch.first_exception)
                    batch.first_exception = std::current_exception();
                batch.claim.store(batch.n, std::memory_order_relaxed);
                return;
            }
        }
    }

    void worker_loop(int tid) {
        name_current_thread("fa-map", tid);
        for (;;) {
            Batch* batch = nullptr;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                for (;;) {
                    if (shutdown_)
                        return;
                    batch = oldest_claimable();
                    if (batch != nullptr) {
                        ++batch->workers;  // pins it against collection
                        break;
                    }
                    work_.wait(lock);
                }
            }
            drain(*batch, tid);
            bool became_quiescent = false;
            {
                std::lock_guard<std::mutex> lock(mutex_);
                --batch->workers;
                became_quiescent = quiescent(*batch);
            }
            // `batch` may already be collected and destroyed.
            if (became_quiescent)
                done_.notify_all();
        }
    }

    mutable std::mutex mutex_;      // guards batches_, workers, first_exception
    std::condition_variable work_;  // a batch was submitted / teardown began
    std::condition_variable done_;  // a batch became quiescent
    std::deque<std::unique_ptr<Batch>> batches_;  // oldest at the front
    std::vector<std::thread> workers_;
    int n_threads_;
    std::size_t window_;
    bool shutdown_ = false;
    // run_one callers inside an item; guarded by mutex_. The destructor waits for zero.
    int helpers_in_flight_ = 0;
};

}}}  // namespace fa::cpu::threading
