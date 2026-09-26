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
#ifndef FA_THREADING_BATCH_WINDOW_H
#define FA_THREADING_BATCH_WINDOW_H

#include "threading/thread_name.h"

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <exception>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

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
    // `tid` is only passed to the body and need not lie in [0, thread_count()). Callers must
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
    int thread_count() const { return n_threads_; }

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

#endif  // FA_THREADING_BATCH_WINDOW_H
