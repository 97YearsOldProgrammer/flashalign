// Ordered three-stage pipeline: produce() on a reader thread, process() on the caller and
// consume() on a writer thread, joined by bounded FIFO queues.
//
// process() returns std::nullopt while it holds a batch back (a stage with a window of
// batches in flight). At end of input, flush() is called until it returns std::nullopt and
// each value it returns goes to the writer in order.
#ifndef FA_THREADING_PIPELINE_H
#define FA_THREADING_PIPELINE_H

#include "threading/thread_name.h"

#include <condition_variable>
#include <cstddef>
#include <exception>
#include <mutex>
#include <optional>
#include <queue>
#include <thread>
#include <utility>

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

#endif  // FA_THREADING_PIPELINE_H
