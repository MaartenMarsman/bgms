#pragma once

#include <atomic>
#include <chrono>
#include <exception>
#include <functional>
#include <stdexcept>
#include <thread>

#if defined(_WIN32)
#include <windows.h>
#include <process.h>
#else
#include <pthread.h>
#endif

#include "utils/progress_manager.h"

// Runs blocking parallel work (e.g. RcppParallel::parallelFor over chains) on
// a helper thread while the calling thread -- the R main thread -- polls the
// progress manager for user interrupts, progress display, and the R callback.
//
// The R API is single-threaded: interrupt checks, console output, and R
// closures are only safe on the R main thread. parallelFor blocks its caller
// and hands sub-ranges to arbitrary scheduler threads, so no code inside it
// may touch the R API. Moving the parallelFor itself onto a helper thread
// frees the main thread to do all R interaction in the poll loop below;
// workers communicate through the manager's atomic counters and exit flag.
//
// The helper thread gets an explicit 8 MB stack: the scheduler runs part of
// the range inline on it, so it needs at least as much stack as any worker.

namespace bgms_threads {

inline constexpr std::size_t helper_stack_size = 8 * 1024 * 1024;

namespace detail {

struct ThreadPayload {
    std::function<void()> work;
};

#if defined(_WIN32)
inline unsigned __stdcall thread_entry(void* arg) {
    ThreadPayload* payload = static_cast<ThreadPayload*>(arg);
    payload->work();
    delete payload;
    return 0;
}
#else
inline void* thread_entry(void* arg) {
    ThreadPayload* payload = static_cast<ThreadPayload*>(arg);
    payload->work();
    delete payload;
    return nullptr;
}
#endif

// Launches `work` on a joinable native thread with an explicit stack size and
// returns a join function.
inline std::function<void()> launch_with_stack(std::function<void()> work) {
    ThreadPayload* payload = new ThreadPayload{std::move(work)};

#if defined(_WIN32)
    uintptr_t handle = _beginthreadex(
        nullptr, static_cast<unsigned>(helper_stack_size), thread_entry,
        payload, 0, nullptr);
    if (handle == 0) {
        delete payload;
        throw std::runtime_error("Failed to launch the chain helper thread.");
    }
    HANDLE h = reinterpret_cast<HANDLE>(handle);
    return [h]() {
        WaitForSingleObject(h, INFINITE);
        CloseHandle(h);
    };
#else
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setstacksize(&attr, helper_stack_size);
    pthread_t handle;
    int rc = pthread_create(&handle, &attr, thread_entry, payload);
    pthread_attr_destroy(&attr);
    if (rc != 0) {
        delete payload;
        throw std::runtime_error("Failed to launch the chain helper thread.");
    }
    return [handle]() { pthread_join(handle, nullptr); };
#endif
}

} // namespace detail

// Runs `work` on the helper thread and polls `pm` from the calling thread
// until the work completes. Exceptions from the helper are rethrown here.
template <typename Work>
inline void run_with_main_thread_progress(ProgressManager& pm, Work&& work) {
    std::atomic<bool> finished{false};
    std::exception_ptr error;

    auto join = detail::launch_with_stack([&]() {
        try {
            work();
        } catch (...) {
            error = std::current_exception();
        }
        finished.store(true, std::memory_order_release);
    });

    while (!finished.load(std::memory_order_acquire)) {
        pm.poll();
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }
    join();

    if (error) std::rethrow_exception(error);
}

} // namespace bgms_threads
