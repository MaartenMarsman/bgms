#ifndef PROGRESS_MANAGER_H
#define PROGRESS_MANAGER_H

#include <Rcpp.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <numeric>

using Clock = std::chrono::steady_clock;

/**
 * Interrupt-check callback passed to R_ToplevelExec.
 *
 * Origin: rcpp_progress interrupt helpers,
 * https://github.com/kforner/rcpp_progress/blob/d851ac62fd0314239e852392de7face5fa4bf48e/inst/include/interrupts.hpp#L24-L31
 */
static void chkIntFn(void *dummy) {
	R_CheckUserInterrupt();
}

/**
 * Check for a pending user interrupt. Runs chkIntFn in a top-level context
 * so the interrupt cannot longjmp out of the caller's context.
 *
 * @return true if the user has requested an interrupt
 */
inline bool checkInterrupt() {
	return (R_ToplevelExec(chkIntFn, NULL) == FALSE);
}

/**
 * @brief Multi-chain progress bar manager for MCMC computations
 *
 * This class provides a thread-safe progress bar that works in both RStudio
 * console and terminal environments. It supports Unicode theming with colored
 * progress indicators and proper cursor positioning.
 *
 * Thread contract: the manager must be constructed on the R main thread. The
 * per-chain counters and the exit flag are atomics, so update() and
 * shouldExit() may be called from any thread. All R API interaction (the
 * interrupt check, console output, and the R callback) happens only on the
 * construction thread: update() drives it from whichever chains that thread
 * executes and skips the display on worker threads, so a worker thread never
 * touches the R interpreter.
 *
 * Key features:
 * - Multi-chain progress tracking with atomic counters
 * - RStudio vs terminal environment detection and adaptation
 * - Unicode and classic theming options
 * - ANSI color support with proper visual length calculations
 * - Console width adaptation and change detection
 * - User interrupt checking, confined to the R main thread
 * - Optional R callback for external progress reporting (e.g., JASP),
 *   invoked as callback(completed, total)
 */
class ProgressManager {

public:

    /**
     * Construct on the R main thread.
     *
     * @param nChains_          Number of parallel chains
     * @param nIter_            Total iterations per chain
     * @param nWarmup_          Warmup iterations per chain
     * @param printEvery_       Print frequency in iterations
     * @param progress_type     Bar style (0 = none, 1 = total, 2 = per-chain)
     * @param useUnicode_       Use Unicode theme instead of ASCII
     * @param progress_callback Optional R function called as callback(completed, total)
     */
    ProgressManager(int nChains_, int nIter_, int nWarmup_, int printEvery_ = 10, int progress_type = 2, bool useUnicode_ = true, SEXP progress_callback = R_NilValue);

    /**
     * Record one completed iteration for a chain. Callable from any thread;
     * printing and interrupt checks run only when called on the R main thread.
     */
    void update(size_t chainId);

    /** Print the final progress state and release the progress lines. */
    void finish();

    /** @return true if a user interrupt was detected and chains should stop. */
    bool shouldExit() const;

private:

    // Runs the R-main-thread work (interrupt check, throttled print, callback).
    // A no-op when called off the construction thread.
    void poll();

    void checkConsoleWidthChange();
    size_t getConsoleWidth() const;
    std::string formatProgressBar(size_t chainId, size_t current, size_t total, double fraction, bool isTotal = false) const;
    std::string formatTimeInfo(double elapsed, double eta) const;
    std::string formatDuration(double seconds) const;
    void setupTheme();

    bool isWarmupPhase() const {
        for (const auto& c : progress)
            if (c.load(std::memory_order_relaxed) < nWarmup)
                return true;
        return false;
    }
    bool isWarmupPhase(const size_t chain_id) const {
      return progress[chain_id].load(std::memory_order_relaxed) < nWarmup;
    }

    size_t totalProgress() const {
        size_t done = 0;
        for (const auto& c : progress)
            done += c.load(std::memory_order_relaxed);
        return done;
    }

    void print();

    void update_prefixes(size_t  width);

    void maybePadToLength(std::string& content) const;

    // set by constructor
    size_t nChains;                    // Number of parallel chains
    size_t nIter;                      // Total Iterations per chain
    size_t nWarmup;                    // Warmup iterations per chain
    size_t printEvery;                 // Print frequency
    size_t progress_type = 2;          // Progress bar style type (0 = "none", 1 = "total", 2 = "per-chain")
    bool useUnicode = true;            // Use Unicode vs ASCII theme
    std::vector<std::atomic<size_t>> progress; // Per-chain progress counters
    std::thread::id main_thread_id;    // Thread the manager was constructed on (the R main thread)
    size_t main_thread_updates_ = 0;   // update() calls seen on the main thread; throttles poll() (main-thread only, no atomic)

    // internal config parameters/ data
    size_t no_spaces_for_total;     // Spacing for total line alignment
    size_t lastPrintedLines = 0;    // Lines printed in last update
    size_t lastPrintedChars = 0;    // Characters printed in last update (RStudio)
    size_t consoleWidth = 80;       // Current console width
    size_t lineWidth = 80;          // Target line width for content
    int prevConsoleWidth = -1;      // Previous console width for change detection

    // Environment and state flags
    bool isRStudio = false;              ///< Whether running in RStudio console
    std::atomic<bool> needsToExit{false}; ///< User interrupt flag (set on the main thread, read by all chains)
    bool widthChanged = false;           ///< Console width changed flag

    // Visual configuration
    size_t barWidth = 40;              // Progress bar width in characters

    // Theme tokens
    std::string lhsToken;           // Left bracket/delimiter
    std::string rhsToken;           // Right bracket/delimiter
    std::string filledToken;        // Filled progress character
    std::string emptyToken;         // Empty progress character
    std::string partialTokenMore;   // Partial progress (>50%)
    std::string partialTokenLess;   // Partial progress (<50%)
    std::string chain_prefix;       // Chain label prefix
    std::string total_prefix;       // Total label prefix
    std::string total_padding;      // Padding for total line alignment

    // Timing
    Clock::time_point start;                   // Start time
    std::chrono::time_point<Clock> lastPrint;  // Last print time

    // Thread synchronization
    std::mutex printMutex;          // Mutex for thread-safe printing

    // R callback (called as callback(completed, total) at throttled intervals)
    Rcpp::Nullable<Rcpp::Function> callback;
};

#endif // PROGRESS_MANAGER_H