#include "utils/progress_manager.h"

ProgressManager::ProgressManager(int nChains_, int nIter_, int nWarmup_, int printEvery_, int progress_type_, bool useUnicode_, SEXP progress_callback)
  : nChains(nChains_), nIter(nIter_ + nWarmup_), nWarmup(nWarmup_), printEvery(printEvery_),
    progress_type(progress_type_), useUnicode(useUnicode_), progress(nChains_), callback(progress_callback) {

  // When a callback is provided, suppress the built-in progress display
  if (callback.isNotNull()) progress_type = 0;

  // R API interaction (interrupt check, printing, callback) is confined to
  // the thread that constructs the manager: the R main thread.
  main_thread_id = std::this_thread::get_id();

  for (size_t i = 0; i < nChains; i++) progress[i] = 0;
  start = Clock::now();
  lastPrint = Clock::now();

  // Check if we're in RStudio
  Rcpp::Environment base("package:base");
  Rcpp::Function getOption = base["getOption"];
  Rcpp::Function Sysgetenv = base["Sys.getenv"];
  SEXP s = Sysgetenv("RSTUDIO");
  isRStudio = Rcpp::as<std::string>(s) == "1";

  no_spaces_for_total = 3 + static_cast<size_t>(std::log10(nChains));
  if (progress_type == 1) no_spaces_for_total = 1; // no need to align, so one space is fine
  total_padding = std::string(no_spaces_for_total, ' ');

  if (isRStudio) {
    consoleWidth = getConsoleWidth();
    lineWidth = std::max(10, std::min(static_cast<int>(consoleWidth) - 25, 70));
  } else {
    // For terminal, use default console width
    consoleWidth = 80;
    lineWidth = 70;
  }

  // cleverly determine barwidth so no line wrapping occurs
  if (lineWidth <= 5) {
    // TODO: we don't want to print anything in this case
    barWidth = 0;
  } else if (lineWidth < 20) {
    barWidth = lineWidth - 10;
  } else if (lineWidth < 40) {
    barWidth = lineWidth - 15;
  } else { // lineWidth > 40
    barWidth = lineWidth > 70 ? 40 : lineWidth - 30;
  }

  if (isRStudio) {
    barWidth = barWidth > 30 ? barWidth - 20 : 10; // minimum bar width of 10 for RStudio

  }

  // Set up theme
  setupTheme();
  update_prefixes(consoleWidth);
}

void ProgressManager::update(size_t chainId) {
  const size_t count = progress[chainId].fetch_add(1, std::memory_order_relaxed) + 1;

  // Worker threads only count. The R API (interrupt check, printing, the
  // callback) is only legal on the R main thread: in a serial run every
  // update lands there and drives the display below; in a parallel run the
  // launcher's poll() loop does instead.
  if (std::this_thread::get_id() != main_thread_id) return;

  if (count % printEvery == 0) poll();
}

void ProgressManager::poll() {
  // Main-thread only: calling this from any other thread is a no-op.
  if (std::this_thread::get_id() != main_thread_id) return;
  if (needsToExit.load(std::memory_order_relaxed)) return;

  // Check for user interrupts
  if (checkInterrupt()) {
    needsToExit.store(true, std::memory_order_relaxed);
    if (progress_type != 0) {
      // This should be immediately-ish visible to the user
      Rcpp::Rcout << "\nUser interrupt detected. Exiting gracefully. It may take a few seconds before all chains are terminated.\n";
    }
    return;
  }

  auto now = Clock::now();
  std::chrono::duration<double> sinceLast = now - lastPrint;

  bool has_output = (progress_type != 0) || callback.isNotNull();

  // Throttle printing to avoid spamming
  if (has_output && sinceLast.count() >= 0.5) {
    if (progress_type != 0) {
      print();
    }
    if (callback.isNotNull()) {
      size_t done = totalProgress();
      size_t totalWork = nChains * nIter;
      Rcpp::Function(callback.get())(done, totalWork);
    }
    lastPrint = now;
  }
}

void ProgressManager::finish() {

  if (progress_type == 0 && callback.isNull()) return;

  if (needsToExit) {
    if (progress_type != 0)
      Rcpp::Rcout << "All chains terminated.\n";
    return;
  }

  // Mark all chains as complete and print one final time
  for (size_t i = 0; i < nChains; i++)
    progress[i] = nIter;

  if (progress_type != 0)
    print();
  if (callback.isNotNull()) {
    size_t totalWork = nChains * nIter;
    Rcpp::Function(callback.get())(totalWork, totalWork);
  }

}

bool ProgressManager::shouldExit() const {
  return needsToExit.load(std::memory_order_relaxed);
}

void ProgressManager::checkConsoleWidthChange() {
  if (!isRStudio) return;

  size_t currentWidth = getConsoleWidth();

  if (prevConsoleWidth == -1) {
    // First time, just store the current width
    prevConsoleWidth = consoleWidth;
    widthChanged = false;
    return;
  }

  if (currentWidth != consoleWidth && currentWidth > 0) {
    // Width has changed
    prevConsoleWidth = consoleWidth;
    consoleWidth = currentWidth;
    widthChanged = true;
  } else {
    widthChanged = false;
  }
}

size_t ProgressManager::getConsoleWidth() const {
  Rcpp::Environment base("package:base");
  Rcpp::Function getOption = base["getOption"];
  SEXP s = getOption("width", 0);
  size_t width = std::max(0, Rcpp::as<int>(s));
  return width + 3; // note: + 3 is not entirely accurate, in reality this is scales with the actual width
}

std::string ProgressManager::formatProgressBar(size_t chainId, size_t current, size_t total, double fraction, bool isTotal) const {
    std::ostringstream builder;

    double exactFilled = fraction * barWidth;
    size_t filled = std::min(static_cast<size_t>(exactFilled), barWidth);

    // Build progress bar with theme
    std::string progressBar = lhsToken;

    // Add filled tokens
    for (size_t i = 0; i < filled; i++)
      progressBar += filledToken;

    // Add partial token if needed
    if (filled < barWidth) {
      double partialAmount = exactFilled - filled;
      if (partialAmount > 0) {
        if (partialAmount > 0.5) {
          progressBar += partialTokenMore;
        } else {
          progressBar += partialTokenLess;
        }
        filled++; // Account for the partial token
      }
    }

    // Add empty tokens
    for (size_t i = filled; i < barWidth; i++)
      progressBar += emptyToken;

    progressBar += rhsToken;

    // store the current length of the progress bar without any additional text
    size_t currentWidth = progressBar.length();

    if (isTotal) {

        std::string warmupOrSampling = isWarmupPhase() ? "(Warmup)" : "(Sampling)";
        builder << total_prefix << total_padding << warmupOrSampling << ": " << progressBar << " " << current << "/" << total
                << " (" << std::fixed << std::setprecision(1) << fraction * 100 << "%)";

    } else {

        std::string warmupOrSampling = isWarmupPhase(chainId - 1) ? " (Warmup)" : " (Sampling)";
        builder << chain_prefix << " " << chainId << warmupOrSampling << ": " << progressBar << " " << current << "/" << total
                << " (" << std::fixed << std::setprecision(1) << fraction * 100 << "%)";
    }

    std::string output = builder.str();

    if (isRStudio && progress_type == 2) {

      currentWidth = output.length() + 2 + barWidth - currentWidth;
      // Pad each line to exactly lineWidth characters (before adding \n)
      if (currentWidth < lineWidth) {
        output += std::string(lineWidth - currentWidth, ' ');
      } else if (currentWidth > lineWidth) {
        output = output.substr(0, lineWidth);
      }
    }

    return output;
}

std::string ProgressManager::formatTimeInfo(double elapsed, double eta) const {
  std::ostringstream builder;
  builder << "Elapsed: " << formatDuration(elapsed) << " | ETA: " << formatDuration(eta);
  return builder.str();
}

// Add this helper function to the class
std::string ProgressManager::formatDuration(double seconds) const {
  if (seconds < 0) {
    return "0s";
  }

  // Convert to different units
  if (seconds < 60) {
    // Less than 1 minute: show seconds
    return std::to_string(static_cast<int>(std::round(seconds))) + "s";
  }
  else if (seconds < 3600) {
    // Less than 1 hour: show minutes and seconds
    size_t mins = static_cast<size_t>(seconds / 60);
    size_t secs = static_cast<size_t>(seconds) % 60;
    if (secs == 0) {
      return std::to_string(mins) + "m";
    } else {
      return std::to_string(mins) + "m " + std::to_string(secs) + "s";
    }
  }
  else if (seconds < 86400) {
    // Less than 1 day: show hours and minutes
    size_t hours = static_cast<size_t>(seconds / 3600);
    size_t mins = static_cast<size_t>((seconds - hours * 3600) / 60);
    if (mins == 0) {
      return std::to_string(hours) + "h";
    } else {
      return std::to_string(hours) + "h " + std::to_string(mins) + "m";
    }
  }
  else {
    // 1 day or more: show days and hours
    size_t days = static_cast<size_t>(seconds / 86400);
    size_t hours = static_cast<size_t>((seconds - days * 86400) / 3600);
    if (hours == 0) {
      return std::to_string(days) + "d";
    } else {
      return std::to_string(days) + "d " + std::to_string(hours) + "h";
    }
  }
}

void ProgressManager::setupTheme() {
  // should be a struct of some kind...
  if (useUnicode) {
    // Unicode theme
    lhsToken         = "⦗";
    rhsToken         = "⦘";
    filledToken      = "\033[38;5;73m━\033[39m";  // Blue filled
    partialTokenMore = "\033[38;5;73m━\033[39m";  // Blue partial (> 0.5)
    partialTokenLess = "\033[37m╺\033[39m";       // Gray partial (< 0.5)
    emptyToken       = "\033[37m━\033[39m";       // Gray empty
  } else {
    // Classic theme
    lhsToken         = "[";
    rhsToken         = "]";
    filledToken      = "=";
    emptyToken       = " ";
    partialTokenMore = " ";
    partialTokenLess = " ";
  }
}

void ProgressManager::print() {
  std::lock_guard<std::mutex> lock(printMutex);

  auto now = Clock::now();
  double elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();

  size_t totalWork = nChains * nIter;
  size_t done = totalProgress();
  double fracTotal = double(done) / totalWork;
  // should actually be the eta of the slowest chain!
  double eta = (fracTotal > 0) ? elapsed / fracTotal - elapsed : 0.0;

  std::ostringstream out;
  // int totalChars = 0;

  // if this is not the first print, delete previous content
  if (progress_type == 2) {

    if (lastPrintedChars > 0) {

      if (isRStudio) {
        out << "\x1b[" << std::to_string((1 + lineWidth) * lastPrintedLines) << "D";
        // out << "\x1b[" << std::to_string(lastPrintedChars) << "D";
      } else {
        // Move cursor up to start of our content and clear everything
        for (size_t i = 0; i < lastPrintedLines; i++) {
          out << "\x1b[1A\x1b[2K"; // Move up one line and clear entire line
        }
      }
    }

    // Print progress for each chain
    for (size_t i = 0; i < nChains; i++) {
      const size_t chain_done = progress[i].load(std::memory_order_relaxed);
      double frac = double(chain_done) / nIter;
      std::string chainProgress = formatProgressBar(i + 1, chain_done, nIter, frac);
      out << chainProgress << "\n";
      // totalChars += chainProgress.length() + 1; // +1 for newline
    }

    // Print total progress if there is more than one chain
    if (nChains > 1) {
      std::string totalProgress = formatProgressBar(0, done, totalWork, fracTotal, true);
      out << totalProgress << "\n";
    }
    // totalChars += totalProgress.length() + 1; // +1 for newline

    // Print time info
    std::string timeInfo = formatTimeInfo(elapsed, eta);
    maybePadToLength(timeInfo);
    out << timeInfo << "\n";
    // totalChars += timeInfo.length() + 1; // +1 for newline

    // Track total lines printed (chains + total + time)

    lastPrintedLines = nChains + (nChains > 1 ? 2 : 1); // used in a generic terminal

    lastPrintedChars = 1;//totalChars;  // used by RStudio

  } else if (progress_type == 1) {

    // Print total progress
    std::string totalProgress = formatProgressBar(0, done, totalWork, fracTotal, true);

    // Print time info
    totalProgress += " " + formatTimeInfo(elapsed, eta);

    if (done < totalWork) {
      out << totalProgress << "\r";
    } else {
      out << totalProgress << "\n";
    }

    // we do not set lastPrintedChars or lastPrintedLines here since we always overwrite the same line

  }

  Rcpp::Rcout << out.str();

}

void ProgressManager::update_prefixes(size_t width) {
  if (width < 20) {
    chain_prefix = "C";
    total_prefix = "T";
  } else if (width < 30) {
    chain_prefix = "Ch";
    total_prefix = "Tot";
  } else {
    chain_prefix = "Chain";
    total_prefix = "Total";
  }
}

void ProgressManager::maybePadToLength(std::string& content) const {
  if (!isRStudio) return;

  // Pad each line to exactly lineWidth characters (before adding \n)
  if (content.length() < lineWidth) {
    content += std::string(lineWidth - content.length(), ' ');
  } else if (content.length() > lineWidth) {
    content = content.substr(0, lineWidth);
  }
}


// Example usage/ test with RcppParallel
// #include <RcppParallel.h>
// struct ChainWorker : public RcppParallel::Worker {
//   int nIter;
//   int delay; // no. milliseconds sleep per loop iteration
//   ProgressManager &pm;

//   ChainWorker(int nIter_, int delay_, ProgressManager &pm_)
//     : nIter(nIter_), delay(delay_), pm(pm_) {}

//   void operator()(std::size_t begin, std::size_t end) {

//     auto chainId = begin;

//     for (int i = 0; i < nIter; i++) {
//       // ---- Simulated work ----
//       std::this_thread::sleep_for(std::chrono::milliseconds(delay));

//       // ---- Update state ----
//       pm.update(chainId);
//       if (pm.shouldExit()) break;
//       // if (Progress::check_abort()) Rcpp::checkUserInterrupt();
//     }
//   }
// };


// // [[Rcpp::export]] // if uncommented, must move .cpp file to src/ for Rcpp to compile
// void runMCMC_parallel(int nChains = 4, int nIter = 100, int nWarmup = 100, int progress_type = 2, bool useUnicode = false,
//   int delay = 20) {

//   ProgressManager pm(nChains, nIter, nWarmup, 10, progress_type, useUnicode);
//   ChainWorker worker(nIter + nWarmup, delay, pm);

//   // Run each chain in parallel
//   RcppParallel::parallelFor(0, nChains, worker);

//   pm.finish();

//   if (pm.shouldExit()) {
//     Rcpp::Rcout << "\nComputation interrupted by user.\n";
//   } else {
//     Rcpp::Rcout << "\nAll chains finished!\n";
//   }
// }
