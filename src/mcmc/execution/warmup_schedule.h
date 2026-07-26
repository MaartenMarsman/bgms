#pragma once

#include <vector>
#include <algorithm>
#include <cmath>
#include <optional>


/**
 * WarmupSchedule - Dynamic warmup schedule with adaptive windows
 *
 * For edge_selection = FALSE:
 *   Stage 1 (init), Stage 2 (doubling windows), Stage 3a (terminal).
 *   total_warmup = user-specified warmup.
 *
 * For edge_selection = TRUE with learn_sd (NUTS):
 *   User warmup is split: 85% for Stage 1-3a, 10% for Stage 3b, 5% for Stage 3c.
 *   - Stage 3b: proposal SD tuning for edge-selection Metropolis moves
 *   - Stage 3c: step size re-adaptation with edge selection active
 *   If Stage 3b would get < 20 iterations, it's skipped (uses default proposal SD).
 *
 * For edge_selection = TRUE with select_during_warmup (Gibbs):
 *   The sampler has no adaptation, so warmup is split into a short full-model
 *   settle window (the first gibbs_settle_fraction of warmup, selection off)
 *   followed by Stage 3c covering the remainder (selection active). Both
 *   windows rescale with the user's warmup budget.
 *
 * For edge_selection = TRUE otherwise (adaptive-Metropolis):
 *   Stage 3c is empty; selection activates at the first sampling iteration.
 *
 * Warning types:
 *   0 = none, 1 = warmup extremely short (< 50),
 *   2 = core stages using proportional fallback,
 *   3 = limited proposal SD tuning (edge_selection && warmup < 300),
 *   4 = Stage 3b skipped (would have < 20 iterations).
 */
struct WarmupSchedule {
  int stage1_end;                 ///< Stage-1 [0 ... stage1_end-1]
  std::vector<int> window_ends;   ///< Stage-2 windows (last index of each)
  int stage3a_start;              ///< First iter in Stage-3a
  int stage3b_start;              ///< First iter in Stage-3b (== stage3c_start if skipped)
  int stage3c_start;              ///< First iter in Stage-3c (== total_warmup if skipped)
  int total_warmup;               ///< User-specified warmup (== first retained iteration)
  bool learn_proposal_sd;         ///< Whether to run the proposal-SD tuner
  bool enable_selection;          ///< Allow edge-indicator moves
  int warning_type;               ///< Warning code (see above)
  bool stage3b_skipped;           ///< True if 3b was skipped due to insufficient budget

  WarmupSchedule(int warmup,
                 bool enable_sel,
                 bool learn_sd,
                 bool select_during_warmup = false)
    : stage1_end(0)
    , window_ends()
    , stage3a_start(0)
    , stage3b_start(0)
    , stage3c_start(0)
    , total_warmup(warmup)
    , learn_proposal_sd(learn_sd)
    , enable_selection(enable_sel)
    , warning_type(0)
    , stage3b_skipped(false)
  {
    // ===== Step 1: Determine budget allocation =====
    int warmup_core;    // Budget for Stage 1-3a (mass matrix + step size)
    int stage3b_budget = 0;
    int selection_warmup_start = -1;  // >= 0: Stage 3c starts here (Gibbs)

    if (enable_sel && learn_sd) {
      // For edge selection models: split warmup as 85%/10%/5%
      warmup_core    = static_cast<int>(0.85 * warmup);
      stage3b_budget = static_cast<int>(0.10 * warmup);
      // Stage 3c gets remainder: warmup - warmup_core - stage3b_budget (~5%)

      // Check if Stage 3b has enough iterations to be meaningful
      if (stage3b_budget < 20) {
        // Skip Stage 3b entirely - will use default proposal SD
        stage3b_skipped = true;
        warning_type = 4;
        warmup_core = warmup;  // Give all warmup to core stages
        stage3b_budget = 0;
      } else if (warmup < 300) {
        // Marginal but runs - warn about limited tuning
        warning_type = 3;
      }
    } else if (enable_sel && select_during_warmup) {
      // Tuning-free sampler (Gibbs): a short full-model settle window, then
      // selection-active warmup (Stage 3c) for the remainder of the budget.
      warmup_core = warmup;
      selection_warmup_start = static_cast<int>(gibbs_settle_fraction * warmup);
    } else {
      // No edge selection: all warmup goes to core stages 1-3a
      warmup_core = warmup;
    }

    // ===== Step 2: Set up core stages (1, 2, 3a) =====
    constexpr int default_init_buffer = 75;   // Stage-1: initial fast adaptation
    constexpr int default_term_buffer = 50;   // Stage-3a: final fast adaptation
    constexpr int default_base_window = 25;   // Stage-2: initial window size

    int init_buffer, term_buffer, base_window;

    if (warmup_core < 20) {
      // Too short for any meaningful adaptation
      if (warning_type == 0) warning_type = 1;  // Don't overwrite more specific warnings
      init_buffer = warmup_core;
      term_buffer = 0;
      base_window = 0;
    } else if (default_init_buffer + default_base_window + default_term_buffer > warmup_core) {
      // Not enough room for fixed buffers; fall back to proportional (15%/75%/10%)
      if (warning_type == 0) warning_type = 2;
      init_buffer = static_cast<int>(0.15 * warmup_core);
      term_buffer = static_cast<int>(0.10 * warmup_core);
      base_window = warmup_core - init_buffer - term_buffer;
    } else {
      // Standard case: use fixed buffers
      init_buffer = default_init_buffer;
      term_buffer = default_term_buffer;
      base_window = default_base_window;
    }

    // Additional warning for extremely short warmup with edge selection
    if (enable_sel && warmup < 50 && warning_type != 1) {
      warning_type = 1;  // Override to most severe
    }

    /* ---------- Stage-1 ---------- */
    stage1_end = init_buffer;

    /* ---------- Stage-3a start ---------- */
    stage3a_start = warmup_core - term_buffer;

    /* ---------- Stage-2: build doubling windows ---------- *
     * Matches Stan's windowed_adaptation::compute_next_window: when the
     * window AFTER the next would overshoot stage3a_start, STRETCH the
     * current next window to absorb the remaining stage-2 budget instead
     * of emitting a small trailing window. Each window_ends entry triggers
     * a mass-matrix update + step-size reinit, so an extra small final
     * window disrupts dual averaging.
     */
    if (base_window > 0 && stage3a_start > stage1_end) {
      int cur = stage1_end;
      int wsize = base_window;
      int next_end = cur + wsize;
      while (next_end <= stage3a_start) {
        window_ends.push_back(next_end);
        cur = next_end;
        if (next_end == stage3a_start) break;
        wsize *= 2;
        next_end = cur + wsize;
        int after_next = next_end + 2 * wsize;
        if (after_next > stage3a_start) {
          next_end = stage3a_start;
        }
      }
    }

    /* ---------- Stage-3b and 3c boundaries ---------- */
    stage3b_start = warmup_core;
    stage3c_start = (selection_warmup_start >= 0)
      ? selection_warmup_start              // Gibbs: settle, then selection
      : warmup_core + stage3b_budget;
  }

  /// Stage query helpers
  bool in_stage1 (int i) const { return i <  stage1_end; }
  bool in_stage2 (int i) const { return i >= stage1_end && i < stage3a_start; }
  bool in_stage3a(int i) const { return i >= stage3a_start && i < stage3b_start; }
  bool in_stage3b(int i) const { return !stage3b_skipped && i >= stage3b_start && i < stage3c_start; }
  bool in_stage3c(int i) const { return enable_selection && !stage3b_skipped && i >= stage3c_start && i < total_warmup; }
  bool sampling (int i) const { return i >= total_warmup; }

  bool has_warning() const { return warning_type > 0; }
  bool warmup_extremely_short() const { return warning_type == 1; }
  bool using_proportional_fallback() const { return warning_type == 2; }
  bool limited_proposal_tuning() const { return warning_type == 3; }
  bool proposal_tuning_skipped() const { return warning_type == 4 || stage3b_skipped; }

  /// Whether indicator moves are enabled (Stage 3c and sampling)
  bool selection_enabled(int i) const {
    return enable_selection && (in_stage3c(i) || sampling(i));
  }

  /// Whether to adapt proposal_sd (Stage-3b only, if not skipped)
  bool adapt_proposal_sd(int i) const {
    return learn_proposal_sd && !stage3b_skipped && in_stage3b(i);
  }

  /// Robbins-Monro decay rate for proposal-SD adaptation. Single source of
  /// truth; every model's tune_proposal_sd consults this.
  static constexpr double proposal_sd_rm_decay = 0.75;

  /// Fraction of the Gibbs warmup spent settling the full model (all edges
  /// included, selection off) before Stage 3c activates the edge moves.
  static constexpr double gibbs_settle_fraction = 0.15;

  /// Robbins-Monro weight for proposal-SD adaptation at the given iteration.
  ///
  /// Returns the RM weight (1-indexed since stage 3b began, decay
  /// `proposal_sd_rm_decay`) iff `adapt_proposal_sd(iter)` is true. Returns
  /// nullopt otherwise. Use this in every `tune_proposal_sd` rather than
  /// computing the weight inline — keeps the policy (which iterations adapt,
  /// what the weight schedule looks like) in one place and makes stage-3b
  /// adaptation sampler-agnostic by construction.
  std::optional<double> rm_weight_for_proposal_sd(int iter) const {
    if (!adapt_proposal_sd(iter)) return std::nullopt;
    const double t = static_cast<double>(iter - stage3b_start + 1);
    return std::pow(t, -proposal_sd_rm_decay);
  }

  /// Current Stage-2 window index (-1 outside Stage-2)
  int current_window(int i) const {
    for (size_t k = 0; k < window_ends.size(); ++k)
      if (i < window_ends[k]) return static_cast<int>(k);
    return -1;
  }
};
