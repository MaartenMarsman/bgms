// -----------------------------------------------------------------------------
// warmup_schedule_test_interface.cpp
//
// R-facing test entry for WarmupSchedule. Exposes the stage boundaries and the
// selection gate for a given (warmup, edge_selection, learn_sd,
// select_during_warmup) configuration, so the per-sampler warmup staging can be
// tested without running a chain.
// -----------------------------------------------------------------------------

#include <Rcpp.h>

#include "mcmc/execution/warmup_schedule.h"

// -----------------------------------------------------------------------------
// test_warmup_schedule:
//   Construct a WarmupSchedule and return its stage boundaries plus the
//   selection_enabled() gate evaluated at each of the supplied iterations.
// -----------------------------------------------------------------------------

// [[Rcpp::export(name = "test_warmup_schedule")]]
Rcpp::List test_warmup_schedule(
    int warmup,
    bool edge_selection,
    bool learn_sd,
    bool select_during_warmup,
    Rcpp::IntegerVector probe_iterations
) {
    WarmupSchedule schedule(warmup, edge_selection, learn_sd, select_during_warmup);

    Rcpp::LogicalVector selection_enabled(probe_iterations.size());
    for (int k = 0; k < probe_iterations.size(); ++k) {
        selection_enabled[k] = schedule.selection_enabled(probe_iterations[k]);
    }

    return Rcpp::List::create(
        Rcpp::Named("stage1_end") = schedule.stage1_end,
        Rcpp::Named("stage3a_start") = schedule.stage3a_start,
        Rcpp::Named("stage3b_start") = schedule.stage3b_start,
        Rcpp::Named("stage3c_start") = schedule.stage3c_start,
        Rcpp::Named("total_warmup") = schedule.total_warmup,
        Rcpp::Named("stage3b_skipped") = schedule.stage3b_skipped,
        Rcpp::Named("selection_enabled") = selection_enabled
    );
}
