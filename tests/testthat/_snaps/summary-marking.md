# an empty group-by-category cell marks its difference rows

    Code
      print(summary(fit))
    Output
      Posterior summaries from Bayesian grouped MRF estimation (bgmCompare):
      
      groups: 1 = 1 (n = 80), 2 = 2 (n = 80)
      
      Category thresholds:
      parameter mean mcse sd n_eff Rhat
      1 A (1)
      2 A (2)
      3 A (3)
      4 B (1)
      5 B (2)
      6 C (1)
      ... (use `summary(fit)$main` to see full output)
      
      Pairwise interactions:
      parameter mean mcse sd n_eff Rhat
      1 A-B
      2 A-C
      3 B-C
      
      Inclusion probabilities:
      parameter mean mcse sd n_eff Rhat n0->1 n1->0
      A (main)
      A-B (pairwise)
      A-C (pairwise)
      B (main)
      B-C (pairwise)
      C (main)
      Note: blank mcse/n_eff/Rhat cells mark indicators whose inclusion draws never
      varied: the per-iteration evidence is so one-sided that the draws round to
      exactly 0 or 1, and no variation is left to estimate precision from. An
      indicator that is almost as certain still shows numbers; the difference is
      rounding, not evidence.
      Main-effect difference rows are blank when main_difference_selection = FALSE
      left their indicators unsampled.
      All computed values remain in `summary(fit)$indicator`.
      
      Group differences (main effects):
      parameter mean mcse sd n_eff share_incl Rhat
      * A (diff1; 1)
      * A (diff1; 2)
      * A (diff1; 3)
      B (diff1; 1)
      B (diff1; 2)
      C (diff1; 1)
      ... (use `summary(fit)$main_diff` to see full output)
      * a group lacks observations in this category or in the reference category; the estimate reflects the prior, not the data
      Note: NA values are suppressed in the print table. They occur for differences
      that were never selected, so the composite ESS and share are undefined;
      `summary(fit)$main_diff` still contains the NA values.
      
      Group differences (pairwise effects):
      parameter mean mcse sd n_eff share_incl Rhat
      A-B (diff1)
      A-C (diff1)
      B-C (diff1)
      Note: NA values are suppressed in the print table. They occur for differences
      that were never selected, so the composite ESS and share are undefined;
      `summary(fit)$pairwise_diff` still contains the NA values.
      
      Use `summary(fit)$<component>` to access full results.
      See the `easybgm` package for other summary and plotting tools.

