# the compare print leaves unselected main differences out of the table

    Code
      print(make_compare_verdicts(main_selected = FALSE))
    Output
      Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):
      presence: log BF > 2.30; absence: log BF < -2.30
      
        presence 1 | undecided 2 | absence 0   (3 indicators)
        Main-effect differences are not under selection (main_difference_selection = FALSE).
      
            parameter  pip log_bf   verdict fragile
       A-B (pairwise) 0.95   2.94  presence   FALSE
       A-C (pairwise) 0.10  -2.20 undecided    TRUE
       B-C (pairwise) 0.52   0.08 undecided   FALSE
      
      1 verdict is Monte-Carlo fragile: a verdict boundary lies within two standard
      errors of the evidence, so the verdict could change on a rerun. Consider a
      longer run.
      
      The fragility flag is not validated for difference indicators: its
      operating point was established on single-network edge indicators only.
      Read it as an indication that a verdict sits near a boundary, not as a
      calibrated error rate.
      
      Difference verdicts are scale-contingent: group differences are priced
      on the association scale through difference_scale, and the calibration
      of that default is under study, so a verdict close to a decision
      threshold can move with the scale.

# the compare print tabulates selected main differences as rows

    Code
      print(make_compare_verdicts(main_selected = TRUE))
    Output
      Edge verdicts at an inclusion Bayes factor of 10 (and 0.1 for absence):
      presence: log BF > 2.30; absence: log BF < -2.30
      
        presence 1 | undecided 4 | absence 1   (6 indicators)
      
            parameter  pip log_bf   verdict fragile
             A (main) 0.05  -2.94   absence   FALSE
       A-B (pairwise) 0.95   2.94  presence   FALSE
       A-C (pairwise) 0.10  -2.20 undecided    TRUE
             B (main) 0.90   2.20 undecided   FALSE
       B-C (pairwise) 0.52   0.08 undecided   FALSE
             C (main) 0.50   0.00 undecided   FALSE
      
      1 verdict is Monte-Carlo fragile: a verdict boundary lies within two standard
      errors of the evidence, so the verdict could change on a rerun. Consider a
      longer run.
      
      The fragility flag is not validated for difference indicators: its
      operating point was established on single-network edge indicators only.
      Read it as an indication that a verdict sits near a boundary, not as a
      calibrated error rate.
      
      Difference verdicts are scale-contingent: group differences are priced
      on the association scale through difference_scale, and the calibration
      of that default is under study, so a verdict close to a decision
      threshold can move with the scale.

