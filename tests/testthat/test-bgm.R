test_that("inclusion probabilities correlate with posterior mode", {
  data("Wenchuan", package = "bgms")
  fit <- bgm(x = Wenchuan, iter = 1e2, burnin = 10)

  posterior_modes = extract_pairwise_interactions(fit)
  posterior_incl_probs = extract_posterior_inclusion_probabilities(fit)

  posterior_modes = posterior_modes[lower.tri(posterior_modes)]
  posterior_incl_probs = posterior_incl_probs[lower.tri(posterior_incl_probs)]

  testthat::expect_gte(cor(abs(posterior_modes), posterior_incl_probs, method = "spearman"), .9)

})
