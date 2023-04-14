test_that("inclusion probabilities correlate with posterior mode", {
  data("Wenchuan", package = "bgms")
  fit <- bgm(x = Wenchuan, iter = 1e2, burnin = 10)

  posterior_modes      <- fit$interactions[lower.tri(fit$interactions)]
  posterior_incl_probs <- fit$gamma[lower.tri(fit$gamma)]

  testthat::expect_gte(cor(abs(posterior_modes), posterior_incl_probs, method = "spearman"), .9)

})
