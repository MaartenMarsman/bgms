// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// sample_omrf_gibbs
IntegerMatrix sample_omrf_gibbs(int no_states, int no_variables, IntegerVector no_categories, NumericMatrix interactions, NumericMatrix thresholds, int iter);
RcppExport SEXP _bgms_sample_omrf_gibbs(SEXP no_statesSEXP, SEXP no_variablesSEXP, SEXP no_categoriesSEXP, SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type no_states(no_statesSEXP);
    Rcpp::traits::input_parameter< int >::type no_variables(no_variablesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_omrf_gibbs(no_states, no_variables, no_categories, interactions, thresholds, iter));
    return rcpp_result_gen;
END_RCPP
}
// sample_bcomrf_gibbs
IntegerMatrix sample_bcomrf_gibbs(int no_states, int no_variables, IntegerVector no_categories, NumericMatrix interactions, NumericMatrix thresholds, StringVector variable_type, IntegerVector reference_category, int iter);
RcppExport SEXP _bgms_sample_bcomrf_gibbs(SEXP no_statesSEXP, SEXP no_variablesSEXP, SEXP no_categoriesSEXP, SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP variable_typeSEXP, SEXP reference_categorySEXP, SEXP iterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type no_states(no_statesSEXP);
    Rcpp::traits::input_parameter< int >::type no_variables(no_variablesSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< StringVector >::type variable_type(variable_typeSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type reference_category(reference_categorySEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    rcpp_result_gen = Rcpp::wrap(sample_bcomrf_gibbs(no_states, no_variables, no_categories, interactions, thresholds, variable_type, reference_category, iter));
    return rcpp_result_gen;
END_RCPP
}
// gibbs_sampler
List gibbs_sampler(arma::imat& observations, arma::imat& indicator, arma::mat& interactions, arma::mat& thresholds, const arma::ivec& num_categories, const double interaction_scale, arma::mat& proposal_sd, arma::mat& proposal_sd_blumecapel, const String& edge_prior, arma::mat& theta, const double beta_bernoulli_alpha, const double beta_bernoulli_beta, const double dirichlet_alpha, const double lambda, const arma::imat& Index, const int iter, const int burnin, arma::imat& num_obs_categories, arma::imat& sufficient_blume_capel, const double threshold_alpha, const double threshold_beta, const bool na_impute, const arma::imat& missing_index, const arma::uvec& is_ordinal_variable, const arma::ivec& reference_category, const bool save, const bool display_progress, bool edge_selection, bool mala);
RcppExport SEXP _bgms_gibbs_sampler(SEXP observationsSEXP, SEXP indicatorSEXP, SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP num_categoriesSEXP, SEXP interaction_scaleSEXP, SEXP proposal_sdSEXP, SEXP proposal_sd_blumecapelSEXP, SEXP edge_priorSEXP, SEXP thetaSEXP, SEXP beta_bernoulli_alphaSEXP, SEXP beta_bernoulli_betaSEXP, SEXP dirichlet_alphaSEXP, SEXP lambdaSEXP, SEXP IndexSEXP, SEXP iterSEXP, SEXP burninSEXP, SEXP num_obs_categoriesSEXP, SEXP sufficient_blume_capelSEXP, SEXP threshold_alphaSEXP, SEXP threshold_betaSEXP, SEXP na_imputeSEXP, SEXP missing_indexSEXP, SEXP is_ordinal_variableSEXP, SEXP reference_categorySEXP, SEXP saveSEXP, SEXP display_progressSEXP, SEXP edge_selectionSEXP, SEXP malaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type indicator(indicatorSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type num_categories(num_categoriesSEXP);
    Rcpp::traits::input_parameter< const double >::type interaction_scale(interaction_scaleSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type proposal_sd(proposal_sdSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type proposal_sd_blumecapel(proposal_sd_blumecapelSEXP);
    Rcpp::traits::input_parameter< const String& >::type edge_prior(edge_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_bernoulli_alpha(beta_bernoulli_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type beta_bernoulli_beta(beta_bernoulli_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type dirichlet_alpha(dirichlet_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type num_obs_categories(num_obs_categoriesSEXP);
    Rcpp::traits::input_parameter< arma::imat& >::type sufficient_blume_capel(sufficient_blume_capelSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_alpha(threshold_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type threshold_beta(threshold_betaSEXP);
    Rcpp::traits::input_parameter< const bool >::type na_impute(na_imputeSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type missing_index(missing_indexSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type is_ordinal_variable(is_ordinal_variableSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type reference_category(reference_categorySEXP);
    Rcpp::traits::input_parameter< const bool >::type save(saveSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type edge_selection(edge_selectionSEXP);
    Rcpp::traits::input_parameter< bool >::type mala(malaSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_sampler(observations, indicator, interactions, thresholds, num_categories, interaction_scale, proposal_sd, proposal_sd_blumecapel, edge_prior, theta, beta_bernoulli_alpha, beta_bernoulli_beta, dirichlet_alpha, lambda, Index, iter, burnin, num_obs_categories, sufficient_blume_capel, threshold_alpha, threshold_beta, na_impute, missing_index, is_ordinal_variable, reference_category, save, display_progress, edge_selection, mala));
    return rcpp_result_gen;
END_RCPP
}
// compare_anova_gibbs_sampler
List compare_anova_gibbs_sampler(arma::imat& observations, const arma::imat& main_effect_indices, const arma::imat& pairwise_effect_indices, const arma::mat& projection, const arma::imat& num_categories, const int num_groups, const arma::imat& group_indices, const double interaction_scale, const double pairwise_difference_scale, const double main_difference_scale, const String& pairwise_difference_prior, const String& main_difference_prior, arma::mat& inclusion_probability_difference, const double pairwise_beta_bernoulli_alpha, const double pairwise_beta_bernoulli_beta, const double main_beta_bernoulli_alpha, const double main_beta_bernoulli_beta, const arma::imat& Index, const int iter, const int burnin, List& num_obs_categories, List& sufficient_blume_capel, const double prior_threshold_alpha, const double prior_threshold_beta, const bool na_impute, const arma::imat& missing_data_indices, const arma::uvec& is_ordinal_variable, const arma::ivec& baseline_category, const bool independent_thresholds, const bool save_main, const bool save_pairwise, const bool save_indicator, const bool display_progress, bool difference_selection);
RcppExport SEXP _bgms_compare_anova_gibbs_sampler(SEXP observationsSEXP, SEXP main_effect_indicesSEXP, SEXP pairwise_effect_indicesSEXP, SEXP projectionSEXP, SEXP num_categoriesSEXP, SEXP num_groupsSEXP, SEXP group_indicesSEXP, SEXP interaction_scaleSEXP, SEXP pairwise_difference_scaleSEXP, SEXP main_difference_scaleSEXP, SEXP pairwise_difference_priorSEXP, SEXP main_difference_priorSEXP, SEXP inclusion_probability_differenceSEXP, SEXP pairwise_beta_bernoulli_alphaSEXP, SEXP pairwise_beta_bernoulli_betaSEXP, SEXP main_beta_bernoulli_alphaSEXP, SEXP main_beta_bernoulli_betaSEXP, SEXP IndexSEXP, SEXP iterSEXP, SEXP burninSEXP, SEXP num_obs_categoriesSEXP, SEXP sufficient_blume_capelSEXP, SEXP prior_threshold_alphaSEXP, SEXP prior_threshold_betaSEXP, SEXP na_imputeSEXP, SEXP missing_data_indicesSEXP, SEXP is_ordinal_variableSEXP, SEXP baseline_categorySEXP, SEXP independent_thresholdsSEXP, SEXP save_mainSEXP, SEXP save_pairwiseSEXP, SEXP save_indicatorSEXP, SEXP display_progressSEXP, SEXP difference_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::imat& >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type main_effect_indices(main_effect_indicesSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type pairwise_effect_indices(pairwise_effect_indicesSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type projection(projectionSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type num_categories(num_categoriesSEXP);
    Rcpp::traits::input_parameter< const int >::type num_groups(num_groupsSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type group_indices(group_indicesSEXP);
    Rcpp::traits::input_parameter< const double >::type interaction_scale(interaction_scaleSEXP);
    Rcpp::traits::input_parameter< const double >::type pairwise_difference_scale(pairwise_difference_scaleSEXP);
    Rcpp::traits::input_parameter< const double >::type main_difference_scale(main_difference_scaleSEXP);
    Rcpp::traits::input_parameter< const String& >::type pairwise_difference_prior(pairwise_difference_priorSEXP);
    Rcpp::traits::input_parameter< const String& >::type main_difference_prior(main_difference_priorSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type inclusion_probability_difference(inclusion_probability_differenceSEXP);
    Rcpp::traits::input_parameter< const double >::type pairwise_beta_bernoulli_alpha(pairwise_beta_bernoulli_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type pairwise_beta_bernoulli_beta(pairwise_beta_bernoulli_betaSEXP);
    Rcpp::traits::input_parameter< const double >::type main_beta_bernoulli_alpha(main_beta_bernoulli_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type main_beta_bernoulli_beta(main_beta_bernoulli_betaSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< const int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< const int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< List& >::type num_obs_categories(num_obs_categoriesSEXP);
    Rcpp::traits::input_parameter< List& >::type sufficient_blume_capel(sufficient_blume_capelSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_threshold_alpha(prior_threshold_alphaSEXP);
    Rcpp::traits::input_parameter< const double >::type prior_threshold_beta(prior_threshold_betaSEXP);
    Rcpp::traits::input_parameter< const bool >::type na_impute(na_imputeSEXP);
    Rcpp::traits::input_parameter< const arma::imat& >::type missing_data_indices(missing_data_indicesSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type is_ordinal_variable(is_ordinal_variableSEXP);
    Rcpp::traits::input_parameter< const arma::ivec& >::type baseline_category(baseline_categorySEXP);
    Rcpp::traits::input_parameter< const bool >::type independent_thresholds(independent_thresholdsSEXP);
    Rcpp::traits::input_parameter< const bool >::type save_main(save_mainSEXP);
    Rcpp::traits::input_parameter< const bool >::type save_pairwise(save_pairwiseSEXP);
    Rcpp::traits::input_parameter< const bool >::type save_indicator(save_indicatorSEXP);
    Rcpp::traits::input_parameter< const bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type difference_selection(difference_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_anova_gibbs_sampler(observations, main_effect_indices, pairwise_effect_indices, projection, num_categories, num_groups, group_indices, interaction_scale, pairwise_difference_scale, main_difference_scale, pairwise_difference_prior, main_difference_prior, inclusion_probability_difference, pairwise_beta_bernoulli_alpha, pairwise_beta_bernoulli_beta, main_beta_bernoulli_alpha, main_beta_bernoulli_beta, Index, iter, burnin, num_obs_categories, sufficient_blume_capel, prior_threshold_alpha, prior_threshold_beta, na_impute, missing_data_indices, is_ordinal_variable, baseline_category, independent_thresholds, save_main, save_pairwise, save_indicator, display_progress, difference_selection));
    return rcpp_result_gen;
END_RCPP
}
// compute_Vn_mfm_sbm
NumericVector compute_Vn_mfm_sbm(int no_variables, double dirichlet_alpha, int t_max, double lambda);
RcppExport SEXP _bgms_compute_Vn_mfm_sbm(SEXP no_variablesSEXP, SEXP dirichlet_alphaSEXP, SEXP t_maxSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type no_variables(no_variablesSEXP);
    Rcpp::traits::input_parameter< double >::type dirichlet_alpha(dirichlet_alphaSEXP);
    Rcpp::traits::input_parameter< int >::type t_max(t_maxSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(compute_Vn_mfm_sbm(no_variables, dirichlet_alpha, t_max, lambda));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bgms_sample_omrf_gibbs", (DL_FUNC) &_bgms_sample_omrf_gibbs, 6},
    {"_bgms_sample_bcomrf_gibbs", (DL_FUNC) &_bgms_sample_bcomrf_gibbs, 8},
    {"_bgms_gibbs_sampler", (DL_FUNC) &_bgms_gibbs_sampler, 29},
    {"_bgms_compare_anova_gibbs_sampler", (DL_FUNC) &_bgms_compare_anova_gibbs_sampler, 34},
    {"_bgms_compute_Vn_mfm_sbm", (DL_FUNC) &_bgms_compute_Vn_mfm_sbm, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bgms(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
