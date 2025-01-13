// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

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
List gibbs_sampler(IntegerMatrix observations, IntegerMatrix indicator, NumericMatrix interactions, NumericMatrix thresholds, IntegerVector no_categories, double interaction_scale, NumericMatrix proposal_sd, NumericMatrix proposal_sd_blumecapel, String edge_prior, NumericMatrix theta, double beta_bernoulli_alpha, double beta_bernoulli_beta, double dirichlet_alpha, double lambda, IntegerMatrix Index, int iter, int burnin, IntegerMatrix n_cat_obs, IntegerMatrix sufficient_blume_capel, double threshold_alpha, double threshold_beta, bool na_impute, IntegerMatrix missing_index, LogicalVector variable_bool, IntegerVector reference_category, bool save, bool display_progress, bool edge_selection);
RcppExport SEXP _bgms_gibbs_sampler(SEXP observationsSEXP, SEXP indicatorSEXP, SEXP interactionsSEXP, SEXP thresholdsSEXP, SEXP no_categoriesSEXP, SEXP interaction_scaleSEXP, SEXP proposal_sdSEXP, SEXP proposal_sd_blumecapelSEXP, SEXP edge_priorSEXP, SEXP thetaSEXP, SEXP beta_bernoulli_alphaSEXP, SEXP beta_bernoulli_betaSEXP, SEXP dirichlet_alphaSEXP, SEXP lambdaSEXP, SEXP IndexSEXP, SEXP iterSEXP, SEXP burninSEXP, SEXP n_cat_obsSEXP, SEXP sufficient_blume_capelSEXP, SEXP threshold_alphaSEXP, SEXP threshold_betaSEXP, SEXP na_imputeSEXP, SEXP missing_indexSEXP, SEXP variable_boolSEXP, SEXP reference_categorySEXP, SEXP saveSEXP, SEXP display_progressSEXP, SEXP edge_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations(observationsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type indicator(indicatorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type interactions(interactionsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories(no_categoriesSEXP);
    Rcpp::traits::input_parameter< double >::type interaction_scale(interaction_scaleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proposal_sd(proposal_sdSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type proposal_sd_blumecapel(proposal_sd_blumecapelSEXP);
    Rcpp::traits::input_parameter< String >::type edge_prior(edge_priorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< double >::type beta_bernoulli_alpha(beta_bernoulli_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta_bernoulli_beta(beta_bernoulli_betaSEXP);
    Rcpp::traits::input_parameter< double >::type dirichlet_alpha(dirichlet_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type n_cat_obs(n_cat_obsSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type sufficient_blume_capel(sufficient_blume_capelSEXP);
    Rcpp::traits::input_parameter< double >::type threshold_alpha(threshold_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type threshold_beta(threshold_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type na_impute(na_imputeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type missing_index(missing_indexSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type variable_bool(variable_boolSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type reference_category(reference_categorySEXP);
    Rcpp::traits::input_parameter< bool >::type save(saveSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type edge_selection(edge_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(gibbs_sampler(observations, indicator, interactions, thresholds, no_categories, interaction_scale, proposal_sd, proposal_sd_blumecapel, edge_prior, theta, beta_bernoulli_alpha, beta_bernoulli_beta, dirichlet_alpha, lambda, Index, iter, burnin, n_cat_obs, sufficient_blume_capel, threshold_alpha, threshold_beta, na_impute, missing_index, variable_bool, reference_category, save, display_progress, edge_selection));
    return rcpp_result_gen;
END_RCPP
}
// compare_gibbs_sampler
List compare_gibbs_sampler(IntegerMatrix observations_gr1, IntegerMatrix observations_gr2, IntegerVector no_categories_gr1, IntegerVector no_categories_gr2, double interaction_scale, double pairwise_difference_scale, double main_difference_scale, String pairwise_difference_prior, String main_difference_prior, NumericMatrix inclusion_probability_difference, double pairwise_beta_bernoulli_alpha, double pairwise_beta_bernoulli_beta, double main_beta_bernoulli_alpha, double main_beta_bernoulli_beta, IntegerMatrix Index, int iter, int burnin, IntegerMatrix n_cat_obs_gr1, IntegerMatrix n_cat_obs_gr2, IntegerMatrix sufficient_blume_capel_gr1, IntegerMatrix sufficient_blume_capel_gr2, double threshold_alpha, double threshold_beta, bool na_impute, IntegerMatrix missing_index_gr1, IntegerMatrix missing_index_gr2, LogicalVector ordinal_variable, IntegerVector reference_category, bool independent_thresholds, bool save, bool display_progress, bool difference_selection);
RcppExport SEXP _bgms_compare_gibbs_sampler(SEXP observations_gr1SEXP, SEXP observations_gr2SEXP, SEXP no_categories_gr1SEXP, SEXP no_categories_gr2SEXP, SEXP interaction_scaleSEXP, SEXP pairwise_difference_scaleSEXP, SEXP main_difference_scaleSEXP, SEXP pairwise_difference_priorSEXP, SEXP main_difference_priorSEXP, SEXP inclusion_probability_differenceSEXP, SEXP pairwise_beta_bernoulli_alphaSEXP, SEXP pairwise_beta_bernoulli_betaSEXP, SEXP main_beta_bernoulli_alphaSEXP, SEXP main_beta_bernoulli_betaSEXP, SEXP IndexSEXP, SEXP iterSEXP, SEXP burninSEXP, SEXP n_cat_obs_gr1SEXP, SEXP n_cat_obs_gr2SEXP, SEXP sufficient_blume_capel_gr1SEXP, SEXP sufficient_blume_capel_gr2SEXP, SEXP threshold_alphaSEXP, SEXP threshold_betaSEXP, SEXP na_imputeSEXP, SEXP missing_index_gr1SEXP, SEXP missing_index_gr2SEXP, SEXP ordinal_variableSEXP, SEXP reference_categorySEXP, SEXP independent_thresholdsSEXP, SEXP saveSEXP, SEXP display_progressSEXP, SEXP difference_selectionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations_gr1(observations_gr1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type observations_gr2(observations_gr2SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories_gr1(no_categories_gr1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type no_categories_gr2(no_categories_gr2SEXP);
    Rcpp::traits::input_parameter< double >::type interaction_scale(interaction_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type pairwise_difference_scale(pairwise_difference_scaleSEXP);
    Rcpp::traits::input_parameter< double >::type main_difference_scale(main_difference_scaleSEXP);
    Rcpp::traits::input_parameter< String >::type pairwise_difference_prior(pairwise_difference_priorSEXP);
    Rcpp::traits::input_parameter< String >::type main_difference_prior(main_difference_priorSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type inclusion_probability_difference(inclusion_probability_differenceSEXP);
    Rcpp::traits::input_parameter< double >::type pairwise_beta_bernoulli_alpha(pairwise_beta_bernoulli_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type pairwise_beta_bernoulli_beta(pairwise_beta_bernoulli_betaSEXP);
    Rcpp::traits::input_parameter< double >::type main_beta_bernoulli_alpha(main_beta_bernoulli_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type main_beta_bernoulli_beta(main_beta_bernoulli_betaSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type Index(IndexSEXP);
    Rcpp::traits::input_parameter< int >::type iter(iterSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type n_cat_obs_gr1(n_cat_obs_gr1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type n_cat_obs_gr2(n_cat_obs_gr2SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type sufficient_blume_capel_gr1(sufficient_blume_capel_gr1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type sufficient_blume_capel_gr2(sufficient_blume_capel_gr2SEXP);
    Rcpp::traits::input_parameter< double >::type threshold_alpha(threshold_alphaSEXP);
    Rcpp::traits::input_parameter< double >::type threshold_beta(threshold_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type na_impute(na_imputeSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type missing_index_gr1(missing_index_gr1SEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type missing_index_gr2(missing_index_gr2SEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type ordinal_variable(ordinal_variableSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type reference_category(reference_categorySEXP);
    Rcpp::traits::input_parameter< bool >::type independent_thresholds(independent_thresholdsSEXP);
    Rcpp::traits::input_parameter< bool >::type save(saveSEXP);
    Rcpp::traits::input_parameter< bool >::type display_progress(display_progressSEXP);
    Rcpp::traits::input_parameter< bool >::type difference_selection(difference_selectionSEXP);
    rcpp_result_gen = Rcpp::wrap(compare_gibbs_sampler(observations_gr1, observations_gr2, no_categories_gr1, no_categories_gr2, interaction_scale, pairwise_difference_scale, main_difference_scale, pairwise_difference_prior, main_difference_prior, inclusion_probability_difference, pairwise_beta_bernoulli_alpha, pairwise_beta_bernoulli_beta, main_beta_bernoulli_alpha, main_beta_bernoulli_beta, Index, iter, burnin, n_cat_obs_gr1, n_cat_obs_gr2, sufficient_blume_capel_gr1, sufficient_blume_capel_gr2, threshold_alpha, threshold_beta, na_impute, missing_index_gr1, missing_index_gr2, ordinal_variable, reference_category, independent_thresholds, save, display_progress, difference_selection));
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
    {"_bgms_gibbs_sampler", (DL_FUNC) &_bgms_gibbs_sampler, 28},
    {"_bgms_compare_gibbs_sampler", (DL_FUNC) &_bgms_compare_gibbs_sampler, 32},
    {"_bgms_compute_Vn_mfm_sbm", (DL_FUNC) &_bgms_compute_Vn_mfm_sbm, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_bgms(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
