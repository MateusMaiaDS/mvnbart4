// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// log_mvn_post_cor_sample
double log_mvn_post_cor_sample(arma::mat y_hat_, arma::mat y_mat_, arma::vec& sigmas, int d, arma::vec& sigmas_0);
RcppExport SEXP _mvnbart4_log_mvn_post_cor_sample(SEXP y_hat_SEXP, SEXP y_mat_SEXP, SEXP sigmasSEXP, SEXP dSEXP, SEXP sigmas_0SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type y_hat_(y_hat_SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_mat_(y_mat_SEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigmas(sigmasSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec& >::type sigmas_0(sigmas_0SEXP);
    rcpp_result_gen = Rcpp::wrap(log_mvn_post_cor_sample(y_hat_, y_mat_, sigmas, d, sigmas_0));
    return rcpp_result_gen;
END_RCPP
}
// sigma_draw_cpp
arma::vec sigma_draw_cpp(int d, arma::vec sigma_0, arma::mat y_mat, arma::mat y_hat, double df);
RcppExport SEXP _mvnbart4_sigma_draw_cpp(SEXP dSEXP, SEXP sigma_0SEXP, SEXP y_matSEXP, SEXP y_hatSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_0(sigma_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_mat(y_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_hat(y_hatSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_draw_cpp(d, sigma_0, y_mat, y_hat, df));
    return rcpp_result_gen;
END_RCPP
}
// sigma_sampler
arma::mat sigma_sampler(int nmcmc, int d, arma::vec sigma_0, arma::mat y_mat, arma::mat y_hat, double df);
RcppExport SEXP _mvnbart4_sigma_sampler(SEXP nmcmcSEXP, SEXP dSEXP, SEXP sigma_0SEXP, SEXP y_matSEXP, SEXP y_hatSEXP, SEXP dfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nmcmc(nmcmcSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_0(sigma_0SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_mat(y_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_hat(y_hatSEXP);
    Rcpp::traits::input_parameter< double >::type df(dfSEXP);
    rcpp_result_gen = Rcpp::wrap(sigma_sampler(nmcmc, d, sigma_0, y_mat, y_hat, df));
    return rcpp_result_gen;
END_RCPP
}
// log_dmvn
double log_dmvn(arma::vec& x, arma::mat& Sigma);
RcppExport SEXP _mvnbart4_log_dmvn(SEXP xSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(log_dmvn(x, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// cppbart
Rcpp::List cppbart(arma::mat x_train, arma::mat y_mat, arma::mat x_test, arma::mat x_cut, int n_tree, int node_min_size, int n_mcmc, int n_burn, arma::mat Sigma_init, arma::vec mu_init, arma::vec sigma_mu, double alpha, double beta, double nu, arma::mat S_0_wish, arma::vec A_j_vec, bool update_Sigma, bool conditional_bool);
RcppExport SEXP _mvnbart4_cppbart(SEXP x_trainSEXP, SEXP y_matSEXP, SEXP x_testSEXP, SEXP x_cutSEXP, SEXP n_treeSEXP, SEXP node_min_sizeSEXP, SEXP n_mcmcSEXP, SEXP n_burnSEXP, SEXP Sigma_initSEXP, SEXP mu_initSEXP, SEXP sigma_muSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP nuSEXP, SEXP S_0_wishSEXP, SEXP A_j_vecSEXP, SEXP update_SigmaSEXP, SEXP conditional_boolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x_train(x_trainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_mat(y_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_test(x_testSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_cut(x_cutSEXP);
    Rcpp::traits::input_parameter< int >::type n_tree(n_treeSEXP);
    Rcpp::traits::input_parameter< int >::type node_min_size(node_min_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type n_mcmc(n_mcmcSEXP);
    Rcpp::traits::input_parameter< int >::type n_burn(n_burnSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_init(Sigma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_mu(sigma_muSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type S_0_wish(S_0_wishSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type A_j_vec(A_j_vecSEXP);
    Rcpp::traits::input_parameter< bool >::type update_Sigma(update_SigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type conditional_bool(conditional_boolSEXP);
    rcpp_result_gen = Rcpp::wrap(cppbart(x_train, y_mat, x_test, x_cut, n_tree, node_min_size, n_mcmc, n_burn, Sigma_init, mu_init, sigma_mu, alpha, beta, nu, S_0_wish, A_j_vec, update_Sigma, conditional_bool));
    return rcpp_result_gen;
END_RCPP
}
// cppbart_CLASS
Rcpp::List cppbart_CLASS(arma::mat x_train, arma::mat y_mat, arma::mat x_test, arma::mat x_cut, int n_tree, int node_min_size, int n_mcmc, int n_burn, arma::mat Sigma_init, arma::vec mu_init, arma::vec sigma_mu, double alpha, double beta);
RcppExport SEXP _mvnbart4_cppbart_CLASS(SEXP x_trainSEXP, SEXP y_matSEXP, SEXP x_testSEXP, SEXP x_cutSEXP, SEXP n_treeSEXP, SEXP node_min_sizeSEXP, SEXP n_mcmcSEXP, SEXP n_burnSEXP, SEXP Sigma_initSEXP, SEXP mu_initSEXP, SEXP sigma_muSEXP, SEXP alphaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x_train(x_trainSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type y_mat(y_matSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_test(x_testSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x_cut(x_cutSEXP);
    Rcpp::traits::input_parameter< int >::type n_tree(n_treeSEXP);
    Rcpp::traits::input_parameter< int >::type node_min_size(node_min_sizeSEXP);
    Rcpp::traits::input_parameter< int >::type n_mcmc(n_mcmcSEXP);
    Rcpp::traits::input_parameter< int >::type n_burn(n_burnSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma_init(Sigma_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type mu_init(mu_initSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type sigma_mu(sigma_muSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(cppbart_CLASS(x_train, y_mat, x_test, x_cut, n_tree, node_min_size, n_mcmc, n_burn, Sigma_init, mu_init, sigma_mu, alpha, beta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_mvnbart4_log_mvn_post_cor_sample", (DL_FUNC) &_mvnbart4_log_mvn_post_cor_sample, 5},
    {"_mvnbart4_sigma_draw_cpp", (DL_FUNC) &_mvnbart4_sigma_draw_cpp, 5},
    {"_mvnbart4_sigma_sampler", (DL_FUNC) &_mvnbart4_sigma_sampler, 6},
    {"_mvnbart4_log_dmvn", (DL_FUNC) &_mvnbart4_log_dmvn, 2},
    {"_mvnbart4_cppbart", (DL_FUNC) &_mvnbart4_cppbart, 18},
    {"_mvnbart4_cppbart_CLASS", (DL_FUNC) &_mvnbart4_cppbart_CLASS, 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_mvnbart4(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
