#include <cmath>  // std::pow

#include <RcppArmadillo.h>
#include "functor.h"
#include "applic.h"
#include "roptim.h"
#include "samin.h"


// Function to create a correlation matrix from correlation coefficients
arma::mat makeSigma(arma::vec sigma, int d) {

     // Checking the inputs
     if(sigma.size()!= (d*(d-1))/2){
          Rcpp::stop("No match among sigma vector and d-dimension");
     }

     arma::mat Sigma(d, d, arma::fill::eye);

     int index = 0;
     for (int i = 0; i < d; i++) {
          for (int j = i + 1; j < d; j++) {
               Sigma(i, j) = sigma(index);
               Sigma(j, i) = sigma(index);
               index++;
          }
     }

     return Sigma;
}


// === Slow version but matches with THE R code // maybe can be slightly improved===
//[[Rcpp::export]]
double log_mvn_post_cor_sample(arma::mat y_hat_, // The number of observations are the column
                               arma::mat y_mat_,
                               arma::vec& sigmas,int d,
                               arma::vec& sigmas_0,
                               arma::mat& Sigma_0){

     // Defining some quantities;
     double likelihood_term = 0.0;
     // int n_sigma_ = (d*(d-1))/2;
     arma::mat Sigma_= makeSigma(sigmas,d);
     arma::mat res_ = (y_mat_-y_hat_).t();
     // cout << "Nrows: " << res_.n_rows << endl;
     for(int ii = 0; ii < res_.n_cols; ii++){
          arma::vec aux_res_ = res_.col(ii);
          likelihood_term = likelihood_term +  arma::as_scalar(aux_res_.t()*arma::inv(Sigma_)*aux_res_);
     }

     // arma::mat Sigma_0_((d*d-d)/2, (d*d-d)/2, arma::fill::eye);
     if((Sigma_0.n_cols!=0.5*(d*d-d)) | (Sigma_0.n_rows!=0.5*(d*d-d))){
             Rcpp::stop("Insert a valid variance prior.");
     }

     arma::mat Sigma_0_ = Sigma_0;

     arma::vec sigmas_diff = (sigmas-sigmas_0);

     // arma::cout << " All good" << arma::endl;
     return -0.5*y_mat_.n_rows*log(det(Sigma_))-0.5*likelihood_term - 0.5*arma::as_scalar(sigmas_diff.t()*arma::inv(Sigma_0_)*sigmas_diff);


     // return likelihood_term;


};

class LogLikePost : public roptim::Functor {
public:

        // Storing the data to initialise the model
        double likelihood_term = 0.0;
        // int n_sigma_ = (d*(d-1))/2;
        arma::mat Sigma;
        arma::mat res;
        int d;
        arma::mat Sigma_0;
        arma::vec sigma_0;
        arma::mat y_mat;
        arma::mat y_hat;
        arma::vec sigmas_diff;

        // Creating the data for the LogLikePost
        LogLikePost(int d_,
              arma::vec sigma_0_,
              arma::mat Sigma_0_,
              arma::mat y_mat_,
              arma::mat y_hat_){

                // Defining some of those elements
                d = d_;
                y_mat = y_mat_;
                y_hat = y_hat_;
                res = (y_mat_-y_hat_).t();
                Sigma_0 = arma::mat((d*d-d)/2, (d*d-d)/2, arma::fill::eye);
                sigma_0 = sigma_0_;



        };
        double operator()(const arma::vec &x) override {

                // Defining the logquantity;
                Sigma = makeSigma(x,d);
                likelihood_term = 0.0;
                // cout << "Nrows: " << res_.n_rows << endl;
                for(int ii = 0; ii < res.n_cols; ii++){
                        arma::vec aux_res_ = res.col(ii);
                        likelihood_term = likelihood_term +  arma::as_scalar(aux_res_.t()*arma::inv(Sigma)*aux_res_);
                }

                // arma::mat Sigma_0 = makeSigma(sigma_0,d);
                arma::vec sigmas_diff = (x-sigma_0);

                double final_value = -0.5*y_mat.n_rows*log(det(Sigma))-0.5*likelihood_term - 0.5*arma::as_scalar(sigmas_diff.t()*arma::inv(Sigma_0)*sigmas_diff);
                // arma::cout << "Final value is: " << final_value << arma::endl;
                // return -0.5*y_mat.n_rows*log(det(Sigma))-0.5*likelihood_term - 0.5*arma::as_scalar(sigmas_diff.t()*arma::inv(Sigma_0)*sigmas_diff);
                return -final_value;

        }


};


// [[Rcpp::export]]
arma::vec sigma_draw_cpp(int d,
                arma::vec sigma_0,
                arma::mat Sigma_0,
                arma::vec sigma_init_optim_,
                arma::mat y_mat,
                arma::mat y_hat,
                double df) {

        // Setting the constructor of the LogLikePost
        LogLikePost optim_obj(d,
                       sigma_0,
                       Sigma_0,
                       y_mat,
                       y_hat);

        roptim::Roptim<LogLikePost> opt("BFGS");
        // opt.control.trace = 1;
        opt.control.maxit = 3;
        opt.set_hessian(true);
        // Setting the obj.
        // arma::vec sigma_init;
        opt.minimize(optim_obj, sigma_init_optim_);
        arma::vec mu_ = opt.par();

        arma::mat S = 1.0*arma::inv(opt.hessian());
        arma::vec r_sample((d*d-d)/2);
        for(int ii = 0; ii < r_sample.size(); ii++){
                r_sample(ii) = R::rnorm(0,1.0);
        }

        // Generating a vec sample
        // arma::cout << S << arma::endl;
        arma::vec y = arma::chol(S).t()*r_sample;
        double u = R::rchisq(df);
        return sqrt(df/u)*y + mu_;

}

//[[Rcpp::export]]
arma::mat sigma_sampler(int nmcmc,
                        int d,
                        arma::vec sigma_0,
                        arma::mat Sigma_0,
                        arma::vec sigma_init_optim,
                        arma::mat y_mat,
                        arma::mat y_hat,
                        double df){

        arma::mat sigma_post((d*d-d)/2,nmcmc,arma::fill::zeros);
        sigma_post.col(0) = sigma_0;

        if(sigma_init_optim.size()!=(d*d-d)/2){
                Rcpp::stop("Insert a valid sigma initialisation");
        }
        // Testing the loglikelihood
        for(int ii=1; ii < nmcmc; ii++){

                // Generating the sampler
                arma::vec draw = sigma_draw_cpp(d,
                                                sigma_0,
                                                Sigma_0,
                                                sigma_init_optim,
                                                y_mat,
                                                y_hat,df);

                double prob_;
                arma::mat new_Sigma0 = makeSigma(draw,d);
                arma::vec eigen_val ;
                arma::mat eigen_vec;
                arma::eig_sym(eigen_val,eigen_vec,new_Sigma0);

                if(arma::min(eigen_val)>0  & max(abs(draw))<1){
                        arma::vec prev_sigma = sigma_post.col(ii-1);
                        prob_ = exp(log_mvn_post_cor_sample(y_hat,y_mat,draw,d,sigma_0,Sigma_0) -
                                log_mvn_post_cor_sample(y_hat,y_mat,prev_sigma,d,sigma_0,Sigma_0));
                } else {
                        prob_ = 0;
                }

                if(R::runif(0,1)<prob_){
                        sigma_post.col(ii) = draw;
                } else {
                        sigma_post.col(ii) = sigma_post.col(ii-1);
                }
        }

        // Returning the posterior
        return sigma_post;
}


