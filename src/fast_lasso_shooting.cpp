//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
SEXP armaCrossprod(arma::mat A, arma::mat B){
  arma::mat C = A.t() * B;
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
List fastLassoShooting(const arma::mat& XX2, const arma::colvec& Xy2, arma::colvec& beta_start,
                       const arma::colvec& lambda, const int max_iter, const double optTol) {
  int p = XX2.n_cols;
  int j;
  int m=0;
  double S0;
  arma::colvec beta_old=beta_start;
  arma::colvec beta_res=beta_start;
  
  while (m < max_iter)
  {
    beta_old = beta_res;
    for (j = 0; j < p; j++) {
      // Compute the Shoot and Update the variable
      S0 = arma::dot(XX2.row(j), beta_res) - XX2(j, j) * beta_res(j) - Xy2(j);
      
      if (S0 > lambda(j))
      {
        beta_res(j) = (lambda(j) - S0) / XX2(j, j);
      }
      else if (S0 < -1 * lambda(j))
      {
        beta_res(j) = (-1 * lambda(j) - S0) / XX2(j, j);
      }
      else if (std::abs(S0) <= lambda(j))
      {
        beta_res(j) = 0;
      }
    }
    if (sum(arma::abs(beta_res - beta_old)) < optTol)
      {
      break;
      }
    m = m + 1;
  }
  
  return List::create(Named("S0") = S0, Named("beta") = beta_res);
}
