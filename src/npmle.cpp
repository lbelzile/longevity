// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


//' Turnbull EM algorithm (low storage implementation)
//'
//' @param x vector of ordered unique elements at which to estimate probabilities
//' @param truncLow integer vector in [0, J-1]; the lowest \code{x} such that \eqn{a_i < x} (left-truncation)
//' @param truncUpp integer vector in [0, J-1]; the largest \code{x} such that \eqn{b_i > x} (right-truncation)
//' @param censLow integer vector in [0, J-1]; the lowest \code{x} value at which the observation could be observed to fail (left censoring)
//' @param censUpp integer vector in [0, J-1]; the largest \code{x} value at which the observation could fail (right-censoring)
//' @param cens logical; if \code{FALSE}, then \code{censUpp = censLow} and a particular update can be avoided in the EM algorithm
//' @param tol tolerance level for terminating the EM algorithm
//' @param zerotol tolerance level for setting a probability to zero
//' @param maxiter maximum number of iteration for the EM algorithm
//' @return a list with the probabilities and the standard errors
//' @keywords internal
// [[Rcpp::export(.turnbull_em)]]
Rcpp::List turnbullem(
    arma::vec x,
    arma::uvec censLow,
    arma::uvec censUpp,
    arma::uvec truncLow,
    arma::uvec truncUpp,
    bool cens = true,
    bool trunc = true,
    double tol = 1e-12,
    int maxiter = 1e5){
  arma::uword J = x.n_elem;
  arma::uword n = censLow.n_elem;
  if(trunc){
    if(truncLow.n_elem != n | truncUpp.n_elem != n){
      Rcpp::stop("All vectors of positions should be of the same length.");
    }
    if(any(truncLow > truncUpp)){
      Rcpp::stop("Invalid arguments: elements are in reverse order.");
    }
    if(any(truncLow < 0) || any(truncUpp > J - 1)){
      Rcpp::stop("Invalid arguments: the vectors of indices are not in [0, J-1].");
    }
  }
  if(cens){
    if(censUpp.n_elem != n){
      Rcpp::stop("All vectors of positions should be of the same length.");
    }
    if(any(censLow > censUpp)){
      Rcpp::stop("Invalid arguments: some censoring intervals have lower limits that exceed the upper limits.");
    }
    if(any(censLow < 0) || any(censUpp > J - 1)){
      Rcpp::stop("Invalid arguments: vectors of indices are not in [0, J-1].");
    }
  }
  arma::vec pCur(J);
  arma::vec pNew(J, arma::fill::zeros);
  pNew.fill(1/(double)J);
  // arma::vec ui(n, arma::fill::zeros);
  arma::vec uiCum(J, arma::fill::zeros);
  double sum_p = 0;
  int niter = 0;
  bool convergence = false;
  while(!convergence && niter < maxiter){
    if(arma::max(arma::abs(pCur - pNew)) < tol){
      convergence = true;
    }
    pCur = pNew;
    Rcpp::checkUserInterrupt();
    uiCum.zeros();
    for(arma::uword i = 0; i < n; ++i){
      // Censoring step - only
      if(cens){
        sum_p = arma::sum(pCur(arma::span(censLow(i), censUpp(i))));
        for(arma::uword j = censLow(i); j <= censUpp(i); ++j){
          uiCum(j) += pCur(j)/sum_p;
        }
      } else{
        uiCum(censLow(i)) += 1;
      }
      if(trunc){
        // Truncation step
        sum_p = arma::sum(pCur(arma::span(truncLow(i), truncUpp(i))));
        for(arma::uword j = 0; j < J; ++j){
          // TODO: check these...
          if(j < truncLow(i) || j > truncUpp(i)){
            uiCum(j) += pCur(j) / sum_p;
          }
        }
      }
    }
    pNew = uiCum / arma::sum(uiCum);
    niter ++;
    // Rcpp::Rcout << "Iteration "<< niter << ": "<< pNew << std::endl;
  }
  return Rcpp::List::create(Rcpp::Named("p") = pNew,
                            Rcpp::Named("conv") = convergence,
                            Rcpp::Named("niter") = niter);
}
