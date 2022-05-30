// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' Turnbull's sets
//'
//' Given truncation and censoring sets,
//' construct disjoint increasing intervals whose
//' left and right endpoints lie in L and R
//' and which contain no other members of L and R
//'
//' @export
//' @keywords internal
// [[Rcpp::export(.turnbull_intervals)]]
arma::mat turnbull_intervals(
    arma::vec Lset,
    arma::vec Rset
){
  // Sort and get unique elements of left and right sets
  Lset = arma::unique(Lset);
  Rset = arma::unique(Rset);
  // Rcpp::Rcout << min(Lset) << " and " << max(Rset) << std::endl;
  int n = std::min(Lset.n_elem, Rset.n_elem);
  // Create container, reduce size latter
  arma::mat turnset(n, 2);
  arma::uword Lind = 0;
  arma::uword Rind = 0;
  arma::uword i = 0;
  // Case of first instance
  // Start with Lind = 0 since min(Rset) > min(Lset)
  for(arma::uword j = Rind; j < Rset.n_elem; ++j){
    // find the next value of R > L
    // stop at first occurrence
    Rind = j;
    if(Rset(j) >= Lset(Lind)){
      // Increment left set if there was a gap
      for(arma::uword k = Lind; k < Lset.n_elem; ++k){
        if(Lset(k) <= Rset(j)){
          Lind = k;
        } else{
          // Stop as soon as we exceed current upper bound
          break;
        }
      }
      break;
    }
  }
  bool keepgoing = Lind < Lset.n_elem - 1 & Rind < Rset.n_elem - 1;
 // Loop as long as there are elements left to add to list
   while(keepgoing){
    // Rcpp::Rcout << Lind << " and " << Rind << std::endl;
    turnset(i,0) = Lset(Lind);
    turnset(i,1) = Rset(Rind);
    i++;
    // If the largest remaining value is left, then stop
    if(Lset(Lset.n_elem - 1) <= Rset(Rind)){
      keepgoing = false;
      break;
    }
    for(arma::uword j = Lind; j < Lset.n_elem; ++j){
      Lind = j;
      if(Lset(j) > Rset(Rind)){
        break;
      }
    }
    if(Rset(Rset.n_elem - 1) < Lset(Lind)){
      keepgoing = false;
      break;
    }
    for(arma::uword j = Rind; j < Rset.n_elem; ++j){
      // find the next value of R > L
      // stop at first occurrence
      Rind = j;
      // Check this is feasible
      if(Rset(j) >= Lset(Lind)){
        // Increment left set if there was a gap
        for(arma::uword k = Lind; k < Lset.n_elem; ++k){
         if(Lset(k) <= Rset(j)){
            Lind = k;
         } else{
           // Stop as soon as we exceed current upper bound
           break;
         }
        }
        break;
      }
    }
  }
  return turnset.head_rows(i);
}

//' Identification sets
//'
//' @param tsets Turnbull's sets
//' @param lcens numeric vector of left censoring
//' @param rcens numeric vector of right censoring
//' @param ltrunc numeric vector of left truncation
//' @param rtrunc numeric vector of right truncation
//' @param trunc logical are observation truncated?
//' @export
//' @keywords internal
// [[Rcpp::export(.censTruncLimits)]]
Rcpp::List censTruncLimits(
  arma::mat tsets,
  int n,
  arma::vec lcens,
  arma::vec rcens,
  arma::vec ltrunc,
  arma::vec rtrunc,
  bool trunc
){
  // Initialize containers
  arma::uvec censLow(n);
  arma::uvec censUpp(n);
  arma::uvec truncLow(n);
  arma::uvec truncUpp(n);
  if(lcens.n_elem != n | rcens.n_elem != n){
    Rcpp::stop("All vectors of censoring intervals should be of the same length.");
  }
  int J = tsets.n_rows;
  if(tsets.n_cols != 2){
    Rcpp::stop("\"tsets\" should be a matrix with two columns.");
  }
  // Create values that are impossible
  censLow.fill(J);
  censUpp.fill(J);
  truncLow.fill(J);
  truncUpp.fill(J);

for(arma::uword i = 0; i < n; ++i){
  for(arma::uword j = 0; j < J; ++j){
   if(lcens(i) <= tsets(j,0) & tsets(j,1) <= rcens(i)){
     // Initialized to J, so pick the smallest integer
     if(j < censLow(i)){
       censLow(i) = j;
     }
     // Since sets are increasing, pick the largest
     censUpp(i) = j;
   }
  }
}
if(!trunc){
  for(arma::uword i = 0; i < n; ++i){
    truncLow(i) = 0;
    truncUpp(i) = J - 1;
  }
} else{
  if(ltrunc.n_elem != n | rtrunc.n_elem != n){
    Rcpp::stop("All vectors of censoring intervals should be of the same length.");
  }
  for(arma::uword i = 0; i < n; ++i){
    for(arma::uword j = 0; j < J; ++j){
    if(ltrunc(i) <= tsets(j,0) & tsets(j,1) <= rtrunc(i)){
        // Initialized to J, so pick the smallest integer
        if(j < truncLow(i)){
         truncLow(i) = j;
        }
        // Since sets are increasing, pick the largest
        truncUpp(i) = j;
      }
    }
  }
}
// Intervals; if the value is J, then intervals are empty
return Rcpp::List::create(Rcpp::Named("truncLow") = truncLow,
                          Rcpp::Named("truncUpp") = truncUpp,
                          Rcpp::Named("censLow") = censLow,
                          Rcpp::Named("censUpp") = censUpp);

}

//' Turnbull EM algorithm (low storage implementation)
//'
//' @param tsets Turnbull's sets
//' @param lcens numeric vector of left censoring
//' @param rcens numeric vector of right censoring
//' @param ltrunc numeric vector of left truncation
//' @param rtrunc numeric vector of right truncation
//' @param cens logical; if \code{FALSE}, then \code{censUpp = censLow} and a particular update can be avoided in the EM algorithm
//' @param tol tolerance level for terminating the EM algorithm
//' @param maxiter maximum number of iteration for the EM algorithm
//' @return a list with the probabilities and the standard errors
//' @keywords internal
// [[Rcpp::export(.turnbull_em)]]
Rcpp::List turnbullem(
    arma::mat tsets,
    int n,
    arma::vec lcens,
    arma::vec rcens,
    arma::vec ltrunc,
    arma::vec rtrunc,
    bool cens = true,
    bool trunc = true,
    double tol = 1e-12,
    double zerotol = 1e-40,
    int maxiter = 1e5){
  arma::uvec censLow(n);
  arma::uvec censUpp(n);
  arma::uvec truncLow(n);
  arma::uvec truncUpp(n);
  if(lcens.n_elem != n | rcens.n_elem != n){
    Rcpp::stop("All vectors of censoring intervals should be of the same length.");
  }
  if(any(lcens > rcens)){
      Rcpp::stop("Invalid arguments: \"lcens\" > \"rcens\" for some elements.");
  }
  int J = tsets.n_rows;
  if(tsets.n_cols != 2){
    Rcpp::stop("\"tsets\" should be a matrix with two columns.");
  }
  // Create values that are impossible
  censLow.fill(J);
  censUpp.fill(J);
  truncLow.fill(J);
  truncUpp.fill(J);

  for(arma::uword i = 0; i < n; ++i){
    for(arma::uword j = 0; j < J; ++j){
      if(lcens(i) <= tsets(j,0) & tsets(j,1) <= rcens(i)){
        // Initialized to J, so pick the smallest integer
        if(j < censLow(i)){
          censLow(i) = j;
        }
        // Since sets are increasing, pick the largest
        censUpp(i) = j;
      }
    }
  }
  if(trunc){
        if(ltrunc.n_elem != n | rtrunc.n_elem != n){
        Rcpp::stop("All vectors of positions should be of the same length.");
      }
      if(any(ltrunc > rtrunc)){
        Rcpp::stop("Invalid arguments: elements are in reverse order.");
      }
      for(arma::uword i = 0; i < n; ++i){
        for(arma::uword j = 0; j < J; ++j){
          if(ltrunc(i) <= tsets(j,0) & tsets(j,1) <= rtrunc(i)){
            // Initialized to J, so pick the smallest integer
            if(j < truncLow(i)){
              truncLow(i) = j;
            }
            // Since sets are increasing, pick the largest
            truncUpp(i) = j;
          }
        }
      }
  } else{
    for(arma::uword i = 0; i < n; ++i){
      truncLow(i) = 0;
      truncUpp(i) = J - 1;
    }
  }
  arma::vec pCur(J);
  arma::vec pNew(J, arma::fill::zeros);
  // Equiprobable initial value
  pNew.fill(1/(double)J);
  arma::vec uiCum(J, arma::fill::zeros);
  double sum_p = 0;
  double abstol = 0;
  int niter = 0;
  bool convergence = false;
  while(!convergence && niter < maxiter){
    abstol = arma::max(arma::abs(pCur - pNew));
    if(abstol < tol){
      convergence = true;
    }
    pCur = pNew;
    Rcpp::checkUserInterrupt();
    uiCum.zeros();
    for(arma::uword i = 0; i < n; ++i){
      // Censoring step - only
      if(censLow(i) < J){
        if(cens){
          sum_p = arma::sum(pCur(arma::span(censLow(i), censUpp(i))));
          for(arma::uword j = censLow(i); j <= censUpp(i); ++j){
            uiCum(j) += pCur(j)/sum_p;
          }
          } else{
          uiCum(censLow(i)) += 1;
        }
      }
      if(trunc & truncLow(i) < J){
        // Truncation step
        sum_p = arma::sum(pCur(arma::span(truncLow(i), truncUpp(i))));
        for(arma::uword j = 0; j < J; ++j){
          if(j < truncLow(i) || j > truncUpp(i)){
            uiCum(j) += pCur(j) / sum_p;
          }
        }
      }
    }
    pNew = uiCum / arma::sum(uiCum);
    // Zero-out probabilities that are negligible
    for(arma::uword j = 0; j < J; ++j){
      if(pNew(j) < zerotol){
        pNew(j) = 0;
      }
    }
    pNew = pNew / arma::sum(pNew);
    niter ++;
    // Rcpp::Rcout << "Iteration "<< niter << ": "<< pNew << std::endl;
  }
  return Rcpp::List::create(Rcpp::Named("p") = pNew,
                            Rcpp::Named("conv") = convergence,
                            Rcpp::Named("abstol") = abstol,
                            Rcpp::Named("niter") = niter);
}
