//////////////////////////////////////////////////////
// rdirich: Dirichlet random number generator ////////
//////////////////////////////////////////////////////
// R. Hillary (2025) /////////////////////////////////
//////////////////////////////////////////////////////

#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include "assert.h"

using namespace Rcpp;

//[[Rcpp::export]]
RcppExport SEXP rdirichlet(SEXP dm_,SEXP alp_)
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_);
  NumericVector alp = as<Rcpp::NumericVector>(alp_);
  int nits = dm[0];
  int m = dm[1];
  NumericMatrix res(nits,m);
  NumericVector tmp(m);

  int i,j;
  double dummy;

  for(i=0;i<nits;i++) {

    for(j=0;j<m;j++) tmp(j) = R::rgamma(alp(j),1.);
    dummy = sum(tmp);
    res.row(i) = tmp/dummy;

  }

  return Rcpp::wrap(res);
}
