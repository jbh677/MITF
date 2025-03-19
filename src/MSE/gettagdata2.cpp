//////////////////////////////////////////////////////
// (partially) cpp'd get.tag.data2 function //////////
//////////////////////////////////////////////////////
// R. Hillary & P. Bessell-Browne CSIRO (2025) ///////
//////////////////////////////////////////////////////

#include <iostream>
#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
#include <mdarrclass.h>
#include "assert.h"

using namespace Rcpp;
using namespace std;

//[[Rcpp::export]]
RcppExport SEXP get_tag_data2(SEXP dm_,SEXP Trel_,SEXP Tmin_,SEXP Tx_,SEXP hx_,SEXP xi_,SEXP prep_,SEXP pret_,SEXP phiOD_,SEXP M_)
{

  IntegerVector dm = as<Rcpp::IntegerVector>(dm_); 
  int nr = dm[0];
  int nrec = dm[1];
  int rrel = dm[2];
  DoubleVector Tv = as<Rcpp::DoubleVector>(Tx_);
  DoubleVector hv = as<Rcpp::DoubleVector>(hx_); 
  DoubleVector pret = as<Rcpp::DoubleVector>(pret_); 
  double xi = as<double>(xi_);
  double prep = as<double>(prep_);
  double phiOD = as<double>(phiOD_);
  double Trel = as<double>(Trel_);
  double Tmin = as<double>(Tmin_);
  double M = as<double>(M_);
  
  arr2d NT;
  arr2d RT;
  arr2d hx;
  arr2d Tx;
  
  NT.arr2(nrec+1,nr);
  RT.arr2(nrec+1,nr); 
  hx.arr2(nrec+1,nr);
  Tx.arr2(nr,nr);

  int elem,r,rr,y; 
  double alpt,bett,ntmp,ptmp,omegaOD,prec;
  DoubleVector Trem(nr);

  for(r=0;r<nr;r++) {
    for(y=0;y<(nrec+1);y++) {

      elem = (nrec+1)*r+y;
      hx(y,r) = hv(elem);

    }
  }

  for(r=0;r<nr;r++) {
    for(rr=0;rr<nr;rr++) {

      elem = nr*r+rr;
      Tx(r,rr) = Tv(elem);
       
    }
  }

  // initialisation

  //omegaOD = (phiOD-1.)/(Trel-phiOD); 
  for(r=0;r<nr;r++) {

    NT(0,r) = r == rrel ? Trel : 0.;
    ptmp = xi*prep*hx(0,r);
    alpt = ((Trel-phiOD)*ptmp)/((1.-ptmp)*(ptmp+(1.-ptmp)*(phiOD-1.)));
    bett = (Trel-phiOD)/(ptmp+(1.-ptmp)*(phiOD-1.)); 
    prec = R::rbeta(alpt,bett);
    RT(0,r) = R::rbinom(NT(0,r),prec);

  }

  // loop across future recaptures

  for(y=1;y<(nrec+1);y++) {
    for(rr=0;rr<nr;rr++) {

      NT(y,rr) = 0.; 
      for(r=0;r<nr;r++) {

        ntmp = NT(y-1,r)-RT(y-1,r)/prep;
        Trem(r) = ntmp > 0. ? ntmp : 0.;
        NT(y,rr) += Trem(r)*pret(y)*exp(-M)*Tx(r,rr);

      }

      ntmp = double(round(NT(y,rr)));
      if(ntmp > Tmin) { 

        ptmp = hx(y,rr)*prep;
        //omegaOD = (phiOD-1.)/(ntmp-phiOD); 
        alpt = ((ntmp-phiOD)*ptmp)/((1.-ptmp)*(ptmp+(1.-ptmp)*(phiOD-1.)));
        bett = (ntmp-phiOD)/(ptmp+(1.-ptmp)*(phiOD-1.)); 
        prec = R::rbeta(alpt,bett); 
        RT(y,rr) = R::rbinom(ntmp,prec); 

      } else {
        
        RT(y,rr) = 0.;

      }
    }
  }

  DoubleVector tv((nrec+1)*nr);
  DoubleVector rv((nrec+1)*nr); 
  for(r=0;r<nr;r++) {
    for(y=0;y<(nrec+1);y++) {

      elem = (nrec+1)*r+y; 
      tv(elem) = NT(y,r);
      rv(elem) = RT(y,r);

    }
  }
  
  List res = Rcpp::List::create(Named("T")=tv,Named("R")=rv);

  NT.del();
  RT.del();
  hx.del();
  Tx.del();

  return Rcpp::wrap(res); 

}
