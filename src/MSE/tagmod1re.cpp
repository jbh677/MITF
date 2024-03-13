// tag estimator 1: no space

#include <TMB.hpp>

template <class Type>
Type square(Type x) {return x*x;}

template <class Type>
Type ilogit(Type x){return Type(1)/(Type(1)+exp(-x));}

template<class Type>
Type objective_function<Type>::operator() ()
{

  DATA_IVECTOR(dms);
  DATA_IVECTOR(relyr);
  DATA_IVECTOR(nrecev);
  DATA_MATRIX(R);

  int ny = dms[0];
  int nrelev = dms[1];
  int nrecmax = dms[2];

  // fixed variables

  DATA_SCALAR(M);    // natural mortality
  DATA_SCALAR(prep); // reporting rate
  DATA_VECTOR(pret); // 1 or 2 tag retention probability

  // parameters

  PARAMETER(muih);         // logit-scale mean harvest rate
  PARAMETER(lnsigmah);     // log-scale RE variance
  PARAMETER_VECTOR(epshy); // logit-scale harvest rate REs

  vector<Type> hy(ny);
  vector<Type> ps(nrecmax);
  vector<Type> pr(nrecmax); 
  array<Type> Rhat(nrelev,nrecmax); 
  array<Type> phat(nrelev,nrecmax); 
  array<Type> pshat(nrelev,nrecmax);
  Type sigmah = exp(lnsigmah);

  int r1,r2,y;
  Type T,NR,htmp,ptmp,psum,pnotag;
  Type nll,pen,objf;

  for(y=0;y<ny;y++) hy(y) = ilogit(muih+epshy(y));

  // loop over release events
  
  nll = pen = Type(0);
  for(r1=0;r1<nrelev;r1++) {

    T = R(r1,0);
    NR = R(r1,1);
    y = relyr(r1);
    ps(0) = pret(0)*exp(-M)*(Type(1)-hy(y));
    pr(0) = hy(y+1)*prep*ps(0);
    psum = pr(0);
    for(r2=1;r2<nrecev(r1);r2++) {
    
      y = relyr(r1)+r2;
      ps(r2) = ps(r2-1)*pret(r2)*exp(-M)*(Type(1)-hy(y));
      pr(r2) = hy(y+1)*prep*ps(r2);
      psum += pr(r2);

    }

    pnotag = Type(1)-psum;

    // likelihood

    for(r2=0;r2<nrecev(r1);r2++) {

      Rhat(r1,r2) = T*pr(r2);
      phat(r1,r2) = pr(r2);
      pshat(r1,r2) = ps(r2);
      nll -= R(r1,r2+2)*log(pr(r2)+1e-8);

    }

    nll -= (T-NR)*pnotag;

    // penalties on overall recapture probability

  }

  // prior on year effects

  nll -= dnorm(epshy,Type(0),sigmah,true).sum(); 

  objf = nll+pen;

  REPORT(hy);
  REPORT(Rhat);
  REPORT(phat);
  REPORT(pshat);
  REPORT(sigmah);

  return(objf);

}
