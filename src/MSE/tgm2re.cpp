// tag estimator 2: no space + RE

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
  DATA_IVECTOR(relr);
  DATA_IVECTOR(nrecev);
  DATA_VECTOR(T);
  DATA_VECTOR(NR);
  DATA_ARRAY(R);

  int ny = dms[0];
  int nr = dms[1];
  int nrelev = dms[2];
  int nrecmax = dms[3];

  // fixed variables

  DATA_SCALAR(M);    // natural mortality
  DATA_VECTOR(prep); // reporting rate
  DATA_VECTOR(pret); // 1 or 2 tag retention probability

  // parameters

  PARAMETER_VECTOR(muih);  // logit-scale mean harvest rate
  PARAMETER(lnsigmah); // RE variance
  PARAMETER_MATRIX(epshy); // logit-scale harvest rate REs
  PARAMETER_MATRIX(xphi);  // movement parameters  

  array<Type> hy(ny,nr);
  array<Type> pr(nrecmax,nr); 
  array<Type> u0(nrecmax+1,nr);
  array<Type> u1(nrecmax+1,nr); 
  matrix<Type> Phi(nr,nr);
  vector<Type> rv(nr-1);
  Type sigmah = exp(lnsigmah);

  int r,rr,t,t2,y;
  Type htmp,ptmp,psum,pnotag,rsum;
  Type nll,pen,objf;

  for(r=0;r<nr;r++) 
    for(y=0;y<ny;y++) hy(y,r) = ilogit(muih(r)+epshy(y,r));

  if(nr == 2) {

    for(r=0;r<nr;r++) {

      rv(0) = exp(xphi(r,0));
      Phi(r,0) = rv(0)/(Type(1)+rv(0));
      Phi(r,1) = Type(1)-rv(0)/(Type(1)+rv(0));

    }

  } else {

    for(r=0;r<nr;r++) {
      
      rv = xphi.row(r);
      rv = exp(rv);
      rsum = rv.sum();
      for(rr=0;rr<nr-1;rr++) Phi(r,rr) = rv(rr)/(Type(1)+rsum);
      rsum = Type(0);
      for(rr=0;rr<nr-1;rr++) rsum += Phi(r,rr); 
      Phi(r,nr-1) = Type(1)-rsum;

    }
  }

  // loop over release events
  
  nll = pen = Type(0);
  for(t=0;t<nrelev;t++) {

      u0.setZero();
      u1.setZero();
      u0(0,relr(t)) = Type(1); // start it in release region
      y = relyr(t);

      // loop across recapture events and regions

      psum = Type(0); 
      for(t2=1;t2<=nrecev(t);t2++) {
  
        y = relyr(t)+t2; 

        for(r=0;r<nr;r++) 
          for(rr=0;rr<nr;rr++) u1(t2,r) += u0(t2-1,rr)*Phi(rr,r)*pret(t2-1)*exp(-M)*(Type(1)-hy(y-1,rr));
          
        for(r=0;r<nr;r++) {

          u0(t2,r) = u1(t2,r);
          pr(t2-1,r) = u0(t2,r)*prep(y)*hy(y,r);
          psum += pr(t2-1,r);

        }
      }

      pnotag = Type(1)-psum;
      
      // likelihood

      for(t2=0;t2<nrecev(t);t2++) 
        for(r=0;r<nr;r++) nll -= R(t,t2,r)*log(pr(t2,r)+1e-8);

      nll -= (T(t)-NR(t))*pnotag;
  }

  // penalties on overall recapture probability

  // prior on year effects

  for(r=0;r<nr;r++)
    for(y=0;y<ny;y++) nll -= dnorm(epshy(y,r),Type(0),sigmah,true); 

  objf = nll+pen;

  REPORT(hy);
  REPORT(Phi);
  REPORT(sigmah);

  return(objf);

}
