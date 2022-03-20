// conditional age-at-length code

#include <TMB.hpp>

template <class Type>
Type square(Type x) {return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  // data & priors

  DATA_ARRAY(lf);
  DATA_ARRAY(alk);
  DATA_VECTOR(lref);
  DATA_VECTOR(aref);
  DATA_VECTOR(alref);
  DATA_ARRAY(A);
  DATA_IVECTOR(dms);
  DATA_VECTOR(priors);

  // parameter decs

  PARAMETER(lnk);
  PARAMETER(logl1);
  PARAMETER(logl2);
  PARAMETER(lncvl);
  PARAMETER_VECTOR(lnmua);
  PARAMETER_VECTOR(lnsda); 

  // useful integers and process variables

  int ny = dms[0];
  int nbins = dms[1];
  int nages = dms[2];
  int a,aa,l,y;

  Type k = exp(lnk);
  Type l1 = exp(logl1);
  Type l2 = exp(logl2);
  Type cvl = exp(lncvl);
  Type linf;
  Type t0;
  Type a1 = alref[0]; // lower reference age
  Type a2 = alref[1]; // upper reference age
  Type sdl;
  Type mul;
  Type eps;
  Type psi;
  Type ldel;
  Type omega;
  Type psum;
  Type nll_alk;
  Type nll_lf;
  Type prior;
  Type objf = 0.;
  matrix<Type> pla(nbins,nages);
  matrix<Type> pa(ny,nages); 
  matrix<Type> pl(ny,nbins);
  array<Type> pal(ny,nbins,nages);
  array<Type> palt(ny,nbins,nages); 

  ///////////////////////////////////////////////////////////////////////
  // generate the probability distros: p(l | a), p(a), p(a | l) & p(l) //
  ///////////////////////////////////////////////////////////////////////

  sdl = sqrt(log(1.+square(cvl)));

  // Schnute parameterisation

  ldel = l2-l1;
  omega = 1.-exp(-k*(a2-a1));
  psi = (l1*omega)/ldel;
  linf = l1+ldel/omega;
  t0 = a1-log(psi+1.)/k; 

  // p(l | a)

  for(a=0;a<nages;a++) {
    for(psum=0.,l=0;l<nbins;l++) {
      mul = linf*(1.-exp(-k*(aref(a)-t0)));
      eps = log(mul/lref[l]);
      pla(l,a) = exp(-square(eps)/(2.*square(sdl)))/(sdl*sqrt(2.*M_PI));
      psum += pla(l,a);
	  }

	  for(l=0;l<nbins;l++) pla(l,a) /= psum;
  } 
   
  // p(a)

  
  for(y=0;y<ny;y++) {
  	for(psum=0.,a=0;a<nages;a++) {
	
      eps = log(aref(a))-lnmua[y];
      sdl = exp(lnsda[y]);
      pa(y,a) = exp(-square(eps)/(2.*square(sdl)))/(sdl*sqrt(2.*M_PI));
      psum += pa(y,a);
  	}

	 for(a=0;a<nages;a++) pa(y,a) /= psum;
  }

  // "true" p(a | l) and p(l)

  for(y=0;y<ny;y++) {
  	for(l=0;l<nbins;l++) {
		  for(psum=0.,a=0;a<nages;a++) {
		
			  palt(y,l,a) = pla(l,a) * pa(y,a);
			  psum += palt(y,l,a);
		
		  }
		  pl(y,l) = psum;
		  for(a=0;a<nages;a++) palt(y,l,a) /= psum;
    }
	  for(psum=0.,l=0;l<nbins;l++) psum += pl(y,l);
	  for(l=0;l<nbins;l++) pl(y,l) /= psum;
  } 
  
  // p(a | l) adjusted for ageing error

  for(y=0;y<ny;y++) {
  	for(l=0;l<nbins;l++) {
      for(a=0;a<nages;a++) {
		
        pal(y,l,a) = 0.;
			  for(aa=0;aa<nages;aa++) pal(y,l,a) += A(a,aa) * palt(y,l,aa);
		
      }
    }
  } 
  
  ///////////////////////////////////////////////////
  // calculate the likelihood, prior and posterior //
  ///////////////////////////////////////////////////

  for(nll_alk=0.,y=0;y<ny;y++) {
  	for(l=0;l<nbins;l++) {
			for(a=0;a<nages;a++) {
			
				if(alk(y,l,a) > 0) nll_alk -= alk(y,l,a)*log(pal(y,l,a));
			}
	  }
  }

  // length frequencies

  for(nll_lf=0.,y=0;y<ny;y++) {
  	for(l=0;l<nbins;l++) { 
      if(lf(y,l) > 0) nll_lf -= lf(y,l)*log(pl(y,l));
	  }
  }

  // prior on mean and SD in age distro

  for(prior=0.,y=0;y<ny;y++) {

	  eps = lnmua[y]-priors[0];
	  prior += log(priors[1])+square(eps)/(2.*square(priors[1]));
	  eps = lnsda[y]-priors[2];
	  prior += log(priors[3])+square(eps)/(2.*square(priors[3]));

  }

  objf = nll_alk + nll_lf + prior; 

  // ADREPORT

  ADREPORT(linf);
  ADREPORT(k);
  ADREPORT(t0);
  ADREPORT(l1);
  ADREPORT(l2);

  // REPORT stuff
  
  REPORT(cvl);
  REPORT(pla);
  REPORT(pal);
  REPORT(pa);
  REPORT(pl);
  REPORT(nll_alk);
  REPORT(nll_lf);
  REPORT(prior);

  return(objf);

}
