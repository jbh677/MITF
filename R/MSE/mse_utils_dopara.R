# load shared objects for estimation models

dyn.load("tgm1.so")
dyn.load("tgm1re.so")
dyn.load("tgm2.so")
dyn.load("tgm2re.so")

# MSE utility code

logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

qt <- function(tau,psi,nu) {return(psi*exp(-nu*tau))}

ptag <- function(ntag,tau,psi,nu) {
  
  if(ntag == 2) return(qt(tau,psi,nu)^2)
  if(ntag == 1) return(2*qt(tau,psi,nu)*(1-qt(tau,psi,nu)))
  if(ntag == 0) return((1-qt(tau,psi,nu))^2)
  
}

p1tag <- function(ntag,tau,psi,nu) {

  if(ntag == 2) return(qt(tau,psi,nu)^2+2*qt(tau,psi,nu)*(1-qt(tau,psi,nu)))
  if(ntag == 1) return(qt(tau,psi,nu))

}

## sum to 1 utils code

# given n * (n-1) movement parameter matrix 
# get n * n transition matrix T where sum(rows(T)) = 1

add.logit <- function(x) {

  nr <- nrow(x)
  if(nr > 2) {

    .p1 <- exp(x)
    .p2 <- .p1/(1+rowSums(.p1))
    .p3 <- 1-rowSums(.p2)
    res <- unname(cbind(.p2,.p3))

  } else {
  
    res <- matrix(nrow=nr,ncol=nr)
    for(i in 1:nr) {
    
      rv <- exp(x[i,1])
      res[i,1] <- rv/(1+rv)
      res[i,2] <- 1-res[i,1]

    }
  }

  return(res)

}

# inverse bit: go from T to n * (n-1) parameters

inv.add.logit <- function(x) {

  nr <- nrow(x)

  if(nr > 2) {

    res <- matrix(nrow=nr,ncol=nr-1)
    for(i in 1:nr) {

      u <- T[i,-nr]
      dx <- length(u)
      A <- matrix(nrow=dx,ncol=dx)
      for(j in 1:dx)
        for(k in 1:dx) A[j,k] <- ifelse(j == k,1-u[j],-u[j])

      v <- as.vector(t(solve(A) %*% u))
      res[i,] <- log(v)
    }
  } else {

    res <- matrix(nrow=nr,ncol=nr-1)
    for(i in 1:nr) res[i,1] <- log(x[i,1]/(1-x[i,1]))

  }

  return(res)

}

### {{{ get.matwt
get.matwt <- function(Biol) {

  mula <- mata <- wta <- array(dim=c(Biol$dms[1],2))

  for(a in 1:Biol$dms[1]) 
    for(s in 1:2) mula[a,s] <- Biol$gpars$linf[s]*(1-exp(-Biol$gpars$k[s]*(Biol$ages[a]-Biol$gpars$t0[s])))

  # maturity and weight-at-age

  for(s in 1:2) {

    sdx <- Biol$gpars$sdla[s]
    for(a in 1:Biol$dms[1]) {
      
      lmin <- max(0,mula[a,s]*(1-sdx*1.96))
      lmax <- mula[a,s]*(1+sdx*1.96)
      lref <- seq(lmin, lmax,length=20)
      dl <- dlnorm(lref,log(mula[a,s]),sdx)
      dl <- dl/sum(dl)
      mlref <- 1/(1+19^(-(lref-Biol$mpars[1])/(Biol$mpars[2]-Biol$mpars[1])))
      wlref <- Biol$wpars[1]*lref^Biol$wpars[2]
      mata[a,s] <- sum(mlref *dl)
      wta[a,s] <- sum(wlref*dl)

    }
  }

  Biol$mula <- mula
  Biol$mata <- mata
  Biol$wta <- wta

  return(Biol)

} #}}}

#{{{ get.pla

get.pla <- function(Biol) {

  nl <- length(Biol$mulbins)
  na <- Biol$dms[1]
  pla <- array(dim=c(nl,na,2))
  for(s in 1:2) {
    for(a in 1:na) {
    
      dx <- dlnorm(Biol$mulbins,log(Biol$mula[a,s]),Biol$gpars$sdla[s],FALSE)
      dx <- dx/sum(dx)
      pla[,a,s] <- dx

    }
  }

  Biol$pla <- pla

  return(Biol)

} #}}}

#{{{ get.tag.age.distro

get.tag.age.distro <- function(Biol,Tag) {

  sdx <- Biol$gpars$sdla
  mula <- Biol$mula
  mulbins <- Biol$mulbins
  mulrel <- Tag$mulrel
  sdlrel <- Tag$sdlrel
  dlrel <- dlnorm(mulbins,mulrel,sdlrel,FALSE)
  dlrel <- dlrel/sum(dlrel)
  paa <- array(dim=c(Biol$dms[1],2))
  for(s in 1:2)
    for(a in 1:Biol$dms[1]) paa[a,s] <- sum(Biol$pla[,a,s]*dlrel)

  return(paa[]/sum(paa))

} #}}}

#{{{ get.sel
get.sel <- function(Fleet,Biol) {

  mula <- Biol$mula
  nff <- Fleet$nf
  sela <- array(dim=c(Biol$dms[1],2,nff))
  sref <- rep(NA,20)

  for(f in 1:nff) {

    selp <- Fleet$selpars[,f]
    for(s in 1:2) {

      sdx <- Biol$gpars$sdla[s] 
      for(a in 1:Biol$dms[1]) {
        
        lmin <- max(0,mula[a,s]*(1-sdx*1.96))
        lmax <- mula[a,s]*(1+sdx*1.96)
        lref <- seq(lmin, lmax,length=20)
        dl <- dlnorm(lref,log(mula[a,s]),sdx)
        dl <- dl/sum(dl) 
        for(i in 1:20) {
          
          if(lref[i] < selp[1]) sref[i] <- 2^{-(lref[i]-selp[1])^2/(selp[2]^2)}
          if(lref[i] >= selp[1]) sref[i] <- 2^{-(lref[i]-selp[1])^2/(selp[3]^2)} 

        }

        sela[a,s,f] <- sum(sref*dl)
      }
    }
  }
  
  Fleet$sela <- sela

  return(Fleet)
} #}}}

#{{{ get.unexploited.eqm
get.unexploited.eqm <- function(Biol) {

  na <- Biol$dms[1]
  nr <- Biol$dms[2]
  M <- Biol$M
  R0 <- Biol$R0
  mat <- Biol$mata[,1]
  wt <- Biol$wta*1e-3
  wtsp <- wt[,1]
  eta <- Biol$eta
  T <- Biol$T
  N <- array(dim=c(na,nr,2))
  for(r in 1:nr) N[1,r,] <- c(0.5,0.5)*eta[r]
  for(s in 1:2) 
    for(a in 2:(na-1)) 
      for(r in 1:nr) N[a,r,s] <- sum(T[,r]*N[a-1,,s]*exp(-M))

  Meqm <- matrix(nrow=nr,ncol=nr)
  veqm <- matrix(nrow=nr,ncol=1)
  for(s in 1:2) {
    for(r in 1:nr) {

      veqm[r,1] <- 0
      for(rr in 1:nr) {

        veqm[r,1] <- veqm[r,1]+T[rr,r]*N[na-1,rr,s]*exp(-M)
        if(r == rr) Meqm[rr,r] <- 1-T[rr,r]*exp(-M)
        if(r != rr) Meqm[rr,r] <- -T[rr,r]*exp(-M)

      }
    }

    Minv <- solve(Meqm)
    veqmt <- t(veqm)
    npgeqm <- veqmt %*% Minv
    N[na,,s] <- npgeqm[1,]
  }

  rho <- sum(apply(N[,,1],2,function(x,wtsp,mat){x <- sum(x*wtsp*mat)},wtsp,mat))

  res <- list(N=N[]*R0,rho=rho)
} #}}}

#{{{ get.exploited.eqm
get.exploited.eqm <- function(Biol,Fleet,hvec) {

  resu <- get.unexploited.eqm(Biol)
  rho0 <- resu$rho

  na <- Biol$dms[1]
  nr <- Biol$dms[2]
  M <- Biol$M
  R0 <- Biol$R0
  hh <- Biol$hh
  mat <- Biol$mata[,1]
  wt <- Biol$wta*1e-3
  wtsp <- wt[,1]
  eta <- Biol$eta
  T <- Biol$T
  N <- H <- array(dim=c(na,nr,2)) 

  # set up H[a,r,s]

  sel <- Fleet$sela
  fref <- Fleet$fref
  for(a in 1:na) {
    for(s in 1:2) {
      for(r in 1:nr) {

        ff <- grep(r,fref)
        H[a,r,s] <- sum(hvec[ff]*sel[a,s,ff])

      }
    }
  }
  
  for(r in 1:nr) N[1,r,] <- c(0.5,0.5)*eta[r]
  for(s in 1:2) 
    for(a in 2:(na-1)) 
      for(r in 1:nr) N[a,r,s] <- sum(T[,r]*N[a-1,,s]*exp(-M)*(1-H[a-1,,s])) 

  Meqm <- matrix(nrow=nr,ncol=nr)
  veqm <- matrix(nrow=nr,ncol=1)
  for(s in 1:2) {
    for(r in 1:nr) {

      veqm[r,1] <- 0
      for(rr in 1:nr) {

        veqm[r,1] <- veqm[r,1]+T[rr,r]*N[na-1,rr,s]*exp(-M)*(1-H[na-1,rr,s])
        if(r == rr) Meqm[rr,r] <- 1-T[rr,r]*exp(-M)*(1-H[na-1,rr,s])
        if(r != rr) Meqm[rr,r] <- -T[rr,r]*exp(-M)*(1-H[na-1,rr,s])

      }
    }

    Minv <- solve(Meqm)
    veqmt <- t(veqm)
    npgeqm <- veqmt %*% Minv
    N[na,,s] <- npgeqm[1,]
  }

  rho <- sum(apply(N[,,1],2,function(x,wtsp,mat){x <- sum(x*wtsp*mat)},wtsp,mat))
  
  # get relative catch-by-fishery

  nff <- Fleet$nf
  cf <- rep(0,nff)
  for(f in 1:nff) {
    for(s in 1:2) {

      cf[f] <- sum(N[,fref[f],s]*sel[,s,f]*wt[a,s]*hvec[f])

    }
  }

  # SSB targets

  sprat <- rho/rho0
  dep <- max(0,(4*hh*sprat+hh-1)/(5*hh-1))
  Rinit <- R0*(4*hh*dep)/(hh*(5*dep-1)+1-dep)
 
  res <- list(N=N[]*Rinit,rho=rho,dep=dep,Rinit=Rinit,C=cf) 
} #}}}


#{{{ objfn.init
objfn.init <- function(theta) {

  hvec <- 1/(1+exp(-theta))
  resx <- get.exploited.eqm(Biol,Fleet,hvec)
  targhat <- logit(c(resx$dep,resx$C/sum(resx$C)))

  return(sum((targhat-targv)^2))

} #}}}

#{{{ get.popdyn
get.popdyn <- function(Biol,Fleet,nyproj=20,CatchControl=FALSE,nits=2) {

  na <- Biol$dms[1]
  nr <- Biol$dms[2]
  ny <- nyproj
  nf <- Fleet$nf
  fref <- Fleet$fref
  sela <- Fleet$sela
  hmult <- Fleet$hmul
  N <- array(dim=c(ny,na,nr,2,nits))
  N[1,,,,1:nits] <- Biol$Ninit[]
  SSB <- array(dim=c(ny,nr,nits))
  SSBtot <- array(dim=c(ny,nits))
  H <- array(0,dim=c(ny,na,nr,2,nits))
  C <- array(dim=c(ny,nf,nits))
  mat <- Biol$mata[,1]
  wt <- Biol$wta*1e-3
  wtsp <- wt[,1]
  eta <- Biol$eta
  M <- Biol$M
  T <- Biol$T
  alp <- Biol$alp
  bet <- Biol$bet
  sigmar <- Biol$sigmar

  # initial SSB

  SSB[1,,] <- apply(N[1,,,1,],c(2,3),function(x,wtsp,mat){x <- sum(x*wtsp*mat)},wtsp,mat)
  SSBtot[1,] <- apply(SSB[1,,],2,sum)
  #SSB[1,,] <- apply(N[1,,,1,],2,function(x,wtsp,mat){x <- sum(x*wtsp*mat)},wtsp,mat)
  #SSBtot[1,] <- sum(SSB[1,,])

  # initial harvest rates

  if(CatchControl) {

    C[1,,] <- Fleet$Cfix
    for(f in 1:nf) {

      X <- apply(N[1,,fref[f],1,]*wt[,1]*sela[,1,f]+N[1,,fref[f],2,]*wt[,2]*sela[,2,f],2,sum)
      htmp <- C[1,f,]/X
      if(any(htmp > 0.9)) {

        htmp[htmp > 0.9] <- 0.9
        C[1,f,] <- X*htmp

      }

      for(s in 1:2) H[1,,fref[f],s,] <- H[1,,fref[f],s,]+sela[,s,f] %o% htmp 

    }

  } else {

    hinit <- Fleet$hinit

    for(a in 1:na) {
      for(s in 1:2) {
       for(r in 1:nr) {

         ff <- grep(r,fref)
         if(is.null(dim(hinit))) {

           H[1,a,r,s,1:nits] <- sum(hinit[ff]*sela[a,s,ff])

         } else {

           H[1,a,r,s,] <- hinit[ff,]*sela[a,s,ff]

         }
        }
      }
    }

    for(f in 1:nf) {

      rx <- fref[f]
      C[1,f,] <- 0
      for(s in 1:2) {

        wts <- wt[,s]
        C[1,f,] <- C[1,f,]+apply(N[1,,rx,s,]*H[1,,rx,s,],2,function(x,wts){x <- sum(x*wts)},wts)

      }
    }
  }

  # loop through the years

  for(y in 2:ny) {

    Rtot <- alp*SSBtot[y-1,]/(1+bet*SSBtot[y-1,])
    for(r in 1:nr) N[y,1,r,,] <- (c(0.5,0.5)*eta[r]) %o% (Rtot*rlnorm(nits,-sigmar^2/2,sigmar))
    for(s in 1:2) {
      for(a in 2:(na-1)) { 
        for(r in 1:nr) {
          
          Tx <- T[,r]
          N[y,a,r,s,] <- apply(N[y-1,a-1,,s,]*(1-H[y-1,a-1,,s,])*exp(-M),2,function(x,Tx){x <- sum(x*Tx)},Tx)

        }
      }
    }

    # plus group

    for(s in 1:2) {
      for(r in 1:nr) {
    
        Tx <- T[,r]
        N[y,na,r,s,] <- apply(N[y-1,na-1,,s,]*(1-H[y-1,na-1,,s,])*exp(-M),2,function(x,Tx){x <- sum(x*Tx)},Tx)+apply(N[y-1,na,,s,]*(1-H[y-1,na,,s,])*exp(-M),2,function(x,Tx){x <- sum(x*Tx)},Tx)

      }
    }

    # harvest rates

    if(CatchControl) {

      C[y,,] <- Fleet$Cfix
      for(f in 1:nf) {

        X <- apply(N[y,,fref[f],1,]*wt[,1]*sela[,1,f]+N[y,,fref[f],2,]*wt[,2]*sela[,2,f],2,sum)
        htmp <- C[y,f,]/X
        if(any(htmp > 0.9)) {

          htmp[htmp > 0.9] <- 0.9
          C[y,f,] <- X*htmp

        }

        for(s in 1:2) H[y,,fref[f],s,] <- H[y,,fref[f],s,]+sela[,s,f] %o% htmp 

      } 

    } else { 

      H[y,,,,] <- hmult[y]*H[y-1,,,,]

      for(f in 1:nf) {

        rx <- fref[f]
        C[y,f,] <- 0
        for(s in 1:2) {

          wts <- wt[,s]
          C[y,f,] <- C[y,f,]+apply(N[y,,rx,s,]*H[y,,rx,s,],2,function(x,wts){x <- sum(x*wts)},wts)

        }
      } 
    }

    # SSB

    SSB[y,,] <- apply(N[y,,,1,],c(2,3),function(x,wtsp,mat){x <- sum(x*wtsp*mat)},wtsp,mat)
    SSBtot[y,] <- apply(SSB[y,,],2,sum)

    # SSB[y,,] <- apply(N[y,,,1,],2,function(x,wtsp,mat){x <- sum(x*wtsp*mat)},wtsp,mat)
    # SSBtot[y,] <- sum(SSB[y,,])
 
  }

  return(list(N=N,H=H,C=C,SSB=SSB,SSBtot=SSBtot))

} #}}}

#{{{ proj.para
proj.para <- function(Biol,Fleet,nyproj,CatchControl,nitsx,ncore,link=TRUE) {

  wrap <- function(nitsx) {

    resx <- get.popdyn(Biol,Fleet,nyproj,CatchControl,nitsx) 
    return(resx)

  }

  nx <- rep(nitsx,ncore)
  restmp <- foreach(i=1:ncore) %dopar% wrap(nx[i])

  # join everything together

  if(link) {

    resz <- list()
    resz[['N']] <- restmp[[1]]$N
    resz[['H']] <- restmp[[1]]$H
    resz[['C']] <- restmp[[1]]$C
    resz[['SSB']] <- restmp[[1]]$SSB
    resz[['SSBtot']] <- restmp[[1]]$SSBtot

    for(i in 2:length(restmp)) {

      resz[['N']] <- abind(resz$N,restmp[[i]]$N,along=5)
      resz[['H']] <- abind(resz$H,restmp[[i]]$H,along=5) 
      resz[['C']] <- abind(resz$C,restmp[[i]]$C,along=3) 
      resz[['SSB']] <- abind(resz$SSB,restmp[[i]]$SSB,along=3) 
      resz[['SSBtot']] <- abind(resz$SSBtot,restmp[[i]]$SSBtot,along=2) 
    
    }

    return(resz)

  } else {

    return(restmp)

  }

} #}}}

#{{{ proj.list - when in MP switched on projection mode
proj.list <- function(prjobj,nyp,X) {

  # X: management quantity i.e. hmult or cfix in the Fleet object 

  wrapp <- function(nchain) {

    Fleett <- Fleet

    iterx <- 
    if(CatchControl) {
      
      Fleett$cfix <- X[iterx]

    } else {

      Fleett$hinit <- rep(X[iterx],nr)
      Fleett$hmult <- rep(1,nyp)

    }

    Bioll <- Biol
    dmn <- dim(prjobj[[nchain]]$N)
    Bioll$Ninit <- prjobj[[nchain]]$N[dmn[1],,,,] 
    resx <- get.popdyn(Bioll,Fleett,nyp,CatchControl,nitsx)

  }

  restmp <- foreach(i=1:length(prjobj)) %dopar% wrapp(i)
 

}

#{{{ get.tag.data2 - follow tagged population instead of Brownie sim approach

get.tag.data2 <- function(xxx) {

  ny <- dim(xxx$N)[1]
  na <- dim(xxx$N)[2]
  nr <- dim(xxx$N)[3]
  pla <- Biol$pla
  nl <- length(Biol$mulbins)
  nf <- Fleet$nf
  fref <- Fleet$fref
  paa <- get.tag.age.distro(Biol,Tag) 
  Ctot <- apply(xxx$C,1,sum)
  M <- Biol$M
  hx <- xxx$H
  nrel <- length(Tag$ytagoff:Tag$ytagon)
  nrec <- Tag$maxrecev
  pret <- p1tag(2,0:nrec,psi,nu)
  phiOD <- Tag$phi # over-dispersion
  Tmin <- floor(phiOD)+1
  xi <- Tag$xi # within-season recapture "discount" probability

  # release numbers by region

  yrng <- Tag$ytagon:Tag$ytagoff
  nyt <- length(yrng)
  Ttot <- round(xxx$C*Tag$ntpt)

  tdf <- expand.grid(yrel=yrng,arel=1:na,srel=1:2,rrel=1:nr,recev=1:(nrec+1),recr=1:nr,T=NA,R=NA) 

  NT <- RT <- HT <- array(0,dim=c(nrec+1,nr)) 
  for(y in yrng) {

    ymaxx <- min(ny,y+nrec)
    yx <- y:ymaxx
    nrecx <- length(yx)-1
 
    # loop starts with release region

    for(r in 1:nr) {

      tmpx <- rmultinom(1,Ttot[y,r],as.vector(paa))[,1]
      Nrel <- array(tmpx,dim=c(na,2)) 

      for(s in 1:2) {
        for(a in 1:na) {
          
          Nx <- Nrel[a,s]

          # check on minimum Nx (imposed via OD factor)

          NT[] <- RT[] <- 0
          if(Nx > Tmin) { 
                          
            NT[1,r] <- Nx
            ptmp <- xi*prep*hx[y,a,,s]
            omegaOD <- (phiOD-1)/(Nx-phiOD) 
            alpt <- ((Nx-phiOD)*ptmp)/((1-ptmp)*(ptmp+(1-ptmp)*(phiOD-1)))
            bett <- (Nx-phiOD)/(ptmp+(1-ptmp)*(phiOD-1))
            prec <- rbeta(nr,alpt,bett)
            HT[1,] <- prec
            RT[1,] <- rbinom(nr,NT[1,],prec)
            aref <- a
            
            # loop across future recaptures

            for(yy in 2:(nrecx+1)) {
              for(rr in 1:nr) {
              
                Tx <- Biol$T[,rr]
                Trem <- NT[yy-1,]-RT[yy-1,]/prep
                Trem[Trem < 0] <- 0
                NT[yy,rr] <- sum(Trem*pret[yy]*exp(-M)*Tx)
                Ttmp <- round(NT[yy,rr])
                if(Ttmp > Tmin) {
                   
                  ptmp <- hx[yx[yy],aref,rr,s]*prep 
                  omegaOD <- (phiOD-1)/(Tmin-phiOD) 
                  alpt <- ((Ttmp-phiOD)*ptmp)/((1-ptmp)*(ptmp+(1-ptmp)*(phiOD-1)))
                  bett <- (Ttmp-phiOD)/(ptmp+(1-ptmp)*(phiOD-1))
                  prec <- rbeta(1,alpt,bett)
                  HT[yy,rr] <- prec
                  RT[yy,rr] <- rbinom(1,Ttmp,prec) 
                }
              } 
              aref <- min(na,aref+1) 

            }

            tdf[tdf$yrel==y & tdf$arel==a & tdf$srel==s & tdf$rrel==r,'T'] <- as.vector(NT)
            tdf[tdf$yrel==y & tdf$arel==a & tdf$srel==s & tdf$rrel==r,'R'] <- as.vector(RT)
          }
        }
      }
    }
  }

  return(subset(tdf,!is.na(T) & (yrel+recev-1<=ny))) 

}
#}}}

# {{{ extend.tag.data - needed for projection mode
extend.tag.data <- function(df,xxx,ycurr,yprj) {

  # selection critieria for what to keep from previous data set:
  # in final year of previous simulations (ycurr) T > Tmin

  Tmin <- floor(Tag$phi)+2
  nrecmax <- Tag$maxrecev
  trem <- subset(df,(df$yrel+df$recev-1) == ycurr & recev < nrecmax+1) 
  xdf <- aggregate(trem$T,by=list(yrel=trem$yrel,arel=trem$arel,srel=trem$srel,rrel=trem$rrel,recev=trem$recev),FUN=sum)
  xdf <- subset(xdf,x>Tmin)
  hx <- xxx$H
  pret <- p1tag(2,0:nrecmax,psi,nu)
  phiOD <- Tag$phi # over-dispersion
  M <- Biol$M

  # calculate remaining recapture events given cap 

  ydel <- yprj-ycurr

  # build up extended data base entry-by-entry

  for(i in 1:dim(xdf)[1]) {

    y <- xdf[i,'yrel']
    a <- xdf[i,'arel']
    s <- xdf[i,'srel']
    r <- xdf[i,'rrel']
    xrecev <- xdf[i,'recev']
    xrecmax <- min(nrecmax+1,xrecev+ydel)
    Tinit <- subset(trem,yrel==y & arel==a & srel==s & rrel==r)$T
    Rinit <- subset(trem,yrel==y & arel==a & srel==s & rrel==r)$R
    tzdf <- expand.grid(yrel=y,arel=a,srel=s,rrel=r,recev=(xrecev+1):xrecmax,recr=1:nr,T=NA,R=NA) 
    nyx <- length(xrecev:xrecmax)
    NT <- RT <- HT <- array(0,dim=c(nyx,nr)) 
    NT[1,] <- Tinit
    RT[1,] <- Rinit
    yx <- ycurr:(ycurr+nyx-1)
    aref <- min(a+xrecev-1,na)
    for(yy in 2:nyx) {
      for(rr in 1:nr) {

        Tx <- Biol$T[,rr]
        Trem <- NT[yy-1,]-RT[yy-1,]/prep
        Trem[Trem < 0] <- 0
        NT[yy,rr] <- sum(Trem*pret[yy]*exp(-M)*Tx)
        Ttmp <- round(NT[yy,rr])
        if(Ttmp > (Tmin-1)) {
         
          ptmp <- hx[yx[yy],aref,rr,s]*prep 
          omegaOD <- (phiOD-1)/(Tmin-phiOD) 
          alpt <- ((Ttmp-phiOD)*ptmp)/((1-ptmp)*(ptmp+(1-ptmp)*(phiOD-1)))
          bett <- (Ttmp-phiOD)/(ptmp+(1-ptmp)*(phiOD-1))
          prec <- rbeta(1,alpt,bett)
          HT[yy,rr] <- prec
          RT[yy,rr] <- rbinom(1,Ttmp,prec) 
        }
      }

      aref <- min(na,aref+1) 

    }

    tzdf$T <- as.vector(NT[-1,])
    tzdf$R <- as.vector(RT[-1,]) 
    
    if(i == 1) {
    
      texdf <- tzdf

    } else {

      texdf <- rbind(texdf,tzdf)
  
    }

  }

  return(texdf)

}
# }}}

#{{{ get.hrates
get.hrates <- function(df=tdf,model='tgm1',mctrl) {

  # get number of repcature regions

  nr <- length(unique(df$rrel))

  if(model=='tgm1' | model=='tgm1re') {

    # aggregate the data for non-spatial models

    ttdf <- subset(df,recev>1) 
    ttdf <- aggregate(ttdf$R,by=list(recev=ttdf$recev,yrel=ttdf$yrel),FUN=sum)
    ntdf <- subset(df,recev==1)
    ntdf$x <- ntdf$T-ntdf$R
    ntdf <- aggregate(ntdf$x,by=list(recev=ntdf$recev,yrel=ntdf$yrel),FUN=sum)
    nrdf <- aggregate(ttdf$x,by=list(yrel=ttdf$yrel),FUN=sum)
    relyr <- ntdf$yrel
    nrel <- ntdf$x
    nrelev <- length(nrel)
    ntot <- nrdf$x
    nrecev <- aggregate(ttdf$recev,by=list(yrel=ttdf$yrel),FUN=max)$x-1
    nrecmax <- max(nrecev)

    R <- array(dim=c(nrelev,nrecmax+2))
    R[,1] <- nrel
    R[,2] <- ntot
    for(y in 1:nrelev) R[y,3:(2+nrecev[y])] <- subset(ttdf,yrel==relyr[y])$x

    nytot <- max(relyr+nrecev)-relyr[1]+1
    dms <- c(nytot,nrelev,max(nrecev))
    pret <- p1tag(2,1:nrecmax,psi,nu)
    relyr <- relyr-min(relyr)+1

    # set up data object

    if(model=='tgm1') {
    
      data <- list(dms=dms,
             relyr=relyr-1,
             nrecev=nrecev,
             R=R,
             M=mctrl$M,
             prep=rep(prep,nytot),
             pret=pret,
             sigmah=mctrl$sigmah)

      pars <- list(muih=logit(0.03),epshy=rep(0,nytot))

      obj1 <- MakeADFun(data=data,parameters=pars,DLL="tgm1",silent=TRUE)
      res1 <- do.call(optim,obj1)
      rep1 <- obj1$rep()  

    } else {

      data <- list(dms=dms,
             relyr=relyr-1,
             nrecev=nrecev,
             R=R,
             M=mctrl$M,
             prep=rep(prep,nytot),
             pret=pret)

      pars <- list(muih=logit(0.03),lnsigmah=log(mctrl$sigmah),epshy=rep(0,nytot)) 
      
      obj1 <- MakeADFun(data=data,parameters=pars,random=c("epshy"),DLL="tgm1re",silent=TRUE)
      res1 <- do.call(optim,obj1)
      rep1 <- obj1$rep()

    }

    return(rep1)

  }

  if(model=='tgm2' | model=='tgm2re') {

    # aggregate the data for spatial models

    tsdf <- subset(df,recev>1) 
    tsdf <- aggregate(tsdf$R,by=list(recev=tsdf$recev,yrel=tsdf$yrel,rrel=tsdf$rrel,recr=tsdf$recr),FUN=sum) 
    ntsdf <- subset(df,recev==1)
    ntsdf$x <- ntsdf$T-ntsdf$R
    ntsdf <- aggregate(ntsdf$x,by=list(recev=ntsdf$recev,yrel=ntsdf$yrel,rrel=ntsdf$rrel),FUN=sum)
    nrsdf <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,rrel=tsdf$rrel),FUN=sum)
    relsyr <- ntsdf$yrel
    nsrel <- ntsdf$x
    nsrelev <- length(nsrel)
    nstot <- nrsdf$x
    nsrecev <- aggregate(tsdf$recev,by=list(yrel=tsdf$yrel,rrel=tsdf$rrel),FUN=max)$x-1
    nsrecmax <- max(nsrecev) 
    nystot <- max(relsyr+nsrecev)-relsyr[1]+1
    psret <- p1tag(2,1:nsrecmax,psi,nu) 
    ysrel <- ntsdf$yrel
    nystot <- max(relsyr+nsrecev)-relsyr[1]+1
    dmss <- c(nystot,nr,nsrelev,nsrecmax)
    RS <- array(dim=c(nsrelev,nsrecmax,dmss[2]))
    for(t in 1:nsrelev) {
      for(r in 1:dmss[2]) {

        yx <- ysrel[t]
        rx <- nrsdf$rrel[t]
        nrecx <- nsrecev[t]
        xdf <- subset(tsdf,yrel==yx & rrel==rx)
        for(rr in 1:dmss[2]) RS[t,1:nrecx,rr] <- subset(xdf,recr==rr)$x
      }
    }
    ysrel <- ysrel-min(ysrel)+1
    
    if(model == 'tgm2') {

      data2 <- list(dms=dmss,
              relyr=ysrel-1,
              relr=nrsdf$rrel-1,
              nrecev=nsrecev,
              T=ntsdf$x,
              NR=nrsdf$x,
              R=RS,
              M=mctrl$M,
              prep=rep(prep,nystot),
              pret=psret,
              sigmah=mctrl$sigmah)

      xphi <- inv.add.logit(Biol$T)
      pars2 <- list(muih=logit(rep(0.03,nr)),
              epshy=matrix(0,nrow=nystot,ncol=nr),
              xphi=xphi)

      obj2 <- MakeADFun(data=data2,parameters=pars2,DLL="tgm2",silent=TRUE)
      res2 <- do.call(optim,obj2)
      rep2 <- obj2$rep()
 
    } else {

      data2 <- list(dms=dmss,
              relyr=ysrel-1,
              relr=nrsdf$rrel-1,
              nrecev=nsrecev,
              T=ntsdf$x,
              NR=nrsdf$x,
              R=RS,
              M=mctrl$M,
              prep=rep(prep,nystot),
              pret=psret) 

      xphi <- inv.add.logit(Biol$T)
      pars2re <- list(muih=logit(rep(0.03,nr)),
                lnsigmah=log(mctrl$sigmah),
                epshy=matrix(0,nrow=nystot,ncol=nr),
                xphi=xphi)

      obj2 <- MakeADFun(data=data2,parameters=pars2re,random=c("epshy"),DLL="tgm2re",silent=TRUE)
      res2 <- do.call(optim,obj2)
      rep2 <- obj2$rep() 

    }

    return(rep2)
  } 


}
#}}}

#{{{ listify.proj - facilitates parallel coding facets
listify.proj <- function(prj) {

  nitsx <- dim(prj$N)[5]
  plst <- list()
  for(nn in 1:nitsx) plst[[nn]] <- list(N=prj$N[,,,,nn],H=prj$H[,,,,nn],C=prj$C[,,nn])

  return(plst)

}
#}}}

#{{{ wrap.get.tag.data2
wrap.get.tag.data2 <- function(k) {

  xxx <- prjlist[[k]]
  tdf <- get.tag.data2(xxx)
  return(tdf)

}
#}}}

#{{{ wrap.extend.tag.data
wrap.extend.tag.data <- function(k) {

  txdf <- extend.tag.data(tdflist[[k]],prjlist[[k]],ycurr,yprj)
  return(rbind(tdflist[[k]],txdf))

} #}}}

#{{{ wrap.get.hrates
wrap.get.hrates <- function(k) {
  
  rep <- get.hrates(tdflist[[k]],model,mctrl)
  return(rep)

}
#}}}

#{{{ get.TAC
get.TAC <- function(tacctrl) {

  if(tacctrl$MP == 'fixedh') {

    ytau <- tacctrl$ytau
    htarg <- tacctrl$htarg
    delmax <- tacctrl$maxChange
    ymax <- dim(hyhat)[1]
    hbar <- apply(hyhat[(ymax-ytau+1):(ymax),],2,mean)
    delta <- htarg/hbar
    delta[delta < 1-delmax] <- 1-delmax
    delta[delta > 1+delmax] <- 1+delmax 
    TACnew <- TACold*delta
    return(TACnew)

  }

}#}}}
