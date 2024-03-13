# MSE utility code

logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

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
  hmult <- Fleet$hmult
  N <- array(dim=c(ny,na,nr,2,nits))
  N[1,,,,1:nits] <- Biol$Ninit
  SSB <- array(dim=c(ny,nr,nits))
  SSBtot <- array(dim=c(ny,nits))
  H <- array(dim=c(ny,na,nr,2,nits))
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

  # initial harvest rates

  if(CatchControl) {

    C[1,,1:nits] <- Fleet$Cfix

  } else {

    hinit <- Fleet$hinit
    for(a in 1:na) {
      for(s in 1:2) {
       for(r in 1:nr) {

         ff <- grep(r,fref)
         H[1,a,r,s,1:nits] <- sum(hinit[ff]*sela[a,s,ff])

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
 
  }

  return(list(N=N,H=H,C=C,SSB=SSB,SSBtot=SSBtot))

} #}}}

#{{{ proj.para
proj.para <- function(Biol,Fleet,nyproj,CatchControl,nitsx,ncore) {

  wrap <- function(nitsx) {

    resx <- get.popdyn(Biol,Fleet,nyproj,CatchControl,nitsx) 
    return(resx)

  }

  restmp <- mclapply(rep(nitsx,ncore),wrap,mc.cores=ncore)

  # join everything together

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

} #}}}

#{{{ get.tag.data
get.tag.data <- function(xxx) {

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
  Ttot <- round(xxx$C*Tag$ntpt)[yrng,]

  # set up release/recapture data frame

  tdf <- expand.grid(yrel=1:nyt,arel=1:na,srel=1:2,rrel=1:nr,recev=1:nrec,recr=1:nr,nrel=NA,nrec=NA)

  for(y in yrng) {

    ymaxx <- min(ny,y+nrec)
    yx <- y:(y+nrec)
    nrecx <- length(yx)-1

    for(r in 1:nr) {

      # number of releases in this group

      tmpx <- rmultinom(1,Ttot[y,r],as.vector(paa))[,1]
      Nrel <- array(tmpx,dim=c(na,2))
      u0 <- prec <- array(dim=c(nrecx+1,nr)) 

      for(s in 1:2) {
        for(a in 1:na) {
          
          Nx <- Nrel[a,s]

          # check on minimum Nx (imposed via OD factor)

          if(Nx > Tmin) { 

            # within season adjustment of release numbers

            prec[1,] <- xi*hx[y,a,,s]
            Rws <- rbinom(1,Nx,xi*hx[y,a,r,s]/prep)
            Nx <- max(Nx-Rws,0)

            # survival + movement + tag shedding

            u0[,] <- 0
            u0[1,r] <- 1
            u1 <- u0
            aref <-a
            for(yy in 2:(nrecx+1)) {
              for(rr in 1:nr) {
              
                Tx <- Biol$T[,rr]
                u1[yy,rr] <- sum(u0[yy-1,]*pret[yy]*exp(-M)*(1-hx[yx[yy]-1,aref,,s])*Tx)

              }

              prec[yy,] <- hx[yx[yy],aref,,s]*prep
              u0 <- u1
              aref <- min(na,aref+1)
            }

            px <- u0*prec # P(tag surviving with geq 1 tag amd recaptured in that region/time)
 
            # simulate from Dirichlet-multinomial

            omegaOD <- (Nx-phiOD)/(phiOD-1)
            ppx <- as.vector(px[-1,])
            pxsum <- sum(ppx)  
            if(pxsum > 0.9) ppx <- 0.9*ppx/pxsum # constraint on P(tag rec >= 1)
            ptmp <- c(ppx,1-pxsum)
            alpx <- ptmp*omegaOD
            psam <- rdirichlet(1,alpx)
            Rx <- rmultinom(1,Nx,psam[1,])[-(nrecx*nr+1)]
            tdf[tdf$yrel==y & tdf$arel==a & tdf$srel==s & tdf$rrel==r,'nrel'] <- Nx
            tdf[tdf$yrel==y & tdf$arel==a & tdf$srel==s & tdf$rrel==r,'nrec'] <- Rx 

          }
        }
      }
    }
  }

  return(subset(tdf,!is.na(nrel)))

}
#}}}

