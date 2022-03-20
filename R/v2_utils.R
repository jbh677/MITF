######################################################
# general purpose fitting plots for Macca ############
######################################################
# R Hillary CSIRO 2019 ###############################
######################################################

library(ggplot2)

# useful functions

logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

# plot length data
# args: report object, fishery, years, length bins

plot.lenfreq <- function(lfdat,lfcov,rep,ff,phiLF) {

  # combine across sexes

  rx <- rf[ff]+1
  fx <- which(ffnm==ff)
  idx <- which(lfcov[,2] == fx-1)
  yx <- sort(lfcov[idx,1])
  sx <- rep$sexrat[yx+1,,,rx]
  px <- rep$ply[yx+1,,,fx]
  phat <- apply(sx*px,c(1,2),sum)
  phat <- t(apply(phat,1,function(x){x <- x/sum(x)}))
  neff <- lfdat[idx,1]/phiLF
  nhat <- apply(phat,2,function(x,neff){x <- x*neff},neff)
  nobs <- lfdat[idx,-1]/phiLF

  # translate to dataframe

  yz <- yx+min(yrs)
  df <- expand.grid(obs=NA,pred=NA,y=range(yz)[1]:range(yz)[2],l=mulbins)
  for(y in 1:length(yz)) {

    df[df$y == yz[y],'obs'] <- nobs[y,]
    df[df$y == yz[y],'pred'] <- nhat[y,]
 
  }

  pp <- ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y)+geom_line(aes(x=l,y=pred),colour='blue')+xlab("length")+ylab("frequency")
  pp

}

# length-conditional age data
# report object, fishery, sex, years, lengths, ages

plot.agelendat <- function(alf,alfcov,rep,ff,sex) {

  # plot mean age (obs and pred) for each length class


  fx <- which(ffnm==ff)
  idx <- which(alfcov[,2] == fx-1)
  alfcovx <- alfcov[idx,]
  alfx <- alf[idx,]
  sx <- ifelse(sex=='female',0,1)
  idx2 <- which(alfcovx[,3] == sx)
  alfcovx <- alfcovx[idx2,]
  alfx <- alfx[idx2,]
  yx <- sort(alfcovx[,1])
  yz <- min(yrs)+unique(yx)
  
  df <- expand.grid(obs=NA,mu=NA,lq=NA,uq=NA,y=range(yz)[1]:range(yz)[2],l=mulbins)

  for(i in 1:dim(alfcovx)[1]) {

    # observed mean age (if there are data in that length bin)     

    nobs <- alfx[i,-1]
    pobs <- nobs/sum(nobs)
    yy <- alfcovx[i,1]+1
    ll <- alfcovx[i,4]+1
    muaobs <- sum(ages*pobs)  
    df[df$y == yrs[yy] & df$l == mulbins[ll],'obs'] <- muaobs 
    phat <- rep$paaly[yy,,ll,sx+1,fx]  
    muahat <- sum(phat*ages)
    sdahat <- sqrt(sum(phat*(ages-muahat)^2))
    df[df$y == yrs[yy] & df$l == mulbins[ll],'mu'] <- muahat
    df[df$y == yrs[yy] & df$l == mulbins[ll],'lq'] <- ifelse(muahat-1.96*sdahat > 0,muahat-1.96*sdahat,0)
    df[df$y == yrs[yy] & df$l == mulbins[ll],'uq'] <- muahat+1.96*sdahat 
    
  }

   # plot it

  pp <- ggplot(df,aes(x=l,y=obs))+geom_point(colour='magenta',size=1)+facet_wrap(~y)+geom_line(aes(x=l,y=mu),colour='blue')+geom_line(aes(x=l,y=uq),linetype='dashed',colour='blue')+geom_line(aes(x=l,y=lq),linetype='dashed',colour='blue')+xlab("length")+ylab("mean age")
  pp 

}

##############################
# plot the tagging data fits #
##############################

# by release event (across length-class and recap area)

plot.tagfits.relyr <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2]
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  totrec <- apply(tagrec,c(1,2),sum)
  prhat <- rep$prhat[,-1,]
  zz <- apply(prhat,c(2,3),function(x,tagrel){x*tagrel},tagrel)
  totrhat <- apply(zz,c(1,2),sum)

  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      Robs <- apply(totrec[idx,],2,sum)
      Rhat <- apply(totrhat[idx,],2,sum)
      xdf[xdf$yrel == y,'obs'] <- Robs[1:tagnrec[idx][1]]
      xdf[xdf$yrel == y,'pred'] <- Rhat[1:tagnrec[idx][1]]

    }

  }

  pp <- ggplot(xdf,aes(x=yrec,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrec,y=pred),colour='magenta',shape=17)+facet_wrap(~yrel)+xlab("year")+ylab("Recaptures")+theme(axis.text.x=element_text(angle=45, hjust=1))
  pp

}

# total by rel year (across rec year, length-class and recap area)

plot.tagfits.relyr.tot <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2]
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  totrec <- apply(tagrec,c(1,2),sum)
  prhat <- rep$prhat[,-1,]
  zz <- apply(prhat,c(2,3),function(x,tagrel){x*tagrel},tagrel)
  totrhat <- apply(zz,c(1,2),sum)

  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      Robs <- apply(totrec[idx,],2,sum)
      Rhat <- apply(totrhat[idx,],2,sum)
      xdf[xdf$yrel == y,'obs'] <- Robs[1:tagnrec[idx][1]]
      xdf[xdf$yrel == y,'pred'] <- Rhat[1:tagnrec[idx][1]]

    }

  }   

  ydf <- aggregate(xdf$obs,by=list(yrel=xdf$yrel),FUN=sum,na.rm=T) 
  zdf <- aggregate(xdf$pred,by=list(yrel=xdf$yrel),FUN=sum,na.rm=T)
  wdf <- data.frame(yrel=ydf$yrel,obs=ydf$x,pred=zdf$x)

  pp <- ggplot(wdf,aes(x=yrel,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrel,y=pred),colour='magenta',shape=17)+xlab("year")+ylab("Recaptures")
  pp

}

# total by rec year (across rel year, length-class and recap area)

plot.tagfits.recyr.tot <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2]
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  totrec <- apply(tagrec,c(1,2),sum)
  prhat <- rep$prhat[,-1,]
  zz <- apply(prhat,c(2,3),function(x,tagrel){x*tagrel},tagrel)
  totrhat <- apply(zz,c(1,2),sum)

  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      Robs <- apply(totrec[idx,],2,sum)
      Rhat <- apply(totrhat[idx,],2,sum)
      xdf[xdf$yrel == y,'obs'] <- Robs[1:tagnrec[idx][1]]
      xdf[xdf$yrel == y,'pred'] <- Rhat[1:tagnrec[idx][1]]

    }

  }

  ydf <- aggregate(xdf$obs,by=list(yrec=xdf$yrec),FUN=sum,na.rm=T) 
  zdf <- aggregate(xdf$pred,by=list(yrec=xdf$yrec),FUN=sum,na.rm=T)
  wdf <- data.frame(yrec=ydf$yrec,obs=ydf$x,pred=zdf$x)

  pp <- ggplot(wdf,aes(x=yrec,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrec,y=pred),colour='magenta',shape=17)+xlab("year")+ylab("Recaptures")
  pp

}

plot.tagfits.recyrarea <- function(rep,tagrel,tagnrec,tagcov,tagrec,nrecmax) {

  nt <- length(tagrel)
  yrng <- range(unique(tagcov[,3]+yrs[1]))
  yrel <- yrng[1]:yrng[2] 
  xdf <- expand.grid(yrel=yrel,yrec=yrel+1,rrec=c('North','South'),obs=NA,pred=NA)
  xdf <- subset(xdf,yrec>yrel)
  xdf <- subset(xdf,yrec-yrel<=nrecmax)
  robs <- tagrec
  rhat <- apply(rep$prhat[,-1,],c(2,3),function(x,tagrel){x <- x*tagrel},tagrel)

  rnm <- c("North","South")
  for(y in yrel) {

    iy <- y-yrs[1]
    idx <- which(tagcov[,3]==iy)
    nT <- sum(tagrel[idx])
    if(nT > 0) {
   
      for(r in 1:nr) {

        Robs <- apply(robs[idx,,r],2,sum)
        Rhat <- apply(rhat[idx,,r],2,sum)
        xdf[xdf$yrel == y & xdf$rrec == rnm[r],'obs'] <- Robs[1:tagnrec[idx][1]]
        xdf[xdf$yrel == y & xdf$rrec == rnm[r],'pred'] <- Rhat[1:tagnrec[idx][1]]

      }
    }

  } 

  ydf <- aggregate(xdf$obs,by=list(yrec=xdf$yrec,rrec=xdf$rrec),FUN=sum,na.rm=T) 
  zdf <- aggregate(xdf$pred,by=list(yrec=xdf$yrec,rrec=xdf$rrec),FUN=sum,na.rm=T)
  wdf <- data.frame(yrec=ydf$yrec,rrec=ydf$rrec,obs=ydf$x,pred=zdf$x)

  pp <- ggplot(wdf,aes(x=yrec,y=obs))+geom_point(colour='blue')+geom_point(aes(x=yrec,y=pred),colour='magenta',shape=17)+facet_wrap(~rrec,scales='free_y')+xlab("year")+ylab("Recaptures")+theme(axis.text.x=element_text(angle=45, hjust=1))
  pp

}

###############################
# over-dispersion diagnostics #
###############################

# LF data

get.phiLFhat <- function(lfdat,lfcov,rep,ff,philf,weighted=T) {

  # combine across sexes

  rx <- rf[ff]+1
  fx <- which(ffnm==ff)
  idx <- which(lfcov[,2] == fx-1)
  yx <- sort(lfcov[idx,1])
  sx <- rep$sexrat[yx+1,,,rx]
  px <- rep$ply[yx+1,,,fx]
  phat <- apply(sx*px,c(1,2),sum)
  phat <- t(apply(phat,1,function(x){x <- x/sum(x)}))
  neff <- lfdat[idx,1]
  nhat <- apply(phat,2,function(x,neff){x <- x*neff},neff)
  nobs <- lfdat[idx,-1]

  phihat <- rep(NA,length(yx))

  for(y in 1:length(yx)) {

    xx <- nhat[y,]
    yy <- nobs[y,]
    wx <- phat[y,]
    epsx <- (yy-xx)/sqrt((philf*neff[y]*wx*(1-wx)))
    if(weighted) phihat[y] <- sum(wx*(epsx-(sum(epsx*wx)))^2)
    else phihat[y] <- var(epsx)
    
  }

  return(phihat)

}

# tagging data

get.phiT <- function(rep,tagrel,tagrec,weighted=T) {

    nt <- length(tagrel)
    wx <- epsx <- rep(NA,nt)

    for(t in 1:nt) {

      tobs <- tagrel[t]
      rsum <- sum(tagrec[t,,])
      psum <- sum(rep$prhat[t,-1,])
      rhat <- tobs*sum(rep$prhat[t,-1,])
      wx[t] <- rsum
      epsx[t] <- (rsum-rhat)/sqrt(tobs*psum*(1-psum))

    }

    wx <- wx/sum(wx)
    if(weighted) phihat <- sum(wx*(epsx-(sum(epsx*wx)))^2)
    else phihat <- var(epsx)

    return(phihat)

}

# reconstruction code

recon <- function(nits=1,parm,parv,parnm1,parnm2,mat.type='new') {
  
  a.lw <- 4.4e-6
  b.lw <- 3.14
  if(mat.type == 'old') {
    
    ml50 <- 139.6
    ml95 <- 185.8
  
  }
  if(mat.type == 'new') {
    
    nx <- data$matpars[1]
    mux <- data$matpars[2]
    
  }
  
  # arrays
  
  N <- array(dim=c(ny,na,ns,nr))
  Nfin <- hyafin <- array(dim=c(nits,na,ns,nr))
  Rec <- array(dim=c(nits,ny))
  SSB <- array(dim=c(nits,ny,ns,nr))
  SSBf <- matrix(nrow=nits,ncol=ny) 
  hyf <- array(dim=c(ny,nf))
  hyfa <- array(dim=c(ny,na,ns,nf))
  hya <- array(dim=c(ny,na,ns,nr))
  sellen <- array(dim=c(nits,nl,nf))
  sel <- array(dim=c(na,ns,nf))
  m <- w <- array(dim=c(na,ns))
  mula <- array(dim=c(na,ns))
  pla <- array(dim=c(ns,nl,na))
  
  # definitely non-MCMC stuff
  
  hh <- data$hh
  M <- data$M
  zeta <- data$zeta
  hmax <- data$hmax
  nages <- length(ages)
  rff <- rf+1
  recymin <- recyr[1]+1
  recymax <- recyr[2]+1
  
  for(i in 1:nits) {

    # parameters
    
    R0 <- ifelse(any(parnm1 == 'logR0'),exp(parm[i,grep('logR0',colnames(parm))]),exp(parv[['logR0']]))
    if(any(parnm1 == 'epsR')) {
      epsR <- unname(parm[i,grep('epsR',colnames(parm))])
    } else {
      epsR <- unname(parv[['epsR']])
    }
    eta1 <- ifelse(any(parnm1 == 'ilogeta'),ilogit(parm[i,grep('ilogeta',colnames(parm))]),ilogit(parv[['ilogeta']])) 
    eta <- c(eta1,1-eta1)  
    if(any(parnm1 == 'ilogpim')) {
      pimv <- unname(ilogit(parm[i,grep('ilogpim',colnames(parm))]))
    } else {
      pimv <- unname(ilogit(parv[['ilogpim']]))
    }
    pim <- matrix(c(pimv[1],1-pimv[1],1-pimv[2],pimv[2]),nrow=2,ncol=2,byrow=T)
    if(any(parnm1 == 'selparsATT')) {
      selparsATT <- unname(exp(parm[i,grep('selparsATT',colnames(parm))]))
    } else {
      selparsATT <- unname(exp(parv[['selparsATT']]))
    }
    if(any(parnm1 == 'selparsNVT')) {
      selparsNVT <- unname(exp(parm[i,grep('selparsNVT',colnames(parm))]))
    } else {
      selparsNVT <- unname(exp(parv[['selparsNVT']]))
    }
    if(any(parnm1 == 'selparsATL')) {
      selparsATL <- unname(exp(parm[i,grep('selparsATL',colnames(parm))]))
    } else {
      selparsATL <- unname(exp(parv[['selparsATL']]))
    }
    if(any(parnm1 == 'selparsNMRL')) {
      selparsNMRL <- unname(exp(parm[i,grep('selparsNMRL',colnames(parm))]))
    } else {
      selparsNMRL <- unname(exp(parv[['selparsNMRL']]))
    }
    if(any(parnm1 == 'selparsSMRL')) {
      selparsSMRL <- unname(exp(parm[i,grep('selparsSMRL',colnames(parm))]))
    } else {
      selparsSMRL <- unname(exp(parv[['selparsSMRL']]))
    }
    if(any(parnm1 == 'logl1')) {
      l1 <- unname(exp(parm[i,grep('logl1',colnames(parm))]))
    } else {
      l1 <- unname(exp(parv[['logl1']]))
    }
    if(any(parnm1 == 'logl2')) {
      l2 <- unname(exp(parm[i,grep('logl2',colnames(parm))]))
    } else {
      l2 <- unname(exp(parv[['logl2']]))
    }
    if(any(parnm1 == 'logk')) {
      k <- unname(exp(parm[i,grep('logk',colnames(parm))]))
    } else {
      k <- unname(exp(parv[['logk']]))
    }
    if(any(parnm1 == 'logsdla')) {
      sdla <- unname(exp(parm[i,grep('logsdla',colnames(parm))]))
    } else {
      sdla <- unname(exp(parv[['logsdla']]))
    }
    
    # mean length-at-age
  
    mula[,1] <- l1[1]+(l2[1]-l1[1])*(1-exp(-k[1]*(ages-aref[1])))/(1-exp(-k[1]*(aref[2]-aref[1])))
    mula[,2] <- l1[2]+(l2[2]-l1[2])*(1-exp(-k[2]*(ages-aref[1])))/(1-exp(-k[2]*(aref[2]-aref[1])))
    
    # P(l | a)
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        psum <-sum(dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE))   
        pla[s,,a] <- dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE)/psum
        
       } 
    }
    
    # estimate weight, maturity and selectivity at age from length
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        # weights dl
        
        ll <- mula[a,s]*exp(-1.96*sdla[s])
        lu <- mula[a,s]*exp(1.96*sdla[s])
        ldel <- (lu-ll)/25
        lref <- ll+(0:24)*ldel
        dl <- dnorm(log(lref),log(mula[a,s]),sdla[s],FALSE)
        dl <- dl/sum(dl)
        
        # maturity and weight-at-age
        
        if(mat.type == 'old') ml <- 1/(1+19^{-(lref-ml50)/(ml95-ml50)})
        if(mat.type == 'new') ml <- lref^nx/(mux^nx+lref^nx)
        wl <- a.lw*lref^b.lw
        m[a,s] <- sum(dl*ml)
        w[a,s] <- sum(dl*wl)/1e+3
        
        ###############
        # selectivity #
        ###############
        
        # 1. ATT (double logistic with dr95 fixed at 3)
        
        sl50 <- selparsATT[1]
        dl95 <- selparsATT[2]
        dr95 <- 3
        sr50 <- sl50+dl95+dr95+selparsATT[3]
        sellen[i,,1] <- (1+19^{-(mulbins-sl50)/dl95})^{-1}*(1-(1+19^{-(mulbins-sr50)/dr95})^{-1})
        seltmp <- (1+19^{-(lref-sl50)/dl95})^{-1}*(1-(1+19^{-(lref-sr50)/dr95})^{-1})
        sel[a,s,1] <- sum(seltmp*dl)
        
        # 2. NVT (gamma with z fixed at 1)
        
        ga <- selparsNVT[1]
        gb <- selparsNVT[2]
        gmax <- (ga/gb)^ga * exp(-ga)
        sellen[i,,2] <- mulbins^ga * exp(-gb*mulbins)/gmax
        seltmp <- lref^ga * exp(-gb*lref)/gmax
        sel[a,s,2] <- sum(seltmp*dl)
          
        # 3. ATL (logistic)
        
        s50 <- selparsATL[1]
        s95 <- s50+selparsATL[2]
        sellen[i,,3] <- 1/(1+19^{-(mulbins-s50)/(s95-s50)})
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,3] <- sum(seltmp*dl)
        
        # 4. NMRL (logistic)
        
        s50 <- selparsNMRL[1]
        s95 <- s50+selparsNMRL[2]
        sellen[i,,4] <- 1/(1+19^{-(mulbins-s50)/(s95-s50)})
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,4] <- sum(seltmp*dl)
        
        # 5. SMRL (logistic)
        
        s50 <- selparsSMRL[1]
        s95 <- s50+selparsSMRL[2]
        sellen[i,,5] <- 1/(1+19^{-(mulbins-s50)/(s95-s50)})
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,5] <- sum(seltmp*dl)
      }
    }
    
    # initial population dynamics
    
    Rec[i,1] <- R0
    for(s in 1:ns) 
      for(r in 1:nr) N[1,1,s,r] <- R0 * eta[r] * zeta[s]
    for(a in 2:(nages-1))
      for(s in 1:ns)
        for(r in 1:nr) N[1,a,s,r] <- sum(pim[,r]*N[1,a-1,s,]*exp(-M))
    
    av <- rep(NA,nr)
    B <- matrix(nrow=nr,ncol=nr)
    for(s in 1:ns) {
      
      for(r in 1:nr) av[r] <- sum(pim[,r]*N[1,nages-1,s,]*exp(-M))
      for(r in 1:nr) {
        for(ss in 1:nr) { 
          B[ss,r] <- pim[ss,r]*exp(-M)
        }
      }
      
      npg <- t(av) %*% solve(diag(1,nrow=nr)-B)
      N[1,nages,s,] <- npg
        
    }
    
    B0 <- sum(apply(N[1,,1,],1,sum)*w[,1]*m[,1])
    SSBf[i,1] <- B0
    for(s in 1:ns)
      for(r in 1:nr) SSB[i,1,s,r] <- sum(N[1,,s,r]*w[,s]*m[,s])
    rho <- B0/R0
    alp <- 4*hh/(rho*(1-hh))
    bet <- (5*hh-1)/(B0*(1-hh))
    
    # initial harvest rates
    
    for(f in 1:nf) {
      
      Xtmp <- 0
      for(s in 1:ns) Xtmp <- Xtmp+sum(N[1,,s,rff[[f]]]*exp(-tau[1,f]*M)*sel[,s,f]*w[,s])
      hyf[1,f] <- C[1,f]/Xtmp
      
      for(s in 1:ns) hyfa[1,,s,f] <- hyf[1,f] * sel[,s,f]
    }
    
    # if there's catch in first year do hmax checks
      
    if(any(C[1,] > 0.01)) {
        
        for(s in 1:ns) {
          for(a in 1:nages) {
            for(r in 1:nr) {
        
              fsub <- names(rff[rff == r])
              csum <- sum(C[1,fsub])
              hsum <- 0
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                hsum <- hsum+sum(hyfa[1,a,s,fx])
              }
              
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                if(hsum > hmax) hyfa[1,a,s,fx] <- hmax*C[1,fx]/csum
            }
          }
        }
      }
    }
    
    # total harvest rates by region
    
    for(r in 1:nr) {
     
      nfr <- length(rff[rff == r]) 
      fv <- rep(NA,nfr)
      for(f in 1:nfr) {
        
        fsub <- names(rff[rff == r])
        fv[f] <- grep(fsub[f],ffnm)
        
      }
      
      hsub <- hyfa[1,,,fv]
      hya[1,,,r] <- apply(hsub,c(1,2),sum)
    }
    
    # annual loop
    
    for(y in 2:ny) {
      
      # recruitment (different for years where recruitment estimated)
    
      Rtot <- alp*SSBf[i,y-1]/(1+bet*SSBf[i,y-1])
      for(s in 1:ns) {
        for(r in 1:nr) {
          
          if(y < recymin | y > recymax) {
            
            N[y,1,s,r] <- Rtot * eta[r] * zeta[s]
            
          } else {
            
            N[y,1,s,r] <- Rtot * eta[r] * zeta[s] *exp(epsR[y-recymin+1]-sigmaR^2/2)
          }
        }
      }
       
      Rec[i,y] <- sum(N[y,1,,]) 
      # post-recruit dynamics
      
      # a = 2 to A+-1
      
      for(s in 1:ns) 
        for(a in 2:(nages-1)) 
          for(r in 1:nr) N[y,a,s,r] <- sum(pim[,r]*N[y-1,a-1,s,]*exp(-M)*(1-hya[y-1,a-1,s,]))
    
      # plus group
          
      for(s in 1:ns)
        for(r in 1:nr) N[y,nages,s,r] <- sum(pim[,r]*N[y-1,nages-1,s,]*exp(-M)*(1-hya[y-1,nages-1,s,]))+sum(pim[,r]*N[y-1,nages,s,]*exp(-M)*(1-hya[y-1,nages,s,]))
                          
      # SSB
      
      SSBf[i,y] <- sum(apply(N[y,,1,],1,sum)*w[,1]*m[,1])
      for(s in 1:ns)
        for(r in 1:nr) SSB[i,y,s,r] <- sum(N[y,,s,r]*w[,s]*m[,s])
      
      # harvest rates
     
      for(f in 1:nf) {
        
        Xtmp <- 0
        for(s in 1:ns) Xtmp <- Xtmp+sum(N[y,,s,rff[[f]]]*exp(-tau[y,f]*M)*sel[,s,f]*w[,s])
        hyf[y,f] <- C[y,f]/Xtmp
        
        for(s in 1:ns) hyfa[y,,s,f] <- hyf[y,f] * sel[,s,f]
      }
      
      # if there's catch in first year do hmax checks
      
      if(any(C[y,] > 0.01)) {
        
        for(s in 1:ns) {
          for(a in 1:nages) {
            for(r in 1:nr) {
              
              fsub <- names(rff[rff == r])
              csum <- sum(C[y,fsub])
              hsum <- 0
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                hsum <- hsum+sum(hyfa[y,a,s,fx])
              }
              
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                if(hsum > hmax) hyfa[y,a,s,fx] <- hmax*C[y,fx]/csum
              }
            }
          }
        }
      }
      
      # total harvest rates by region
      
      for(r in 1:nr) {
        
        nfr <- length(rff[rff == r]) 
        fv <- rep(NA,nfr)
        for(f in 1:nfr) {
          
          fsub <- names(rff[rff == r])
          fv[f] <- grep(fsub[f],ffnm)
          
        }
        
        hsub <- hyfa[y,,,fv]
        hya[y,,,r] <- apply(hsub,c(1,2),sum)
      } 
    }
    
    Nfin[i,,,] <- N[ny,,,]
    hyafin[i,,,] <- hya[ny,,,]

    if(i %% 50 == 0) cat("iteration",i,"of",nits,"\n")
    
  }
  
  return(list(SSB=SSB,SSBf=SSBf,Rec=Rec,Nfin=Nfin,hyafin=hyafin,sellen=sellen))
}

# projection code

proj <- function(nits = 1,nyp = 1,parm,parv,parnm1,parnm2,B0,nfin,hyafin,taufin,tacsplit,recdev='parametric',mat.type='new') {
  
  a.lw <- 4.4e-6
  b.lw <- 3.14
  if(mat.type == 'old') {
    
    ml50 <- 139.6
    ml95 <- 185.8
    
  }
  if(mat.type == 'new') {
    
    nx <- data$matpars[1]
    mux <- data$matpars[2]
    
  }
  
  # arrays
  
  N <- array(dim=c(nyp,na,ns,nr))
  Rec <- array(dim=c(nits,nyp))
  SSBf <- Rec <- matrix(nrow=nits,ncol=nyp) 
  hyf <- array(dim=c(nyp,nf))
  hyfa <- array(dim=c(nyp,na,ns,nf))
  hya <- array(dim=c(nyp,na,ns,nr))
  sel <- array(dim=c(na,ns,nf))
  m <- w <- array(dim=c(na,ns))
  mula <- array(dim=c(na,ns))
  pla <- array(dim=c(ns,nl,na))
  
  # definitely non-MCMC stuff
  
  hh <- data$hh
  M <- data$M
  zeta <- data$zeta
  sigmaR <- data$sigmaR
  hmax <- data$hmax
  nages <- length(ages)
  rff <- rf+1
  Cf <- tacsplit
  
  for(i in 1:nits) {
    
    # parameters
    
    R0 <- ifelse(any(parnm1 == 'logR0'),exp(parm[i,grep('logR0',colnames(parm))]),exp(parv[['logR0']]))
    if(any(parnm1 == 'epsR')) {
      epsR <- unname(parm[i,grep('epsR',colnames(parm))])
    } else {
      epsR <- unname(parv[['epsR']])
    }
    eta1 <- ifelse(any(parnm1 == 'ilogeta'),ilogit(parm[i,grep('ilogeta',colnames(parm))]),ilogit(parv[['ilogeta']])) 
    eta <- c(eta1,1-eta1)  
    if(any(parnm1 == 'ilogpim')) {
      pimv <- unname(ilogit(parm[i,grep('ilogpim',colnames(parm))]))
    } else {
      pimv <- unname(ilogit(parv[['ilogpim']]))
    }
    pim <- matrix(c(pimv[1],1-pimv[1],1-pimv[2],pimv[2]),nrow=2,ncol=2,byrow=T)
    if(any(parnm1 == 'selparsATT')) {
      selparsATT <- unname(exp(parm[i,grep('selparsATT',colnames(parm))]))
    } else {
      selparsATT <- unname(exp(parv[['selparsATT']]))
    }
    if(any(parnm1 == 'selparsNVT')) {
      selparsNVT <- unname(exp(parm[i,grep('selparsNVT',colnames(parm))]))
    } else {
      selparsNVT <- unname(exp(parv[['selparsNVT']]))
    }
    if(any(parnm1 == 'selparsATL')) {
      selparsATL <- unname(exp(parm[i,grep('selparsATL',colnames(parm))]))
    } else {
      selparsATL <- unname(exp(parv[['selparsATL']]))
    }
    if(any(parnm1 == 'selparsNMRL')) {
      selparsNMRL <- unname(exp(parm[i,grep('selparsNMRL',colnames(parm))]))
    } else {
      selparsNMRL <- unname(exp(parv[['selparsNMRL']]))
    }
    if(any(parnm1 == 'selparsSMRL')) {
      selparsSMRL <- unname(exp(parm[i,grep('selparsSMRL',colnames(parm))]))
    } else {
      selparsSMRL <- unname(exp(parv[['selparsSMRL']]))
    }
    if(any(parnm1 == 'logl1')) {
      l1 <- unname(exp(parm[i,grep('logl1',colnames(parm))]))
    } else {
      l1 <- unname(exp(parv[['logl1']]))
    }
    if(any(parnm1 == 'logl2')) {
      l2 <- unname(exp(parm[i,grep('logl2',colnames(parm))]))
    } else {
      l2 <- unname(exp(parv[['logl2']]))
    }
    if(any(parnm1 == 'logk')) {
      k <- unname(exp(parm[i,grep('logk',colnames(parm))]))
    } else {
      k <- unname(exp(parv[['logk']]))
    }
    if(any(parnm1 == 'logsdla')) {
      sdla <- unname(exp(parm[i,grep('logsdla',colnames(parm))]))
    } else {
      sdla <- unname(exp(parv[['logsdla']]))
    }
    
    # mean length-at-age
    
    mula[,1] <- l1[1]+(l2[1]-l1[1])*(1-exp(-k[1]*(ages-aref[1])))/(1-exp(-k[1]*(aref[2]-aref[1])))
    mula[,2] <- l1[2]+(l2[2]-l1[2])*(1-exp(-k[2]*(ages-aref[1])))/(1-exp(-k[2]*(aref[2]-aref[1])))
    
    # P(l | a)
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        psum <-sum(dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE))   
        pla[s,,a] <- dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE)/psum
        
      } 
    }
    
    # estimate weight, maturity and selectivity at age from length
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        # weights dl
        
        ll <- mula[a,s]*exp(-1.96*sdla[s])
        lu <- mula[a,s]*exp(1.96*sdla[s])
        ldel <- (lu-ll)/25
        lref <- ll+(0:24)*ldel
        dl <- dnorm(log(lref),log(mula[a,s]),sdla[s],FALSE)
        dl <- dl/sum(dl)
        
        # maturity and weight-at-age
    
        if(mat.type == 'old') ml <- 1/(1+19^{-(lref-ml50)/(ml95-ml50)})
        if(mat.type == 'new') ml <- lref^nx/(mux^nx+lref^nx)
        wl <- a.lw*lref^b.lw
        m[a,s] <- sum(dl*ml)
        w[a,s] <- sum(dl*wl)/1e+3
        
        ###############
        # selectivity #
        ###############
        
        # 1. ATT (double logistic with dr95 fixed at 3)
        
        sl50 <- selparsATT[1]
        dl95 <- selparsATT[2]
        dr95 <- 3
        sr50 <- sl50+dl95+dr95+selparsATT[3]
        seltmp <- (1+19^{-(lref-sl50)/dl95})^{-1}*(1-(1+19^{-(lref-sr50)/dr95})^{-1})
        sel[a,s,1] <- sum(seltmp*dl)
        
        # 2. NVT (gamma with z fixed at 1)
        
        ga <- selparsNVT[1]
        gb <- selparsNVT[2]
        gmax <- (ga/gb)^ga * exp(-ga)
        seltmp <- lref^ga * exp(-gb*lref)/gmax
        sel[a,s,2] <- sum(seltmp*dl)
        
        # 3. ATL (logistic)
        
        s50 <- selparsATL[1]
        s95 <- s50+selparsATL[2]
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,3] <- sum(seltmp*dl)
        
        # 4. NMRL (logistic)
        
        s50 <- selparsNMRL[1]
        s95 <- s50+selparsNMRL[2]
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,4] <- sum(seltmp*dl)
        
        # 5. SMRL (logistic)
        
        s50 <- selparsSMRL[1]
        s95 <- s50+selparsSMRL[2]
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,5] <- sum(seltmp*dl)
      }
    }
    
    # initial numbers-at-age
    
    N[1,,,] <- Nfin[i,,,]
    hya[1,,,] <- hyafin[i,,,]
    SSBf[i,1] <- sum(apply(N[1,,1,],1,sum)*w[,1]*m[,1])
    Rec[i,1] <- sum(N[1,1,,])
    
    # project forwards
    
    for(y in 2:nyp) {
      
      # recruitment (different for years where recruitment estimated)
      
      rho <- B0[i]/R0
      alp <- 4*hh/(rho*(1-hh))
      bet <- (5*hh-1)/(B0[i]*(1-hh))
      Rtot <- alp*SSBf[i,y-1]/(1+bet*SSBf[i,y-1])
      for(s in 1:ns) {
        for(r in 1:nr) {
          
          if(recdev == 'parametric') N[y,1,s,r] <- Rtot * eta[r] * zeta[s] * rlnorm(1,-sigmaR^2/2,sigmaR)
          if(recdev == 'nonparametric') {
             
            epsr <- sample(epsR,1,replace=T)
            N[y,1,s,r] <- Rtot * eta[r] * zeta[s] * exp(epsr-sigmaR^2/2)
            
          }
        }
      }
      Rec[i,y] <- sum(N[y,1,,]) 
      
      # post-recruit dynamics
      
      # a = 2 to A+-1
      
      for(s in 1:ns) 
        for(a in 2:(nages-1)) 
          for(r in 1:nr) N[y,a,s,r] <- sum(pim[,r]*N[y-1,a-1,s,]*exp(-M)*(1-hya[y-1,a-1,s,]))
      
      # plus group
      
      for(s in 1:ns)
        for(r in 1:nr) N[y,nages,s,r] <- sum(pim[,r]*N[y-1,nages-1,s,]*exp(-M)*(1-hya[y-1,nages-1,s,]))+sum(pim[,r]*N[y-1,nages,s,]*exp(-M)*(1-hya[y-1,nages,s,]))
      
      # SSB
      
      SSBf[i,y] <- sum(apply(N[y,,1,],1,sum)*w[,1]*m[,1])
      #for(s in 1:ns)
      #  for(r in 1:nr) SSB[i,y,s,r] <- sum(N[y,,s,r]*w[,s]*m[,s])
      
      # harvest rates
      
      for(f in 1:nf) {
        
        Xtmp <- 0
        for(s in 1:ns) Xtmp <- Xtmp+sum(N[y,,s,rff[[f]]]*exp(-taufin[f]*M)*sel[,s,f]*w[,s])
        hyf[y,f] <- Cf[f]/Xtmp
        
        for(s in 1:ns) hyfa[y,,s,f] <- hyf[y,f] * sel[,s,f]
      }
      
      # if there's catch in first year do hmax checks
      
      if(any(Cf > 0.01)) {
        
        for(s in 1:ns) {
          for(a in 1:nages) {
            for(r in 1:nr) {
              
              fsub <- names(rff[rff == r])
              csum <- sum(Cf[fsub])
              hsum <- 0
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                hsum <- hsum+sum(hyfa[y,a,s,fx])
              }
              
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                if(hsum > hmax) hyfa[y,a,s,fx] <- hmax*Cf[fx]/csum
              }
            }
          }
        }
      }
      
      # total harvest rates by region
      
      for(r in 1:nr) {
        
        nfr <- length(rff[rff == r]) 
        fv <- rep(NA,nfr)
        for(f in 1:nfr) {
          
          fsub <- names(rff[rff == r])
          fv[f] <- grep(fsub[f],ffnm)
          
        }
        
        hsub <- hyfa[y,,,fv]
        hya[y,,,r] <- apply(hsub,c(1,2),sum)
      } 
    }
   
    if(i %% 50 == 0) cat("iteration",i,"of",nits,"\n")
  }
  
  return(list(SSBf=SSBf,Rec=Rec))
}

# get tac split given total TAC, ATT share and N-S share

get.tacsplit <- function(tac,tacatl,nssplit,ffnm) {
  
  tsplit <- rep(NA,length(ffnm))
  names(tsplit) <- ffnm
  tsplit['ATT'] <- tsplit['NVT'] <- 0
  tsplit['ATL'] <- tacatl
  trem <- tac-tacatl
  tsplit['NMRL'] <- trem * nssplit[1]
  tsplit['SMRL'] <- trem * nssplit[2]
  return(tsplit)
  
}

# create transition matrices

get.tmat <- function(gdf,lbins) {

  nbins <- length(lbins)-1
  T <- array(dim=c(2,nbins,nbins)) 
  tau <- 1 
  Linf <- gdf$Linf
  k <- gdf$k
  ns <- 2

  for(s in 1:ns) {
    for(i in 1:nbins) {

      lx <- lbins[i]
      ly <- lbins[i+1]

      for(j in 1:nbins) {

        llj <- lbins[j]
        luj <- lbins[j+1]
          
        lli <- lx + lvbinc(lx,tau,k[s],Linf[s]) 
        lui <- ly + lvbinc(ly,tau,k[s],Linf[s]) 
      
        # need to work out Lebesgue measure of intersection 
        # of image and actual length bin / length bin

        if(lli > llj & lui < luj) {
          ptmp <- 1
        } else {
          tmp <- c(max(llj,lli),min(luj,lui))
          mu <- ifelse(tmp[1] < tmp[2],tmp[2]-tmp[1],0)
          nu <- lui-lli
          ptmp <- mu/nu 
        }
        
        T[s,i,j] <- ptmp
      }
    }
  } 
 
  return(T)

}

# growth interval

lvbinc <- function(lrel,tau,k,Linf) {

  return(max((Linf-lrel)*(1-exp(-k*tau)),1e-9*lrel))
}

# hybrid profile likelihood/Laplace approx. estimation of sigmaR

est.sigmaR <- function(ddata,parss,dll,vsig,sub=TRUE) {

  nsig <- length(vsig)
  mnlogl <- rep(NA,nsig)

  mapp <- list(lsigmaR = factor(NA))

  for(i in 1:nsig) {

    # create AD object for given sigmaR value

    parss$lsigmaR <- log(vsig[i])
    objj <- MakeADFun(data=ddata,parameters=parss,map=mapp,DLL=dll)

    # estimate parameters

    ress <- do.call("optim",objj)

    # get Hessian at MLE

    hess <- objj$he()

    if(sub) {
      
      # subset for epsR only

      idx <- grep("epsR",names(ress$par))
      xhess <- hess[idx,idx]

    } else {

      xhess <- hess

    }

    # approximate full marginal likelihood

    nrec <- length(recyr[1]:recyr[2])
    mnlogl[i] <- -nrec*log(sqrt(2*pi))+0.5*log(det(xhess))+ress$value

    cat(i,"of",nsig,"\n")

  }

  return(mnlogl)

}

# v2 MCMC reconstruction code

recon2 <- function(nits=1,parm,parv,parnm1,parnm2,mat.type='new') {
  
  a.lw <- 4.4e-6
  b.lw <- 3.14
  if(mat.type == 'old') {
    
    ml50 <- 139.6
    ml95 <- 185.8
  
  }
  if(mat.type == 'new') {
    
    nx <- data$matpars[1]
    mux <- data$matpars[2]
    
  }
  
  # arrays
  
  N <- array(dim=c(ny,na,ns,nr))
  Nfin <- hyafin <- array(dim=c(nits,na,ns,nr))
  Rec <- array(dim=c(nits,ny))
  SSB <- array(dim=c(nits,ny,ns,nr))
  SSBf <- matrix(nrow=nits,ncol=ny) 
  hyf <- array(dim=c(ny,nf))
  hyfa <- array(dim=c(ny,na,ns,nf))
  hya <- array(dim=c(ny,na,ns,nr))
  sellen <- array(dim=c(nits,nl,nf))
  sel <- array(dim=c(na,ns,nf))
  m <- w <- array(dim=c(na,ns))
  mula <- array(dim=c(na,ns))
  pla <- array(dim=c(ns,nl,na))
  
  # definitely non-MCMC stuff
  
  hh <- data$hh
  M <- data$M
  zeta <- data$zeta
  hmax <- data$hmax
  nages <- length(ages)
  rff <- rf+1
  recymin <- recyr[1]+1
  recymax <- recyr[2]+1
  
  for(i in 1:nits) {

    # parameters
    
    R0 <- ifelse(any(parnm1 == 'logR0'),exp(parm[i,grep('logR0',colnames(parm))]),exp(parv[['logR0']]))
    if(any(parnm1 == 'epsR')) {
      epsR <- unname(parm[i,grep('epsR',colnames(parm))])
    } else {
      epsR <- unname(parv[['epsR']])
    }
    eta1 <- ifelse(any(parnm1 == 'ilogeta'),ilogit(parm[i,grep('ilogeta',colnames(parm))]),ilogit(parv[['ilogeta']])) 
    eta <- c(eta1,1-eta1)  
    if(any(parnm1 == 'ilogpim')) {
      pimv <- unname(ilogit(parm[i,grep('ilogpim',colnames(parm))]))
    } else {
      pimv <- unname(ilogit(parv[['ilogpim']]))
    }
    pim <- matrix(c(pimv[1],1-pimv[1],1-pimv[2],pimv[2]),nrow=2,ncol=2,byrow=T)
    if(any(parnm1 == 'selparsATT')) {
      selparsATT <- unname(exp(parm[i,grep('selparsATT',colnames(parm))]))
    } else {
      selparsATT <- unname(exp(parv[['selparsATT']]))
    }
    if(any(parnm1 == 'selparsNVT')) {
      selparsNVT <- unname(exp(parm[i,grep('selparsNVT',colnames(parm))]))
    } else {
      selparsNVT <- unname(exp(parv[['selparsNVT']]))
    }
    if(any(parnm1 == 'selparsATL')) {
      selparsATL <- unname(exp(parm[i,grep('selparsATL',colnames(parm))]))
    } else {
      selparsATL <- unname(exp(parv[['selparsATL']]))
    }
    if(any(parnm1 == 'selparsNMRL')) {
      selparsNMRL <- unname(exp(parm[i,grep('selparsNMRL',colnames(parm))]))
    } else {
      selparsNMRL <- unname(exp(parv[['selparsNMRL']]))
    }
    if(any(parnm1 == 'selparsSMRL')) {
      selparsSMRL <- unname(exp(parm[i,grep('selparsSMRL',colnames(parm))]))
    } else {
      selparsSMRL <- unname(exp(parv[['selparsSMRL']]))
    }
    sigmaR <- ifelse(any(parnm1 == 'lsigmaR'),exp(parm[i,grep('lsigmaR',colnames(parm))]),exp(parv[['lsigmaR']]))

    # growth handled separately now

    l1 <- data$l1
    l2 <- data$l2
    k <- data$k
    sdla <- data$sdla

    # mean length-at-age
  
    mula[,1] <- l1[1]+(l2[1]-l1[1])*(1-exp(-k[1]*(ages-aref[1])))/(1-exp(-k[1]*(aref[2]-aref[1])))
    mula[,2] <- l1[2]+(l2[2]-l1[2])*(1-exp(-k[2]*(ages-aref[1])))/(1-exp(-k[2]*(aref[2]-aref[1])))
    
    # P(l | a)
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        psum <-sum(dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE))   
        pla[s,,a] <- dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE)/psum
        
       } 
    }
    
    # estimate weight, maturity and selectivity at age from length
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        # weights dl
        
        ll <- mula[a,s]*exp(-1.96*sdla[s])
        lu <- mula[a,s]*exp(1.96*sdla[s])
        ldel <- (lu-ll)/25
        lref <- ll+(0:24)*ldel
        dl <- dnorm(log(lref),log(mula[a,s]),sdla[s],FALSE)
        dl <- dl/sum(dl)
        
        # maturity and weight-at-age
        
        if(mat.type == 'old') ml <- 1/(1+19^{-(lref-ml50)/(ml95-ml50)})
        if(mat.type == 'new') ml <- lref^nx/(mux^nx+lref^nx)
        wl <- a.lw*lref^b.lw
        m[a,s] <- sum(dl*ml)
        w[a,s] <- sum(dl*wl)/1e+3
        
        ###############
        # selectivity #
        ###############
        
        # 1. ATT (double logistic with dr95 fixed at 3)
        
        sl50 <- selparsATT[1]
        dl95 <- selparsATT[2]
        dr95 <- 3
        sr50 <- sl50+dl95+dr95+selparsATT[3]
        sellen[i,,1] <- (1+19^{-(mulbins-sl50)/dl95})^{-1}*(1-(1+19^{-(mulbins-sr50)/dr95})^{-1})
        seltmp <- (1+19^{-(lref-sl50)/dl95})^{-1}*(1-(1+19^{-(lref-sr50)/dr95})^{-1})
        sel[a,s,1] <- sum(seltmp*dl)
        
        # 2. NVT (gamma with z fixed at 1)
        
        ga <- selparsNVT[1]
        gb <- selparsNVT[2]
        gmax <- (ga/gb)^ga * exp(-ga)
        sellen[i,,2] <- mulbins^ga * exp(-gb*mulbins)/gmax
        seltmp <- lref^ga * exp(-gb*lref)/gmax
        sel[a,s,2] <- sum(seltmp*dl)
          
        # 3. ATL (logistic)
        
        s50 <- selparsATL[1]
        s95 <- s50+selparsATL[2]
        sellen[i,,3] <- 1/(1+19^{-(mulbins-s50)/(s95-s50)})
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,3] <- sum(seltmp*dl)
        
        # 4. NMRL (logistic)
        
        s50 <- selparsNMRL[1]
        s95 <- s50+selparsNMRL[2]
        sellen[i,,4] <- 1/(1+19^{-(mulbins-s50)/(s95-s50)})
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,4] <- sum(seltmp*dl)
        
        # 5. SMRL (logistic)
        
        s50 <- selparsSMRL[1]
        s95 <- s50+selparsSMRL[2]
        sellen[i,,5] <- 1/(1+19^{-(mulbins-s50)/(s95-s50)})
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,5] <- sum(seltmp*dl)
      }
    }
    
    # initial population dynamics
    
    Rec[i,1] <- R0
    for(s in 1:ns) 
      for(r in 1:nr) N[1,1,s,r] <- R0 * eta[r] * zeta[s]
    for(a in 2:(nages-1))
      for(s in 1:ns)
        for(r in 1:nr) N[1,a,s,r] <- sum(pim[,r]*N[1,a-1,s,]*exp(-M))
    
    av <- rep(NA,nr)
    B <- matrix(nrow=nr,ncol=nr)
    for(s in 1:ns) {
      
      for(r in 1:nr) av[r] <- sum(pim[,r]*N[1,nages-1,s,]*exp(-M))
      for(r in 1:nr) {
        for(ss in 1:nr) { 
          B[ss,r] <- pim[ss,r]*exp(-M)
        }
      }
      
      npg <- t(av) %*% solve(diag(1,nrow=nr)-B)
      N[1,nages,s,] <- npg
        
    }
    
    B0 <- sum(apply(N[1,,1,],1,sum)*w[,1]*m[,1])
    SSBf[i,1] <- B0
    for(s in 1:ns)
      for(r in 1:nr) SSB[i,1,s,r] <- sum(N[1,,s,r]*w[,s]*m[,s])
    rho <- B0/R0
    alp <- 4*hh/(rho*(1-hh))
    bet <- (5*hh-1)/(B0*(1-hh))
    
    # initial harvest rates
    
    for(f in 1:nf) {
      
      Xtmp <- 0
      for(s in 1:ns) Xtmp <- Xtmp+sum(N[1,,s,rff[[f]]]*exp(-tau[1,f]*M)*sel[,s,f]*w[,s])
      hyf[1,f] <- C[1,f]/Xtmp
      
      for(s in 1:ns) hyfa[1,,s,f] <- hyf[1,f] * sel[,s,f]
    }
    
    # if there's catch in first year do hmax checks
      
    if(any(C[1,] > 0.01)) {
        
        for(s in 1:ns) {
          for(a in 1:nages) {
            for(r in 1:nr) {
        
              fsub <- names(rff[rff == r])
              csum <- sum(C[1,fsub])
              hsum <- 0
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                hsum <- hsum+sum(hyfa[1,a,s,fx])
              }
              
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                if(hsum > hmax) hyfa[1,a,s,fx] <- hmax*C[1,fx]/csum
            }
          }
        }
      }
    }
    
    # total harvest rates by region
    
    for(r in 1:nr) {
     
      nfr <- length(rff[rff == r]) 
      fv <- rep(NA,nfr)
      for(f in 1:nfr) {
        
        fsub <- names(rff[rff == r])
        fv[f] <- grep(fsub[f],ffnm)
        
      }
      
      hsub <- hyfa[1,,,fv]
      hya[1,,,r] <- apply(hsub,c(1,2),sum)
    }
    
    # annual loop
    
    for(y in 2:ny) {
      
      # recruitment (different for years where recruitment estimated)
    
      Rtot <- alp*SSBf[i,y-1]/(1+bet*SSBf[i,y-1])
      for(s in 1:ns) {
        for(r in 1:nr) {
          
          if(y < recymin | y > recymax) {
            
            N[y,1,s,r] <- Rtot * eta[r] * zeta[s]
            
          } else {
            
            N[y,1,s,r] <- Rtot * eta[r] * zeta[s] *exp(epsR[y-recymin+1]-sigmaR^2/2)
          }
        }
      }
       
      Rec[i,y] <- sum(N[y,1,,]) 
      # post-recruit dynamics
      
      # a = 2 to A+-1
      
      for(s in 1:ns) 
        for(a in 2:(nages-1)) 
          for(r in 1:nr) N[y,a,s,r] <- sum(pim[,r]*N[y-1,a-1,s,]*exp(-M)*(1-hya[y-1,a-1,s,]))
    
      # plus group
          
      for(s in 1:ns)
        for(r in 1:nr) N[y,nages,s,r] <- sum(pim[,r]*N[y-1,nages-1,s,]*exp(-M)*(1-hya[y-1,nages-1,s,]))+sum(pim[,r]*N[y-1,nages,s,]*exp(-M)*(1-hya[y-1,nages,s,]))
                          
      # SSB
      
      SSBf[i,y] <- sum(apply(N[y,,1,],1,sum)*w[,1]*m[,1])
      for(s in 1:ns)
        for(r in 1:nr) SSB[i,y,s,r] <- sum(N[y,,s,r]*w[,s]*m[,s])
      
      # harvest rates
     
      for(f in 1:nf) {
        
        Xtmp <- 0
        for(s in 1:ns) Xtmp <- Xtmp+sum(N[y,,s,rff[[f]]]*exp(-tau[y,f]*M)*sel[,s,f]*w[,s])
        hyf[y,f] <- C[y,f]/Xtmp
        
        for(s in 1:ns) hyfa[y,,s,f] <- hyf[y,f] * sel[,s,f]
      }
      
      # if there's catch in first year do hmax checks
      
      if(any(C[y,] > 0.01)) {
        
        for(s in 1:ns) {
          for(a in 1:nages) {
            for(r in 1:nr) {
              
              fsub <- names(rff[rff == r])
              csum <- sum(C[y,fsub])
              hsum <- 0
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                hsum <- hsum+sum(hyfa[y,a,s,fx])
              }
              
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                if(hsum > hmax) hyfa[y,a,s,fx] <- hmax*C[y,fx]/csum
              }
            }
          }
        }
      }
      
      # total harvest rates by region
      
      for(r in 1:nr) {
        
        nfr <- length(rff[rff == r]) 
        fv <- rep(NA,nfr)
        for(f in 1:nfr) {
          
          fsub <- names(rff[rff == r])
          fv[f] <- grep(fsub[f],ffnm)
          
        }
        
        hsub <- hyfa[y,,,fv]
        hya[y,,,r] <- apply(hsub,c(1,2),sum)
      } 
    }
    
    Nfin[i,,,] <- N[ny,,,]
    hyafin[i,,,] <- hya[ny,,,]

    if(i %% 50 == 0) cat("iteration",i,"of",nits,"\n")
    
  }
  
  return(list(SSB=SSB,SSBf=SSBf,Rec=Rec,Nfin=Nfin,hyafin=hyafin,sellen=sellen))
}

# v2 projection code

proj2 <- function(nits = 1,nyp = 1,parm,parv,parnm1,parnm2,B0,nfin,hyafin,taufin,tacsplit,recdev='parametric',mat.type='new') {
  
  a.lw <- 4.4e-6
  b.lw <- 3.14
  if(mat.type == 'old') {
    
    ml50 <- 139.6
    ml95 <- 185.8
    
  }
  if(mat.type == 'new') {
    
    nx <- data$matpars[1]
    mux <- data$matpars[2]
    
  }
  
  # arrays
  
  N <- array(dim=c(nyp,na,ns,nr))
  Rec <- array(dim=c(nits,nyp))
  SSBf <- Rec <- matrix(nrow=nits,ncol=nyp) 
  hyf <- array(dim=c(nyp,nf))
  hyfa <- array(dim=c(nyp,na,ns,nf))
  hya <- array(dim=c(nyp,na,ns,nr))
  sel <- array(dim=c(na,ns,nf))
  m <- w <- array(dim=c(na,ns))
  mula <- array(dim=c(na,ns))
  pla <- array(dim=c(ns,nl,na))
  
  # definitely non-MCMC stuff
  
  hh <- data$hh
  M <- data$M
  zeta <- data$zeta
  sigmaR <- data$sigmaR
  hmax <- data$hmax
  nages <- length(ages)
  rff <- rf+1
  Cf <- tacsplit
  
  for(i in 1:nits) {
    
    # parameters
    
    R0 <- ifelse(any(parnm1 == 'logR0'),exp(parm[i,grep('logR0',colnames(parm))]),exp(parv[['logR0']]))
    if(any(parnm1 == 'epsR')) {
      epsR <- unname(parm[i,grep('epsR',colnames(parm))])
    } else {
      epsR <- unname(parv[['epsR']])
    }
    eta1 <- ifelse(any(parnm1 == 'ilogeta'),ilogit(parm[i,grep('ilogeta',colnames(parm))]),ilogit(parv[['ilogeta']])) 
    eta <- c(eta1,1-eta1)  
    if(any(parnm1 == 'ilogpim')) {
      pimv <- unname(ilogit(parm[i,grep('ilogpim',colnames(parm))]))
    } else {
      pimv <- unname(ilogit(parv[['ilogpim']]))
    }
    pim <- matrix(c(pimv[1],1-pimv[1],1-pimv[2],pimv[2]),nrow=2,ncol=2,byrow=T)
    if(any(parnm1 == 'selparsATT')) {
      selparsATT <- unname(exp(parm[i,grep('selparsATT',colnames(parm))]))
    } else {
      selparsATT <- unname(exp(parv[['selparsATT']]))
    }
    if(any(parnm1 == 'selparsNVT')) {
      selparsNVT <- unname(exp(parm[i,grep('selparsNVT',colnames(parm))]))
    } else {
      selparsNVT <- unname(exp(parv[['selparsNVT']]))
    }
    if(any(parnm1 == 'selparsATL')) {
      selparsATL <- unname(exp(parm[i,grep('selparsATL',colnames(parm))]))
    } else {
      selparsATL <- unname(exp(parv[['selparsATL']]))
    }
    if(any(parnm1 == 'selparsNMRL')) {
      selparsNMRL <- unname(exp(parm[i,grep('selparsNMRL',colnames(parm))]))
    } else {
      selparsNMRL <- unname(exp(parv[['selparsNMRL']]))
    }
    if(any(parnm1 == 'selparsSMRL')) {
      selparsSMRL <- unname(exp(parm[i,grep('selparsSMRL',colnames(parm))]))
    } else {
      selparsSMRL <- unname(exp(parv[['selparsSMRL']]))
    }
    sigmaR <- ifelse(any(parnm1 == 'lsigmaR'),exp(parm[i,grep('lsigmaR',colnames(parm))]),exp(parv[['lsigmaR']]))
    
    
    # growth handled separately now

    l1 <- data$l1
    l2 <- data$l2
    k <- data$k
    sdla <- data$sdla

    # mean length-at-age
    
    mula[,1] <- l1[1]+(l2[1]-l1[1])*(1-exp(-k[1]*(ages-aref[1])))/(1-exp(-k[1]*(aref[2]-aref[1])))
    mula[,2] <- l1[2]+(l2[2]-l1[2])*(1-exp(-k[2]*(ages-aref[1])))/(1-exp(-k[2]*(aref[2]-aref[1])))
    
    # P(l | a)
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        psum <-sum(dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE))   
        pla[s,,a] <- dnorm(log(mulbins),log(mula[a,s]),sdla[s],FALSE)/psum
        
      } 
    }
    
    # estimate weight, maturity and selectivity at age from length
    
    for(s in 1:ns) {
      for(a in 1:na) {
        
        # weights dl
        
        ll <- mula[a,s]*exp(-1.96*sdla[s])
        lu <- mula[a,s]*exp(1.96*sdla[s])
        ldel <- (lu-ll)/25
        lref <- ll+(0:24)*ldel
        dl <- dnorm(log(lref),log(mula[a,s]),sdla[s],FALSE)
        dl <- dl/sum(dl)
        
        # maturity and weight-at-age
    
        if(mat.type == 'old') ml <- 1/(1+19^{-(lref-ml50)/(ml95-ml50)})
        if(mat.type == 'new') ml <- lref^nx/(mux^nx+lref^nx)
        wl <- a.lw*lref^b.lw
        m[a,s] <- sum(dl*ml)
        w[a,s] <- sum(dl*wl)/1e+3
        
        ###############
        # selectivity #
        ###############
        
        # 1. ATT (double logistic with dr95 fixed at 3)
        
        sl50 <- selparsATT[1]
        dl95 <- selparsATT[2]
        dr95 <- 3
        sr50 <- sl50+dl95+dr95+selparsATT[3]
        seltmp <- (1+19^{-(lref-sl50)/dl95})^{-1}*(1-(1+19^{-(lref-sr50)/dr95})^{-1})
        sel[a,s,1] <- sum(seltmp*dl)
        
        # 2. NVT (gamma with z fixed at 1)
        
        ga <- selparsNVT[1]
        gb <- selparsNVT[2]
        gmax <- (ga/gb)^ga * exp(-ga)
        seltmp <- lref^ga * exp(-gb*lref)/gmax
        sel[a,s,2] <- sum(seltmp*dl)
        
        # 3. ATL (logistic)
        
        s50 <- selparsATL[1]
        s95 <- s50+selparsATL[2]
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,3] <- sum(seltmp*dl)
        
        # 4. NMRL (logistic)
        
        s50 <- selparsNMRL[1]
        s95 <- s50+selparsNMRL[2]
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,4] <- sum(seltmp*dl)
        
        # 5. SMRL (logistic)
        
        s50 <- selparsSMRL[1]
        s95 <- s50+selparsSMRL[2]
        seltmp <- 1/(1+19^{-(lref-s50)/(s95-s50)})
        sel[a,s,5] <- sum(seltmp*dl)
      }
    }
    
    # initial numbers-at-age
    
    N[1,,,] <- nfin[i,,,]
    hya[1,,,] <- hyafin[i,,,]
    SSBf[i,1] <- sum(apply(N[1,,1,],1,sum)*w[,1]*m[,1])
    Rec[i,1] <- sum(N[1,1,,])
    
    # project forwards
    
    for(y in 2:nyp) {
      
      # recruitment (different for years where recruitment estimated)
      
      rho <- B0[i]/R0
      alp <- 4*hh/(rho*(1-hh))
      bet <- (5*hh-1)/(B0[i]*(1-hh))
      Rtot <- alp*SSBf[i,y-1]/(1+bet*SSBf[i,y-1])
      for(s in 1:ns) {
        for(r in 1:nr) {
          
          if(recdev == 'parametric') N[y,1,s,r] <- Rtot * eta[r] * zeta[s] * rlnorm(1,-sigmaR^2/2,sigmaR)
          if(recdev == 'nonparametric') {
             
            epsr <- sample(epsR,1,replace=T)
            N[y,1,s,r] <- Rtot * eta[r] * zeta[s] * exp(epsr-sigmaR^2/2)
            
          }
        }
      }
      Rec[i,y] <- sum(N[y,1,,]) 
      
      # post-recruit dynamics
      
      # a = 2 to A+-1
      
      for(s in 1:ns) 
        for(a in 2:(nages-1)) 
          for(r in 1:nr) N[y,a,s,r] <- sum(pim[,r]*N[y-1,a-1,s,]*exp(-M)*(1-hya[y-1,a-1,s,]))
      
      # plus group
      
      for(s in 1:ns)
        for(r in 1:nr) N[y,nages,s,r] <- sum(pim[,r]*N[y-1,nages-1,s,]*exp(-M)*(1-hya[y-1,nages-1,s,]))+sum(pim[,r]*N[y-1,nages,s,]*exp(-M)*(1-hya[y-1,nages,s,]))
      
      # SSB
      
      SSBf[i,y] <- sum(apply(N[y,,1,],1,sum)*w[,1]*m[,1])
      #for(s in 1:ns)
      #  for(r in 1:nr) SSB[i,y,s,r] <- sum(N[y,,s,r]*w[,s]*m[,s])
      
      # harvest rates
      
      for(f in 1:nf) {
        
        Xtmp <- 0
        for(s in 1:ns) Xtmp <- Xtmp+sum(N[y,,s,rff[[f]]]*exp(-taufin[f]*M)*sel[,s,f]*w[,s])
        hyf[y,f] <- Cf[f]/Xtmp
        
        for(s in 1:ns) hyfa[y,,s,f] <- hyf[y,f] * sel[,s,f]
      }
      
      # if there's catch in first year do hmax checks
      
      if(any(Cf > 0.01)) {
        
        for(s in 1:ns) {
          for(a in 1:nages) {
            for(r in 1:nr) {
              
              fsub <- names(rff[rff == r])
              csum <- sum(Cf[fsub])
              hsum <- 0
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                hsum <- hsum+sum(hyfa[y,a,s,fx])
              }
              
              for(f in 1:length(fsub)) {
                
                fx <- grep(fsub[f],ffnm)
                if(hsum > hmax) hyfa[y,a,s,fx] <- hmax*Cf[fx]/csum
              }
            }
          }
        }
      }
      
      # total harvest rates by region
      
      for(r in 1:nr) {
        
        nfr <- length(rff[rff == r]) 
        fv <- rep(NA,nfr)
        for(f in 1:nfr) {
          
          fsub <- names(rff[rff == r])
          fv[f] <- grep(fsub[f],ffnm)
          
        }
        
        hsub <- hyfa[y,,,fv]
        hya[y,,,r] <- apply(hsub,c(1,2),sum)
      } 
    }
   
    if(i %% 50 == 0) cat("iteration",i,"of",nits,"\n")
  }
  
  return(list(SSBf=SSBf,Rec=Rec))
}
