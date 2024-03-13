######################################################
# testing MR estimators for toothie MSE ##############
######################################################
# R. Hillary CSIRO 2024 ##############################
######################################################

library(TMB)
library(parallel)
source("test_sumto1.R")

load("filestore/dep50.rda")

logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

# load up test data set

tdf <- read.table("test_tags.dat",header=TRUE)

## 1. aggregated across sex, age, and space

compile("tagmod1.cpp")
dyn.load(dynlib("tagmod1"))

ttdf <- aggregate(tdf$nrec,by=list(recev=tdf$recev,yrel=tdf$yrel),FUN=sum)

# sorting out the nrel "issue"

ntdf <- aggregate(tdf$nrel,by=list(recev=tdf$recev,yrel=tdf$yrel),FUN=sum)
ntdf <- aggregate(ntdf$x,by=list(yrel=ntdf$yrel),FUN=function(x){x <- mean(x)/2})
nrdf <- aggregate(ttdf$x,by=list(yrel=ttdf$yrel),FUN=sum)

relyr <- ntdf$yrel
nrel <- ntdf$x
nrelev <- length(nrel)
ntot <- nrdf$x
nrecev <- aggregate(ttdf$recev,by=list(yrel=ttdf$yrel),FUN=max)$x
nrecmax <- max(nrecev)

R <- array(dim=c(nrelev,nrecmax+2))
R[,1] <- nrel
R[,2] <- ntot
for(y in 1:nrelev) R[y,3:(2+nrecev[y])] <- subset(ttdf,yrel==relyr[y])$x

nytot <- max(relyr+nrecev)-relyr[1]+1
dms <- c(nytot,nrelev,max(nrecev))
pret <- p1tag(2,1:nrecmax,psi,nu)

data <- list(dms=dms,
             relyr=relyr-1,
             nrecev=nrecev,
             R=R,
             M=0.13,
             prep=prep,
             pret=pret,
             sigmah=0.1)

pars <- list(muih=logit(0.03),
             epshy=rep(0,nytot))

obj1 <- MakeADFun(data=data,parameters=pars,DLL="tagmod1")
obj1$fn()

# run the mofo

res1 <- do.call(optim,obj1)
rep1 <- obj1$rep()

# RE version

compile("tagmod1re.cpp")
dyn.load(dynlib("tagmod1re"))

parsre <- list(muih=logit(0.03),lnsigmah=log(0.1),epshy=rep(0,nytot))

obj1re <- MakeADFun(data=data,parameters=parsre,random=c("epshy"),DLL="tagmod1re")

res1re <- do.call(optim,obj1re)
rep1re <- obj1re$rep()
sdrep1 <- sdreport(obj1re)
summary(sdrep1)

## 2. aggregated across sex and age

compile("tagmod2.cpp")
dyn.load(dynlib("tagmod2"))

# prep the data with spatial recapture info 

tsdf <- aggregate(tdf$nrec,by=list(recev=tdf$recev,yrel=tdf$yrel,rrel=tdf$rrel,recr=tdf$recr),FUN=sum)

# sorting out the nrel "issue"

ntsdf <- aggregate(tdf$nrel,by=list(recev=tdf$recev,yrel=tdf$yrel,rrel=tdf$rrel),FUN=sum)
ntsdf <- aggregate(ntsdf$x,by=list(yrel=ntsdf$yrel,rrel=ntsdf$rrel),FUN=function(x){x <- mean(x)/2})
nrsdf <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,rrel=tsdf$rrel),FUN=sum)
relsyr <- ntsdf$yrel
nsrel <- ntsdf$x
nsrelev <- length(nsrel)
nstot <- nrsdf$x
nsrecev <- aggregate(tsdf$recev,by=list(yrel=tsdf$yrel,rrel=tsdf$rrel),FUN=max)$x
nsrecmax <- max(nsrecev)

# modified tag retention

psret <- p1tag(2,0:(nrecmax-1),psi,nu)

# data object construction

dmss <- c(nytot,length(unique(tsdf$rrel)),nsrelev,max(nsrecev))

ysrel <- ntsdf$yrel
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

data2 <- list(dms=dmss,
              relyr=ysrel-1,
              relr=nrsdf$rrel-1,
              nrecev=nsrecev,
              T=ntsdf$x,
              NR=nrsdf$x,
              R=RS,
              M=0.13,
              prep=prep,
              pret=psret,
              sigmah=0.1)

xphi <- inv.add.logit(Biol$T)
pars2 <- list(muih=logit(rep(0.03,nr)),
              epshy=matrix(0,nrow=nytot,ncol=nr),
              xphi=xphi)

obj2 <- MakeADFun(data=data2,parameters=pars2,DLL="tagmod2")
obj2$fn()
rep2 <- obj2$rep()

res2 <- do.call(optim,obj2)
rep2 <- obj2$rep()

# RE version of spatial estimator

compile("tagmod2re.cpp")
dyn.load(dynlib("tagmod2re"))

pars2re <- list(muih=logit(rep(0.03,nr)),
                lnsigmah=log(0.1),
                epshy=matrix(0,nrow=nytot,ncol=nr),
                xphi=xphi)

obj2re <- MakeADFun(data=data2,parameters=pars2re,random=c("epshy"),DLL="tagmod2re")

res2re <- do.call(optim,obj2re)
rep2re <- obj2re$rep()
sdrep2 <- sdreport(obj2re)
summary(sdrep2)
