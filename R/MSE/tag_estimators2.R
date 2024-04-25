######################################################
# testing MR estimators for toothie MSE ##############
######################################################
# R. Hillary CSIRO 2024 ##############################
######################################################

library(TMB)
library(parallel)
source("mse_utils.R")

load("filestore/dep50.rda")

logit <- function(x){return(log(x/(1-x)))}
ilogit <- function(x){return(1/(1+exp(-x)))}

# load up test data set

tdf <- read.table("test_tags2.dat",header=TRUE)

## 1. aggregated across sex, age, and space

compile("tagmod1.cpp")
dyn.load(dynlib("tagmod1"))

t2df <- subset(tdf,recev>1)
ttdf <- aggregate(t2df$R,by=list(recev=t2df$recev,yrel=t2df$yrel),FUN=sum)

# sorting out the nrel "issue"

t3df <- subset(tdf,recev==1)
ntdf <- aggregate(t3df$T-t3df$R,by=list(yrel=t3df$yrel),FUN=sum)
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

data <- list(dms=dms,
             relyr=relyr-1,
             nrecev=nrecev,
             R=R,
             M=0.13,
             prep=rep(prep,nytot),
             pret=pret,
             sigmah=0.25)

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
ntsdf <- aggregate(ntsdf$x,by=list(yrel=ntsdf$yrel,rrel=ntsdf$rrel),FUN=function(x){x <- mean(x)})
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
              prep=rep(prep,nytot),
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

#############################################
# running on the actual Macca 2021 tag data #
#############################################

mtag <- read.table("tags21.dat",header=TRUE)
maxrecev <- 9
mtag <- subset(mtag,recyr-relyr<=maxrecev)

# aggregate across space and length

midf <- aggregate(mtag$R,by=list(yrel=mtag$relyr,yrec=mtag$recyr),FUN=sum)

# maximum number of recapture events 

tdf <- aggregate(mtag$T,by=list(yrel=mtag$relyr,yrec=mtag$recyr),FUN=sum)
ntdf <- aggregate(tdf$x,by=list(yrel=tdf$yrel),FUN=function(x){x <- mean(x)/2})
nrdf <- aggregate(midf$x,by=list(yrel=midf$yrel),FUN=sum)

relyr <- ntdf$yrel
nrel <- ntdf$x
nrelev <- length(nrel)
ntot <- nrdf$x
nrecev <- aggregate(tdf$yrec-tdf$yrel,by=list(yrel=tdf$yrel),FUN=max)$x
nrecmax <- max(nrecev)

R <- array(dim=c(nrelev,nrecmax+2))
R[,1] <- nrel
R[,2] <- ntot
for(y in 1:nrelev) R[y,3:(2+nrecev[y])] <- subset(midf,yrel==relyr[y])$x

nytot <- max(relyr+nrecev)-relyr[1]+1
dms <- c(nytot,nrelev,max(nrecev))
pret <- p1tag(2,1:nrecmax,psi,nu)
prepdf <- read.table("tagrep2021.dat",header=TRUE)

data <- list(dms=dms,
             relyr=relyr-relyr[1],
             nrecev=nrecev,
             R=R,
             M=0.13,
             prep=subset(prepdf,Year>=min(midf$yrec)-1 & Year<=max(midf$yrec))$RR,
             pret=pret,
             sigmah=0.3)

pars <- list(muih=logit(0.05),
             epshy=rep(0,nytot))

obj1 <- MakeADFun(data=data,parameters=pars,DLL="tagmod1")
obj1$fn()

## v1

res1 <- do.call(optim,obj1)
rep1 <- obj1$rep()

## v1re

parsre <- list(muih=logit(0.05),lnsigmah=log(0.5),epshy=rep(0,nytot))

obj1re <- MakeADFun(data=data,parameters=parsre,random=c("epshy"),DLL="tagmod1re")

res1re <- do.call(optim,obj1re)
rep1re <- obj1re$rep()
sdrep1 <- sdreport(obj1re)
summary(sdrep1)

## v2

tsdf <- aggregate(mtag$R,by=list(yrel=mtag$relyr,recyr=mtag$recyr,relr=mtag$relarea,recr=mtag$recarea),FUN=sum)
ntsdf <- aggregate(mtag$T,by=list(yrel=mtag$relyr,recyr=mtag$recyr,rrel=mtag$relarea),FUN=sum)
ntsdf <- aggregate(ntsdf$x,by=list(yrel=ntsdf$yrel,rrel=ntsdf$rrel),FUN=function(x){x <- mean(x)/2})
nrsdf <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,rrel=tsdf$relr),FUN=sum)
relsyr <- ntsdf$yrel
nsrel <- ntsdf$x
nsrelev <- length(nsrel)
nstot <- nrsdf$x
nsrecev <- aggregate(tsdf$recyr-tsdf$yrel,by=list(yrel=tsdf$yrel,tsdf$relr),FUN=max)$x
nsrecmax <- max(nsrecev)

dmss <- c(nytot,length(unique(tsdf$relr)),nsrelev,max(nsrecev))

ysrel <- ntsdf$yrel
RS <- array(dim=c(nsrelev,nsrecmax,dmss[2]))
for(t in 1:nsrelev) {
  for(r in 1:dmss[2]) {

      yx <- ysrel[t]
      rx <- nrsdf$rrel[t]
      nrecx <- nsrecev[t]
      xdf <- subset(tsdf,yrel==yx & relr==rx)
      for(rr in 1:dmss[2]) RS[t,1:nrecx,rr] <- subset(xdf,recr==rr)$x
  }
}

psret <- p1tag(2,0:(nrecmax-1),psi,nu)

data2 <- list(dms=dmss,
              relyr=ysrel-ysrel[1],
              relr=nrsdf$rrel-1,
              nrecev=nsrecev,
              T=ntsdf$x,
              NR=nrsdf$x,
              R=RS,
              M=0.13,
              prep=subset(prepdf,Year>=min(midf$yrec)-1 & Year<=max(midf$yrec))$RR,
              pret=psret,
              sigmah=0.8)


Ttmp <- matrix(c(1-0.05,0.05,0.02,1-0.02),nrow=2,ncol=2,byrow=T)
xphi <- inv.add.logit(Ttmp)
pars2 <- list(muih=logit(c(0.05,0.01)),
              epshy=matrix(0,nrow=nytot,ncol=nr),
              xphi=xphi)

obj2 <- MakeADFun(data=data2,parameters=pars2,DLL="tagmod2")
obj2$fn()
rep2 <- obj2$rep()

res2 <- do.call(optim,obj2)
eep2 <- obj2$rep()
sdrep2 <- sdreport(obj2)
summary(sdrep2)

## v2re

pars2re <- list(muih=logit(c(0.05,0.01)),
                lnsigmah=log(0.8),
                epshy=matrix(0,nrow=nytot,ncol=nr),
                xphi=xphi)

obj2re <- MakeADFun(data=data2,parameters=pars2re,random=c("epshy"),DLL="tagmod2re")

res2re <- do.call(optim,obj2re)
rep2re <- obj2re$rep()
sdrep2 <- sdreport(obj2re)
summary(sdrep2)

############################
# fitting summary to Macca #
############################

library(ggplot2)

# spatially aggregated

Robs <- R[,-c(1,2)]                    
midf$xhat <- NA
yx <- midf$yrel-min(midf$yrel)+1
for(i in 1:dim(midf)[1]) {

  rex <- midf$yrec[i]-midf$yrel[i]  
  midf$xhat[i] <- rep1re$Rhat[yx[i],rex]

}

ggplot(midf)+geom_point(aes(x=yrec,y=x),colour='brown')+geom_point(aes(x=yrec,y=xhat),colour='dark green')+facet_wrap(~yrel)+theme_bw()+xlab("Year of recapture")+ylab("Recaptures")+ggtitle("Non-spatial: recaptures given releases")+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1)) 

# spatially explicit

RSobs <- data2$R
tsdf$xhat <- NA
for(i in 1:dim(ntsdf)[1]) {

  rx <- ntsdf$rrel[i]
  for(rr in 1:nr) {

    tsdf[tsdf$yrel==ntsdf$yrel[i] & tsdf$relr==rx & tsdf$recr==rr,'xhat'] <- rep2re$Rhat[i,1:nsrecev[i],rr] 

  }
}

# aggregated across release and recapture regions

tmp1 <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,recyr=tsdf$recyr),FUN=sum)
tmp2 <- aggregate(tsdf$xhat,by=list(yrel=tsdf$yrel,recyr=tsdf$recyr),FUN=sum)
tmp1 <- cbind(tmp1,tmp2$x)
names(tmp1)[[4]] <- "xhat"

ggplot(tmp1)+geom_point(aes(x=recyr,y=x),colour='brown')+geom_point(aes(x=recyr,y=xhat),colour='dark green')+facet_wrap(~yrel)+theme_bw()+xlab("Year of recapture")+ylab("Recaptures")+ggtitle("Spatial: recaptures given releases")+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

# aggregated across release region and recapture year

tmp1 <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,recr=tsdf$recr),FUN=sum)
tmp2 <- aggregate(tsdf$xhat,by=list(yrel=tsdf$yrel,recr=tsdf$recr),FUN=sum)
tmp1 <- cbind(tmp1,tmp2$x)
names(tmp1)[[4]] <- "xhat"

ggplot(tmp1)+geom_point(aes(x=yrel,y=x),colour='brown')+geom_point(aes(x=yrel,y=xhat),colour='dark green')+facet_wrap(~recr)+theme_bw()+xlab("Year of release")+ylab("Recaptures")+ggtitle("Spatial: year of release and recapture region")+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

# aggregated across release region and release year

tmp1 <- aggregate(tsdf$x,by=list(yrec=tsdf$recyr,recr=tsdf$recr),FUN=sum)
tmp2 <- aggregate(tsdf$xhat,by=list(yrec=tsdf$recyr,recr=tsdf$recr),FUN=sum)
tmp1 <- cbind(tmp1,tmp2$x)
names(tmp1)[[4]] <- "xhat"

ggplot(tmp1)+geom_point(aes(x=yrec,y=x),colour='brown')+geom_point(aes(x=yrec,y=xhat),colour='dark green')+facet_wrap(~recr)+theme_bw()+xlab("Year of recapture")+ylab("Recaptures")+ggtitle("Spatial: year of recapture and recapture region")+theme(axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))

# aggregated across release region

tmp1 <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,recyr=tsdf$recyr,recr=tsdf$recr),FUN=sum)
tmp2 <- aggregate(tsdf$xhat,by=list(yrel=tsdf$yrel,recyr=tsdf$recyr,recr=tsdf$recr),FUN=sum)
tmp1 <- cbind(tmp1,tmp2$x)
names(tmp1)[[5]] <- "xhat"

ggplot(tmp1)+geom_point(aes(x=recyr,y=x),colour='brown')+geom_point(aes(x=recyr,y=xhat),colour='dark green')+facet_grid(recr~yrel)+theme_bw()+xlab("Year of recapture")+ylab("Recaptures")

# aggregated across recapture region

tmp1 <- aggregate(tsdf$x,by=list(yrel=tsdf$yrel,recyr=tsdf$recyr,relr=tsdf$relr),FUN=sum)
tmp2 <- aggregate(tsdf$xhat,by=list(yrel=tsdf$yrel,recyr=tsdf$recyr,relr=tsdf$relr),FUN=sum)
tmp1 <- cbind(tmp1,tmp2$x)
names(tmp1)[[5]] <- "xhat"

ggplot(tmp1)+geom_point(aes(x=recyr,y=x),colour='brown')+geom_point(aes(x=recyr,y=xhat),colour='dark green')+facet_grid(relr~yrel)+theme_bw()+xlab("Year of recapture")+ylab("Recaptures")

# save it

save.image("filestore/testing_estimators.rda")

