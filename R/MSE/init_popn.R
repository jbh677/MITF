######################################################
# initial popn of OM #################################
######################################################
# R. Hillary CSIRO 2023 ##############################
######################################################

library(parallel)
library(abind)
library(MCMCpack)
library(TMB)
source("mse_utils.R")                  

###########
# Biology #
###########

Biol <- list()

# dimensions

na <- 50
nr <- 2 # spatial regions
Biol$dms <- c(na,nr)
Biol$ages <- 1:na
Biol$R0 <- 1e+6
Biol$M <- 0.13
Biol$hh <- 0.75
Biol$eta <- rep(1/nr,nr) # relative spatial recruitment
sm50 <- 110
sm95 <- 150
a.wl <- 4.4e-6
b.wl <- 3.14
linf <- c(170,170)
k <- c(0.045,0.045)
t0 <- c(-1.5,-1.5)
sdla <- c(0.1,0.1)
Biol$lbins <- seq(20,190,by=10)
Biol$mulbins <- 0.5*(Biol$lbins[-1]+Biol$lbins[-length(Biol$lbins)])
Biol$mpars <- c(sm50,sm95)
Biol$gpars <- list(linf=linf,k=k,t0=t0,sdla=sdla)
Biol$wpars <- c(a.wl,b.wl)
Biol <- get.matwt(Biol)
Biol <- get.pla(Biol)
pstick <- 0.9 # fraction staying in same region
pmove <- (1-pstick)/(nr-1)
T <- matrix(nrow=nr,ncol=nr)
for(i in 1:nr)
  for(j in 1:nr) T[i,j] <- ifelse(i == j,pstick,pmove)
Biol$T <- T

##########
# Fleets #
##########

Fleet <- list()

# dimensions

nf <- 2
Fleet$nf <- nf
Fleet$fref <- c(1,2)
# double normal parameterisation:
# 1. smax; 2. sl; 3. sr
selpars <- array(dim=c(3,nf))
selpars[,] <- c(80,2,50)
Fleet$selpars <- selpars
Fleet <- get.sel(Fleet,Biol)

########################################
# initial population & fishery targets #
########################################

# overall female SSB depletion

deltarg <- 0.5
rhotarg <- (deltarg*(5*Biol$hh-1)+1-Biol$hh)/(4*Biol$hh)

# target relative catch by fishery 

pctarg <- rep(1/nf,nf)

# target vector

targv <- logit(c(deltarg,pctarg))

theta.init <- logit(c(0.04,0.04))
res.init <- optim(theta.init,objfn.init,method=c("L-BFGS-B"),control=list(trace=1))
hinit <- ilogit(res.init$par)
popinit <- get.exploited.eqm(Biol,Fleet,hinit)

# wee check

tmp <- get.unexploited.eqm(Biol)
B0 <- Biol$R0*tmp$rho
hh <- Biol$hh
alp <- 4*hh/(tmp$rho*(1-hh))
bet <- (5*hh-1)/(B0*(1-hh))
alp*B0/(1+bet*B0) # better equal R0
Biol$R0

########################
# move forward in time #
########################

#load("filestore/dep50.rda")

# total number of years

nyproj <- 20                           

# additional biology

Biol$sigmar <- 0.4
Biol$alp <- alp
Biol$bet <- bet
Biol$Ninit <- tmp$N # start at unexploited eqm

# additional fleet information

Fleet$hinit <- hinit
Fleet$hmult <- rep(1,nyproj)

# type of control

CatchControl <- FALSE
Fleet$Cfix <- c(240,240)

# iterations

nits <- 1000

# multi-core

nitsx <- 100
ncore <- 10

system.time(prj <- get.popdyn(Biol,Fleet,nyproj,CatchControl,nits)) 

boxplot(t(prj$SSBtot),outline=FALSE,col='magenta',ylim=c(0,B0))


system.time(prj <- proj.para(Biol,Fleet,nyproj,CatchControl,nitsx,ncore,TRUE))
#system.time(prj <- proj.para(Biol,Fleet,nyproj,CatchControl,nitsx,ncore,FALSE))

boxplot(t(prj$SSBtot)/B0,outline=FALSE,col='magenta',ylim=c(0,1.05))
abline(h=0.5,lty=2)
par(mfrow=c(nf,1))
for(f in 1:nf) boxplot(t(prj$C[,f,]),outline=FALSE,col='magenta',ylim=c(0,400))

################
# tag dynamics #
################

# tags per tonne of catch

ntpt <- 3

# tag mortality (probability of surviving tag release)

ptmort <- 1

# tag shedding parameters

psi <- 0.9 # instantaneous retention probability
nu <- 0.015 # annual rate of shedding post loss
qt <- function(tau,psi,nu) {return(psi*exp(-nu*tau))}
# P(ntag = 2, 1, 0 | tau) = {qt^2, 2*qt*(1-qt), (1-qt)^2}

# tag reporting rate (different for number of attached tags)

prep <- 0.95

# create tag object

Tag <- list(ntpt=ntpt,prep=prep,ptmort=ptmort,mulrel=log(72),sdlrel=0.3,ytagon=10,ytagoff=nyproj-1,maxrecev=7,phi=1.05,xi=0.25)

# testing 

xxx <- list(N=prj$N[,,,,1],H=prj$H[,,,,1],C=prj$C[,,1])

system.time(tdf <- get.tag.data2(xxx))

# generate tag data using parallel functionality

prjlist <- listify.proj(prj)
system.time(tdflist <- mclapply(1:nits,wrap.get.tag.data2,mc.cores=ncore))

mctrl <- list(M=Biol$M,sigmah=0.05)
model <- 'tgm1re'
system.time(replist <- mclapply(1:nits,wrap.get.hrates,mc.cores=ncore))

hyhat <- array(dim=c(length(replist[[1]]$hy),nits))
for(nn in 1:nits) hyhat[,nn] <- replist[[nn]]$hy
hbar <- mean(unlist(lapply(prjlist,function(x){mean(x$H[1:dim(hyhat)[1],,,])})))
boxplot(t(hyhat),outline=FALSE,col='green',ylim=c(0,0.05))
abline(h=hbar,lty=2)

# spatial options

model <- 'tgm2re'
system.time(replist <- mclapply(1:nits,wrap.get.hrates,mc.cores=ncore))

hyshat <- array(dim=c(dim(replist[[1]]$hy),nits))
for(nn in 1:nits) hyshat[,,nn] <- replist[[nn]]$hy
hbar <- mean(unlist(lapply(prjlist,function(x){mean(x$H[1:dim(hyhat)[1],,,])})))
par(mfrow=c(2,1))
boxplot(t(hyshat[,1,]),outline=FALSE,col='green',ylim=c(0,0.05))
abline(h=hbar,lty=2)
boxplot(t(hyshat[,2,]),outline=FALSE,col='green',ylim=c(0,0.05))
abline(h=hbar,lty=2)

# spatial estimates

Phi.mc <- array(dim=c(nits,nr,nr))
for(nn in 1:nits) Phi.mc[nn,,] <- replist[[nn]]$Phi
apply(Phi.mc,c(2,3),quantile,c(0.025,0.5,0.975))
Biol$T

#######################
# running a simple MP #
#######################

# MSE control file

mctrl <- list(M=Biol$M,sigmah=0.05)
msectrl <- list(initmod='tgm1re',prjmod='tgm1re',mctrl=mctrl,yrng=1:25,nyinit=20,firstdec=20,tactau=5)

# MP definition

tacctrl <- list(MP='fixedh',ytau=4,htarg=0.026,maxChange=0.2,catchSplit=rep(1/nf,nf))

save.image(file='filestore/dep50.rda')

