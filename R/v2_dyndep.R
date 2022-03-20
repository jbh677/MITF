#####################################################
# v2_dyndep dynamic depletion ##############################
#####################################################
# R. Hillary, J. Day CSIRO 2021 #####################
#####################################################

require(TMB)

compile("v2_dyndep.cpp")
dyn.load("v2_dyndep.so")
source("v2_utils.R")

# load up the data

load("MIdata_maxrecev_9.rda")

##########################################
# set up the data and fixed input object #
##########################################

A <- matrix(scan("aerr_mat.dat"),ncol=na,byrow=T)
aref <- c(5,20) # reference ages for Schnute growth parameterisation
dms <- c(ny,na,nl,nf,nr) # key dims: y, a, l, f, r
datflag <- c(1,1,1) # LF, ALF, tags
recyr <- c(1,ny-5) # min and max years for recruitment estimation
nobs <- c(dim(LF)[1],dim(ALF)[1],dim(tagrec)[1],nrecmax) # number of obs.
#phiLF <- rep(2,nf) # fisheries-specifc OD factor for LF data
#phiT <- 1.5 # OD factor for tagging data
hh <- 0.75 # steepness
#sigmaR <- 0.25 # stk-rec variability
M <- 0.13 # natural mortality
zeta <- c(0.5,0.5) # birth sex-ratio
hmax <- 0.8 # maximum harvest-rate per age-class and area
prtotmax <- 0.9 # maximum total chance of recovering a tag release
psi <- 1 # immediate tag shedding probability
nu <- 0 # continuous rate of tag shedding
sdl <- 0.05 # random variability (SD) in growth increments
matdf <- subset(read.table("maturity.tab",header=T),Sex=='F')
matpars <- c(matdf$nu,matdf$mu)  # female maturity-at-length parameters
gdf <- read.table("growth.tab",header=T)
l1 <- c(subset(gdf,Sex=='F')$l1,subset(gdf,Sex=='M')$l1)
l2 <- c(subset(gdf,Sex=='F')$l2,subset(gdf,Sex=='M')$l2)
k <- c(subset(gdf,Sex=='F')$k,subset(gdf,Sex=='M')$k)
sdla <- c(subset(gdf,Sex=='F')$sdl,subset(gdf,Sex=='M')$sdl)
sdilogeta <- 1.5
zz <- rnorm(10000,0,sdilogeta)
round(quantile(1/(1+exp(-zz)),c(0.025,0.975)),2) # prior 95% CI
Gamma <- get.tmat(gdf,lbins) # female/male growth transition matrices
phiLFmin <- 1.05 # min over-dispersion
phiLFmax <- min(LF[,1])-0.05 #max over-dispersion
phiTmin <- 1.05
phiTmax <- min(tagrel)-0.05

data <- list(C=unname(C),
             tau=unname(tau),
             LFcov=LFcov,
             LF=LF,
             ALFcov=ALFcov,
             ALF=ALF,
             tagrel=tagrel,
             tagnrec=tagnrec,
             tagcov=tagcov,
             tagrec=tagrec,
             A=A,
             Gamma=Gamma,
             aref=aref,
             phiLFmin = phiLFmin,
             phiLFmax = phiLFmax,
             phiTmin = phiTmin,
             phiTmax = phiTmax, 
             dms=dms,
             datflag=datflag,
             recyr=recyr,
             nobs=nobs,
             rf=unname(rf),
             mulbins=mulbins,
             lbins=lbins,
             hh=hh,
             M=M,
             zeta=zeta,
             hmax=hmax,
             prtotmax=prtotmax,
             psi=psi,
             nu=nu,
             prep=rrdf$RR,
             sdl=sdl,
             matpars=matpars,
             l1=l1,
             l2=l2,
             k=k,
             sdla=sdla,
             sdilogeta=sdilogeta)

######################################
# set up the initial parameter guess #
######################################

R0 <- 8e+5
lsigmaR <- log(0.3)
ilogeta <- logit(0.3)
ilogpim <- logit(c(0.99,0.94))
selparsATT <- c(55,5,30) # dr95 is fixed at 3
selparsNVT <- c(20,0.2857)
selparsATL <- c(60,10)
selparsNMRL <- c(70,20)
selparsSMRL <- c(90,30)
xphiLF <- rep(log((2-phiLFmin)/(phiLFmax-2)),nf)
xphiT <- log((1.85-phiTmin)/(phiTmax-1.85))
nyrec <- length(recyr[1]:recyr[2])

pars1 <- list(logR0 = log(R0),
             epsR = rep(0,nyrec),
             lsigmaR = lsigmaR,
             ilogeta = ilogeta,
             ilogpim = ilogpim,
             selparsATT = log(selparsATT),
             selparsNVT = log(selparsNVT), 
             selparsATL = log(selparsATL), 
             selparsNMRL = log(selparsNMRL), 
             selparsSMRL = log(selparsSMRL),
             xphiLF = xphiLF,
             xphiT = xphiT)

###############################
# start the estimation phases #
###############################

# phase 1: R0

map1 <- list(epsR = rep(factor(NA),nyrec),
             lsigmaR = factor(NA), 
             ilogeta = factor(NA),
             ilogpim = rep(factor(NA),nr),
             selparsATT = rep(factor(NA),3), 
             selparsNVT = rep(factor(NA),2),
             selparsATL = rep(factor(NA),2), 
             selparsNMRL = rep(factor(NA),2), 
             selparsSMRL = rep(factor(NA),2), 
             xphiLF = rep(factor(NA),nf),
             xphiT = factor(NA)
             )

obj1 <- MakeADFun(data=data,parameters=pars1,map=map1,DLL="v2_dyndep")
res1 <- do.call("optim",obj1)
rep1 <- obj1$rep()
sdrep1 <- sdreport(obj1)
summary(sdrep1)

# phase 2: selectivity

map2 <- list(epsR = rep(factor(NA),nyrec),
             lsigmaR = factor(NA), 
             ilogeta = factor(NA),
             ilogpim = rep(factor(NA),nr),
             xphiLF = rep(factor(NA),nf),
             xphiT = factor(NA) )


pars2 <- list(logR0 = res1$par[['logR0']],
             epsR = rep(0,nyrec),
             lsigmaR = lsigmaR,
             ilogeta = ilogeta,
             ilogpim = ilogpim,
             selparsATT = log(selparsATT),
             selparsNVT = log(selparsNVT), 
             selparsATL = log(selparsATL), 
             selparsNMRL = log(selparsNMRL), 
             selparsSMRL = log(selparsSMRL),
             xphiLF = xphiLF,
             xphiT = xphiT
             )

obj2 <- MakeADFun(data=data,parameters=pars2,map=map2,DLL="v2_dyndep")
res2 <- do.call("optim",obj2)
rep2 <- obj2$rep()
sdrep2 <- sdreport(obj2)
summary(sdrep2)

# phase 3: spatial recruitment 

map3 <- list(epsR = rep(factor(NA),nyrec),
             lsigmaR = factor(NA), 
             ilogpim = rep(factor(NA),nr),
             xphiLF = rep(factor(NA),nf),
             xphiT = factor(NA) )


pars3 <- list(logR0 = res2$par[['logR0']],
             epsR = rep(0,nyrec),
             lsigmaR = lsigmaR,
             ilogeta = ilogeta,
             ilogpim = ilogpim,
             selparsATT = unname(res2$par[grep('selparsATT',names(res2$par))]),
             selparsNVT = unname(res2$par[grep('selparsNVT',names(res2$par))]), 
             selparsATL = unname(res2$par[grep('selparsATL',names(res2$par))]), 
             selparsNMRL = unname(res2$par[grep('selparsNMRL',names(res2$par))]), 
             selparsSMRL = unname(res2$par[grep('selparsSMRL',names(res2$par))]), 
             xphiLF = xphiLF,
             xphiT = xphiT 
             )

obj3 <- MakeADFun(data=data,parameters=pars3,map=map3,DLL="v2_dyndep")
res3 <- do.call("optim",obj3)
rep3 <- obj3$rep()
sdrep3 <- sdreport(obj3)
summary(sdrep3)

# phase 4: movement 

map4 <- list(epsR = rep(factor(NA),nyrec),
             lsigmaR = factor(NA), 
             xphiLF = rep(factor(NA),nf),
             xphiT = factor(NA) )


pars4 <- list(logR0 = res3$par[['logR0']],
             epsR = rep(0,nyrec),
             lsigmaR = lsigmaR, 
             ilogeta = res3$par[['ilogeta']],
             ilogpim = ilogpim,
             selparsATT = unname(res3$par[grep('selparsATT',names(res3$par))]),
             selparsNVT = unname(res3$par[grep('selparsNVT',names(res3$par))]), 
             selparsATL = unname(res3$par[grep('selparsATL',names(res3$par))]), 
             selparsNMRL = unname(res3$par[grep('selparsNMRL',names(res3$par))]), 
             selparsSMRL = unname(res3$par[grep('selparsSMRL',names(res3$par))]), 
             xphiLF = xphiLF,
             xphiT = xphiT 
             )

obj4 <- MakeADFun(data=data,parameters=pars4,map=map4,DLL="v2_dyndep")
res4 <- do.call("optim",obj4)
rep4 <- obj4$rep()
sdrep4 <- sdreport(obj4)
summary(sdrep4)

# phase 5: recruitment 

map5 <- list(lsigmaR = factor(NA), 
             xphiLF = rep(factor(NA),nf),
             xphiT = factor(NA) )

pars5 <- list(logR0 = res4$par[['logR0']],
             epsR = rep(0,nyrec),
             lsigmaR = lsigmaR,
             ilogeta = res4$par[['ilogeta']],
             ilogpim =unname(res4$par[grep('ilogpim',names(res4$par))]),
             selparsATT = unname(res4$par[grep('selparsATT',names(res4$par))]),
             selparsNVT = unname(res4$par[grep('selparsNVT',names(res4$par))]), 
             selparsATL = unname(res4$par[grep('selparsATL',names(res4$par))]), 
             selparsNMRL = unname(res4$par[grep('selparsNMRL',names(res4$par))]), 
             selparsSMRL = unname(res4$par[grep('selparsSMRL',names(res4$par))]), 
             xphiLF = xphiLF,
             xphiT = xphiT 
             )

obj5 <- MakeADFun(data=data,parameters=pars5,map=map5,DLL="v2_dyndep")
res5 <- do.call("optim",obj5)
rep5 <- obj5$rep()
sdrep5 <- sdreport(obj5)
summary(sdrep5)

# phase 6: over-dispersion parameters

map6 <- list(lsigmaR = factor(NA))

pars6 <- list(logR0 = res5$par[['logR0']],
             epsR = unname(res5$par[grep('epsR',names(res5$par))]),
             lsigmaR = lsigmaR,
             ilogeta = res5$par[['ilogeta']],
             ilogpim =unname(res5$par[grep('ilogpim',names(res5$par))]),
             selparsATT = unname(res5$par[grep('selparsATT',names(res5$par))]),
             selparsNVT = unname(res5$par[grep('selparsNVT',names(res5$par))]), 
             selparsATL = unname(res5$par[grep('selparsATL',names(res5$par))]), 
             selparsNMRL = unname(res5$par[grep('selparsNMRL',names(res5$par))]), 
             selparsSMRL = unname(res5$par[grep('selparsSMRL',names(res5$par))]), 
             xphiLF = xphiLF,
             xphiT = xphiT 
             )

obj6 <- MakeADFun(data=data,parameters=pars6,map=map6,DLL="v2_dyndep")
res6 <- do.call("optim",obj6)
rep6 <- obj6$rep()
sdrep6 <- sdreport(obj6)
summary(sdrep6)

# phase 7: recruitment SD (random effects bit)

pars7 <- list(logR0 = res6$par[['logR0']],
             epsR = unname(res6$par[grep('epsR',names(res5$par))]),
             lsigmaR = lsigmaR,
             ilogeta = res6$par[['ilogeta']],
             ilogpim =unname(res6$par[grep('ilogpim',names(res6$par))]),
             selparsATT = unname(res6$par[grep('selparsATT',names(res6$par))]),
             selparsNVT = unname(res6$par[grep('selparsNVT',names(res6$par))]), 
             selparsATL = unname(res6$par[grep('selparsATL',names(res6$par))]), 
             selparsNMRL = unname(res6$par[grep('selparsNMRL',names(res6$par))]), 
             selparsSMRL = unname(res6$par[grep('selparsSMRL',names(res6$par))]), 
             xphiLF = unname(res6$par[grep('xphiLF',names(res6$par))]),
             xphiT = res6$par[['xphiT']] 
             )

obj7 <- MakeADFun(data=data,parameters=pars7,map=map6,DLL="v2_dyndep",random="epsR")
res7 <- do.call("optim",obj7)
rep7 <- obj7$rep()
sdrep7 <- sdreport(obj7)
summary(sdrep6)



# plots

phiLF <- rep6$phiLF
names(phiLF) <- colnames(C)
repp <- rep6
ff <- 'ATT'
plot.lenfreq(LF,LFcov,repp,ff,phiLF[ff]) 
plot.agelendat(ALF,ALFcov,repp,'SMRL','female')
plot.tagfits.relyr(repp,tagrel,tagnrec,tagcov,tagrec,nrecmax)
plot.tagfits.relyr.tot(repp,tagrel,tagnrec,tagcov,tagrec,nrecmax)
plot.tagfits.recyr.tot(repp,tagrel,tagnrec,tagcov,tagrec,nrecmax)
plot.tagfits.recyrarea(repp,tagrel,tagnrec,tagcov,tagrec,nrecmax)

######################
# MCMC using tmbstan #
######################a

# for now fix eta at MLE (MCMC issues...)

map7 <- list(ilogeta = factor(NA),
             logl1 = rep(factor(NA),ns),
             logl2 = rep(factor(NA),ns),
             logk = rep(factor(NA),ns),
             logsdla = rep(factor(NA),ns))

pars7 <- list(logR0 = res5rw$par[['logR0']],
             epsR = rep(0,nyrec),
             ilogeta = res5rw$par[['ilogeta']],
             ilogpim =unname(res5rw$par[grep('ilogpim',names(res5rw$par))]),
             selparsATT = unname(res5rw$par[grep('selparsATT',names(res5rw$par))]),
             selparsNVT = unname(res5rw$par[grep('selparsNVT',names(res5rw$par))]), 
             selparsATL = unname(res5rw$par[grep('selparsATL',names(res5rw$par))]), 
             selparsNMRL = unname(res5rw$par[grep('selparsNMRL',names(res5rw$par))]), 
             selparsSMRL = unname(res5rw$par[grep('selparsSMRL',names(res5rw$par))]), 
             logl1 = log(l1),
             logl2 = log(l2),
             logk = log(k),
             logsdla = log(sdla))

obj7 <- MakeADFun(data=data,parameters=pars7,map=map7,DLL="v2_dyndep")
res7 <- do.call("optim",obj7)
rep7 <- obj7$rep()
sdrep7 <- sdreport(obj7)
summary(sdrep7)

library(tmbstan)
options(mc.cores = parallel::detectCores())

mcmc <- tmbstan(obj7, init = 'last.par.best', iter = 1200, warmup = 200, chains=1,control=list(max_treedepth=12))

traceplot(mcmc, pars=names(obj7$par), inc_warmup=TRUE)

traceplot(mcmc, pars=names(obj7$par), inc_warmup=FALSE)

pairs(mcmc, pars=names(obj7$par)[1])

parlist <- extract(mcmc)


# get a sample of ilogeta using multi-variate normal approximation
# does for parameters 30:npar

library(mvtnorm)
npar <- length(res5rw$par)
Sigmax <- sdrep5rw$cov.fixed[30:npar,30:npar]
mu1 <- res5rw$par[['ilogeta']]
mu2 <- unname(res5rw$par[31:npar])
Sigma22 <- Sigmax[-1,]
Sigma22 <- Sigma22[,-1]
Sigma12 <- Sigmax[1,-1]
Sigma21 <- Sigmax[-1,1]
Sigmabar <- Sigmax[1,1]-t(Sigma12) %*% solve(Sigma22) %*% Sigma21

# create named matrix of MCMC iterations

parmat <- matrix(nrow=1000,ncol=npar-1)
colnames(parmat) <- names(res5rw$par[-30])
for(i in 1:length(parlist)) {

  xnm <- names(parlist)[i]
  parmat[,grep(xnm,colnames(parmat))] <- parlist[[i]]

}

etamc <- rep(NA,1000)
for(i in 1:1000) {

  a <- parmat[i,-c(1:29)]
  mubar <- mu1+t(Sigma12) %*% solve(Sigma22) %*% (a-mu2)
  etamc[i] <- rnorm(1,mubar,Sigmabar)

}

parmat2 <- cbind(parmat[,1:29],etamc,parmat[,-c(1:29)])
colnames(parmat2)[[30]] <- 'ilogeta'

library(coda)
gdiag <- geweke.diag(parmat2,frac1 = 0.45, frac2 = 0.5)
plot(gdiag[['z']],xlab='Parameter',ylab='Z-value',pch=19,col='blue',cex=0.75)
abline(h=1.96,lty=2,lwd=1.5,col='magenta')
abline(h=-1.96,lty=2,lwd=1.5,col='magenta')

save.image("run1_9tagev_newmat_dyndep.rda")

