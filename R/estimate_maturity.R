######################################################
# estimate maturity ##################################
######################################################
# R. Hillary CSIRO 2021 O & A ########################
######################################################


# get the data

mdf <- read.csv("../../data/toothfish2021/Data/2021/20210212/CSIRO-GF-FishMaturity.csv",header=T)

library(ggplot2)

ggplot(mdf,aes(x=TL..mm.,y=Gonad.State))+geom_point(colour='blue',size=0.25)+facet_wrap(~Sex)

ggplot(subset(mdf,Sex == 'M' | Sex == 'F'),aes(x=TL..mm./1e+3,y=Gonad.State))+geom_point(colour='blue',size=0.25)+facet_wrap(~Sex)+xlab("Length (m)")+ylab("Maturity stage")

mdf <- subset(mdf,!is.na(TL..mm.))
mdf <- subset(mdf,!is.na(Gonad.State))

ggplot(subset(mdf,Sex == 'M' | Sex == 'F'),aes(TL..mm./1e+3))+geom_histogram(colour='blue')+facet_grid(Gonad.State~Sex,scales='free_y')+xlab("Length (m)")+ylab("Number of samples")

# look for sus ones

weird <- subset(mdf,Gonad.State >= 4 & TL..mm. <= 700) # 

# get the length data

awl <- 4.4e-6
bwl <- 3.14
ldf <- read.csv("../../data/toothfish2021/Data/2021/20210212/CSIRO-GF-FishLengths.csv",header=T)
lx <- match(weird$Fish.Serial.Number,ldf$Fish.Serial.Number)
ldf$TL..mm.[lx]-weird$TL..mm. # length right?
ldf$Weight..g.[lx] 
awl*(ldf$TL..mm.[lx]/10)^bwl*1e+3
(1-ldf$Weight..g.[lx]/(awl*(ldf$TL..mm.[lx]/10)^bwl*1e+3))*100

###########################################################################
# generate the actual data (number of fish stage 2+) for given length bin #
###########################################################################

# choose the sex

sex <- 'F'
mdf$Sex <- as.character(mdf$Sex)
mdf <- subset(mdf,Sex==sex)
mdf$TL..cm. <- mdf$TL..mm./10

# length bins

gdf <- data.frame(gstate=mdf$Gonad.State,l=mdf$TL..cm.)
tplus <- 2 # number of years to add growth increment
grdf <- subset(read.table("growth.tab",header=T),Sex==sex)
ginc <- (grdf$Linf-gdf$l)*(1-exp(-grdf$k*tplus))

range(gdf$l)
ll <- 30 # cm not mm
lu <- 200 # cm not mm
del <- 5
lbins <- c(30,40,50,60,seq(65,140,by=5),150,160,200) # female
#lbins <- ll+(0:((lu-ll)/del))*del
mulbins <- 0.5*(lbins[-1]+lbins[-length(lbins)])
nbins <- length(lbins)-1
gcrit <- 2 # critical maturity state
gX <- expand.grid(l=mulbins,n=NA,k=NA)
gXX <- expand.grid(l=mulbins,n=NA,k2=NA,k3p=NA)
for(i in 1:nbins) {

  zdf <- subset(gdf,l >= lbins[i] & l < lbins[i+1])
  d1 <- dim(zdf)[1]
  d2 <- dim(subset(zdf,gstate >= gcrit))[1]
  d2only <- dim(subset(zdf,gstate == gcrit))[1]
  d3plus <- dim(subset(zdf,gstate > gcrit))[1]
  gX[i,'k'] <- d2
  gX[i,'n'] <- d1
  gXX[i,'n'] <- d1
  gXX[i,'k2'] <- d2only
  gXX[i,'k3p'] <- d3plus

}
gXX <- subset(gXX,k2+k3p>0)

px <- gX$k/gX$n
plot(mulbins,px,xlab='length',ylab='P[mstate = 2+]',type='l',col='blue',lwd=1.5)

phimax <- min(gXX$n)-0.05
phimin <- 1.05

# beta-binomial likelihood for proportion mature by size-class

matnloglbb <- function(theta) {

  #phi <- 1.05+exp(theta[1]) # OD (minimum 1.05)
  nx <- exp(theta[1]) # power parameter
  mux <- exp(theta[2]) # half-saturation (effectively lm50)
  phi <- (phimax*exp(theta[3])+phimin)/(1+exp(theta[3]))

  nlogl <- 0
  for(i in 1:dim(gXX)[1]) {
  
    l <- gXX$l[i]
    k2 <- gXX$k2[i]
    k3 <- gXX$k3p[i]
    k <- k2+k3
    n <- gXX$n[i]
    pk <- c(k2,k3)/sum(c(k2,k3))
    l2 <- l+(grdf$Linf-l)*(1-exp(-grdf$k*tplus))
    l3 <- l
    lref <- sum(pk*c(l2,l3))
    pmu <- lref^nx/(mux^nx+lref^nx)
    alpx <- ((n-phi)*pmu)/((1-pmu)*(pmu+(1-pmu)*(phi-1)))
    betx <- (n-phi)/((1-pmu)*(pmu+(1-pmu)*(phi-1)))
    nlogl <- nlogl+lgamma(k+alpx)+lgamma(n-k+betx)+lgamma(alpx+betx)-lgamma(n+alpx+betx)-lgamma(alpx)-lgamma(betx)

  }

  return(-nlogl)
}

matnloglbin <- function(theta) {

  nx <- exp(theta[1]) # power parameter
  mux <- exp(theta[2]) # half-saturation (effectively lm50)

  nlogl <- 0
  for(i in 1:dim(gXX)[1]) {
  
    l <- gXX$l[i]
    k2 <- gXX$k2[i]
    k3 <- gXX$k3p[i]
    if(k2 > 0 | k3 > 0) {
      
      n <- gXX$n[i]
      k <- k2+k3 
      pk <- c(k2,k3)/sum(c(k2,k3))
      l2 <- l+(grdf$Linf-l)*(1-exp(-grdf$k*tplus))
      l3 <- l
      lref <- sum(pk*c(l2,l3))
      pmu <- lref^nx/(mux^nx+lref^nx) 
      nlogl <- nlogl+dbinom(k,n,pmu,log=TRUE)   
      #cat(nlogl,"\n")
    }
  }

  return(-nlogl)
}

# binomial

theta1 <- log(c(6,97))
res <- optim(theta1,matnloglbin,method="L-BFGS-B",hessian=T)
# estimates
exp(res$par)
# CV
sqrt(diag(solve(res$hessian)))/exp(res$par)/exp(res$par)

# beta-binomial

theta2 <- log(c(6,97,(2-phimin)/(phimax-2)))
resbb <- optim(theta2,matnloglbb,method="L-BFGS-B",hessian=T,control=list(trace=1))
# estimates
exp(resbb$par)
# CV
sqrt(diag(solve(resbb$hessian)))/exp(resbb$par)/exp(resbb$par)

# get 50% and 90% maturity

ntmp <- exp(res$par[1]) 
l50 <- mutmp <- exp(res$par[2])
psi <- 0.95
l95 <- (psi*mutmp^ntmp/(1-psi))^{1/ntmp}

# Predicted vs. observed

nx <- exp(res$par[1])
mux <- exp(res$par[2])
khat <- gXX$k2
for(i in 1:length(khat)) {

    l <- gXX$l[i]
    k2 <- gXX$k2[i]
    k3 <- gXX$k3p[i]
    n <- gXX$n[i]
    k <- k2+k3 
    pk <- c(k2,k3)/sum(c(k2,k3))
    l2 <- l+(grdf$Linf-l)*(1-exp(-grdf$k*tplus))
    l3 <- l
    lref <- sum(pk*c(l2,l3))
    pmu <- lref^nx/(mux^nx+lref^nx) 
    khat[i] <- n*lref^nx/(mux^nx+lref^nx) 
}

ymax <- max(c(max(khat),max(gXX$k2+gXX$k3p)))
plot(gXX$l,khat,type='l',xlab='length',ylab='No. > gcrit',col='blue',ylim=c(0,ymax))
points(gXX$l,gXX$k2+gXX$k3p,pch=19,col='magenta')

# standardised residuals

phat <- khat/gXX$n
kres <- (gXX$k2+gXX$k3p-khat)/(khat*(1-phat)) 
sd(kres)

# actual maturity-at-length

lx <- seq(0,lu,length=100)
matl <- lx^nx/(mux^nx+lx^nx)
plot(lx,matl,xlab='length',ylab='maturity',type='l',col='blue',lwd=1.5)
abline(h=0.5,lty=2)
abline(v=mux,lty=2)

save.image(paste(sex,"_matl_2021.rda",sep=""))

# both

load("F_matl_2021.rda")
lz <- seq(0,210,length=100)
matlf <- lz^nx/(mux^nx+lz^nx)
lf50 <- mux
load("M_matl_2021.rda")
lm50 <- mux
matlm <- lz^nx/(mux^nx+lz^nx)
plot(lz,matlf,xlab='length',ylab='maturity',type='l',col='blue',lwd=1.5)
lines(lz,matlm,lty=1,lwd=1.5,col='black')
abline(h=0.5,lty=2,col='red')
abline(v=lf50,lty=2,col='blue')
abline(v=lm50,lty=2,col='black')
legend("topleft",lty=c(1,1,2),lwd=c(1.5,1.5,1.5),col=c("blue","black","red"),legend=c("Female","Male","50% maturity"),bty='n')

########################################################################
# comparing female maturity-at-age for current and alternative options #
########################################################################

lm50 <- 1.396
lm95 <- 1.855
lbins <- seq(0.1,2.2,by=0.1)

load("gpars.rda")
fg <- subset(gpars,sex=='F')
k <- fg[fg$par=='k','val']
linf <- fg[fg$par=='linf','val']
t0 <- fg[fg$par=='t0','val']
cvl <- fg[fg$par=='cvl','val']
sdla <- sqrt(log(1+cvl^2))
ages <- 1:50
nages <- length(ages)
mula <- linf*(1-exp(-k*(ages-t0)))
mata.old <- mata.new <- rep(NA,length(ages))
for(a in 1:nages) {

  ll <- max(1e-3,mula[a]-2*mula[a]*cvl)
  lu <- mula[a]+2*mula[a]*cvl
  lref <- seq(ll,lu,length=50)
  dl <- dlnorm(lref,log(mula[a]),sdla,F)
  dl <- dl/sum(dl)
  mold <- 1/(1+19^{-(lref-lm50)/(lm95-lm50)})
  mnew <- lref^nx/(mux^nx+lref^nx)
  mata.old[a] <- sum(dl*mold)
  mata.new[a] <- sum(dl*mnew)

}

plot(ages,mata.old,xlab='age',ylab='Expected maturity',ylim=c(0,1),main='',type='l',lwd=1.5)
lines(ages,mata.new,lty=1,col='blue',lwd=1.5)
legend('topleft',lty=c(1,1),lwd=c(1.5,1.5),col=c('black','blue'),legend=c('Old','Revised'),bty='n')

