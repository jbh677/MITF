# TMB-based conditional age-at-length estimator
# schnute LVB parameterisation

require(TMB)

# compile it; load it

compile("condAatL_schnute.cpp")
dyn.load(dynlib("condAatL_schnute"))

# unload in debug mode

#dyn.unload(dynlib("condAatL"))

# get the data

g.df <- read.csv("../../../data/toothfish2021/Data/2021/20210212/CSIRO-GF-AgeDataFirstReading.csv",header=T)
names(g.df)[[3]] <- 'Length'
names(g.df)[[5]] <- 'Age'
names(g.df)[[6]] <- 'Year'
g.df <- subset(g.df,Sex == 'M' | Sex == 'F')
g.df$Length <- g.df$Length/10

library(ggplot2)
ggplot(g.df,aes(x=Age,y=Length))+geom_point(colour='blue',size=1)+facet_wrap(~Sex)+xlab("age")+ylab("length")

# which sex?

sex <- 'F'
g2.df <- subset(g.df,Sex==sex)

# get rid of NAs

g2.df <- subset(g2.df,!is.na(Length))
g2.df <- subset(g2.df,!is.na(Age))
g2.df <- subset(g2.df,!is.na(Year))

# plots

ggplot(g2.df,aes(x=Age,y=Length))+geom_point(colour='blue',size=1)+xlab("age")+ylab("length")

ggplot(g2.df,aes(Length))+geom_histogram(colour='blue')+facet_wrap(~Year)

# get it into format for estimation

lrng <- range(g2.df$Length)
arng <- range(g2.df$Age)
yrng <- range(g2.df$Year)
yrs <- sort(unique(g2.df$Year))
nysam <- length(yrs)
ages <- 1:arng[2]
nages <- length(ages)
lrng
lbins <- seq(20,200,by=10)
lref <- 0.5*(lbins[-c(1)]+lbins[-c(length(lbins))])  
nbins <- length(lref)

lf <- matrix(0,ncol=nbins,nrow=nysam)
lf.df <- expand.grid(Length=lref,Year=yrs,n=NA)
for(y in 1:nysam) {
  lsam <- subset(g2.df,Year==yrs[y])$Length
  for(l in 1:nbins) {
      lu <- lbins[l+1]
      ll <- lbins[l]
      lf[y,l] <- length(lsam[lsam>=ll & lsam<lu]) 
  }
  lf.df[(1+(y-1)*nbins):(nbins+(y-1)*nbins),'n'] <- lf[y,]
}
apply(lf,1,sum)

ggplot(lf.df,aes(x=Length,y=n))+geom_line(colour='blue')+facet_wrap(~Year)

alk <- array(0,dim=c(nysam,nbins,nages))

for(y in 1:nysam) {
    
    tmp1.df <- subset(g2.df,Year==yrs[y])
    tmp2.df <- data.frame(l=tmp1.df$Length,a=tmp1.df$Age)
    for(l in 1:nbins) {

      lu <- lbins[l+1]
      ll <- lbins[l]
      df <- subset(tmp2.df,l>=ll & l<lu)
      if(dim(df)[1] > 0) {
         
         tb <- apply(table(df),2,sum)
         arng <- as.numeric(names(tb))
         anm <- names(tb)
         for(a in 1:length(arng)) {

           alk[y,l,arng[a]] <- tb[[anm[a]]]

         }
      }
     }
}

apply(alk,1,sum)

# ageing error matrix (diagnonal => perfect)

#A <- diag(1,length(ages))

aerr.df <- read.table(file="age_err.dat",header=T)
A <- diag(NA,length(ages))
a.df <- subset(aerr.df,True <= ages[length(ages)])
for(j in 1:length(ages)) {

  atrue <- ages[j]
  a.sd <- subset(a.df,True==atrue)$SD
  pa <- dnorm(ages,atrue,a.sd)
  pa <- pa/sum(pa)
  A[,j] <- pa
}

apply(A,2,sum)

# pick reference ages for Schnute application

a1 <- 5
a2 <- 20

# priors

mupr <- log(mean(ages))
cvpr <- 1
sdpr <- sqrt(log(1+cvpr^2))
plnorm(max(ages),log(mupr),sdpr)
musigma <- 0
sdsigma <- 0.25

# data object

data <- list(lf=lf,alk=alk,lref=lref,aref=ages,alref=c(a1,a2),dms=dim(alk),A=A,priors=c(mupr,sdpr,musigma,sdsigma))

# initial estimates

linf <- 165
k <- 0.055
t0 <- -0.78
cvl <- 0.13
sdl <- sqrt(log(1+cvl^2))
l1 <- linf*(1-exp(-k*(a1-t0)))
l2 <- linf*(1-exp(-k*(a2-t0)))

pl <- t(apply(lf,1,function(x){x/sum(x)}))
mul <- apply(pl,1,function(x,lref){x <- sum(x*lref)},lref)
mua <- t0-1/k*log(1-mul/linf)
cva <- rep(0.25,nysam)

# parameter list

pars <- list(
                lnk = log(k),
                logl1 = log(l1),
                logl2 = log(l2),
                lncvl = log(cvl),
                lnmua = log(mua),
                lnsda = log(sqrt(log(1+cva^2)))
              )

# create AD object

obj <- MakeADFun(data=data,parameters=pars,DLL="condAatL_schnute")

# run it

res <- do.call("optim",obj)

exp(res$par)

sdrep <- sdreport(obj)
zz <- summary(sdrep,"fixed")

rep <- obj$report()

# tmbstan MCMC

library(tmbstan)

mcmc <- tmbstan(obj, init ='last.par.best' , iter = 2000, chains=1,control=list(max_treedepth=12,adapt_delta=0.8))

traceplot(mcmc, pars=names(obj$par)[1:4], inc_warmup=TRUE)

traceplot(mcmc, pars=names(obj$par)[1:4], inc_warmup=FALSE)

pairs(mcmc, pars=names(obj$par)[1:4])

parlist <- extract(mcmc)

###################
# fitting summary #
###################

# LF data

lf.hat <- rep$pl
ntot <- apply(lf,1,sum)
lf.df$nhat <- rep(NA,prod(nysam*nbins))

for(y in 1:nysam) lf.df[(1+(y-1)*nbins):(nbins+(y-1)*nbins),'nhat'] <- lf.hat[y,]*ntot[y]

ggplot(lf.df,aes(x=Length,y=n))+geom_point(colour='magenta')+geom_line(aes(x=Length,y=nhat),colour='blue')+facet_wrap(~Year)

# A | L data

pal.hat <- rep$pal

alk.df <- expand.grid(Length=lref,Year=yrs,mua=NA,mua.hat=NA,lq=NA,uq=NA)
for(y in 1:nysam) {
  for(l in 1:nbins) {

    aobs <- alk[y,l,]
    if(sum(aobs)>0) pobs <- aobs/sum(aobs)
    phat <- pal.hat[y,l,]
    ahat <- phat*lf[y,l]
    mua <- ifelse(sum(aobs) > 0,sum(ages*pobs),NA)
    mua.hat <- sum(ages*phat)
    sd.hat <- sqrt(sum(phat*(ages-mua.hat)^2))
    lq <- max(0,mua.hat-1.96*sd.hat)
    uq <- mua.hat+1.96*sd.hat
    alk.df[l+(y-1)*nbins,'mua'] <- mua
    alk.df[l+(y-1)*nbins,'mua.hat'] <- mua.hat
    alk.df[l+(y-1)*nbins,'lq'] <- lq
    alk.df[l+(y-1)*nbins,'uq'] <- uq 
  }
}

ggplot(alk.df,aes(x=Length,y=mua))+geom_point(colour='magenta')+geom_line(aes(x=Length,y=mua.hat))+geom_line(aes(x=Length,y=uq),linetype='dashed',colour='blue')+geom_line(aes(x=Length,y=lq),linetype='dashed',colour='blue')+facet_wrap(~Year)

# plot p(a)

pa.hat <- rep$pa
rownames(pa.hat) <- yrs
colnames(pa.hat) <- ages
pa.df <- expand.grid(age=ages,year=yrs,p=NA)
for(y in 1:nysam) pa.df[pa.df$year == yrs[y],'p'] <- pa.hat[y,]

ggplot(pa.df,aes(x=age,y=p))+geom_line(colour='blue')+facet_wrap(~year)

# plot p(l | a)

pla.hat <- rep$pla
rownames(pla.hat) <- lref
colnames(pla.hat) <- ages
levelplot(t(pla.hat),xlab='Age',ylab='Length',main='P(length | age)',scales=list(y=list(rot=45), x=list(rot=45)),cex=0.5)

# plot p(a | l) over years using levelplot

pal.df <- expand.grid(age=ages,year=yrs,length=lref,p=NA)
for(l in 1:nbins)
  for(y in 1:nysam) 
    pal.df[pal.df$length == lref[l] & pal.df$year == yrs[y],'p'] <- pal.hat[y,l,]

levelplot(p~length*age|as.factor(year),pal.df,main='P(age | length)',cex=0.5)

fname <- paste(paste(sex,"schnute",sep="_"),"rda",sep=".")
save.image(fname)


