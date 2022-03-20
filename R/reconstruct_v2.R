###########################################################
# reconstruction of dynamics from MCMC samples ############
###########################################################
# R. Hillary CSIRO O & A 2019 #############################
###########################################################

library(TMB)
load("run1_9tagev_newmat.rda")
source("v2_utils.R")
nits <- dim(parmat)[1]

# "active" parameters

parnm1 <- colnames(parmat)

# "fixed" parameters

parnm2 <- names(map8)

summres <- recon2(nits,parmat,pars8,parnm1,parnm2,'new') 
save.image("recon_9tagev_newmat.rda")

#################
# summary plots #
#################

library(ggplot2)

# female SSB

sq <- apply(summres$SSBf,2,quantile,c(0.025,0.5,0.975))
df <- data.frame(y=yrs,med=sq[2,],lq=sq[1,],uq=sq[3,])
pp1 <- ggplot()+geom_errorbar(data=df,mapping=aes(x=ISOdate(y,1,1),ymin=lq,ymax=uq),width=0,col='blue')
pp1 <- pp1+geom_point(data=df, mapping=aes(x=ISOdate(y,1,1), y=med), size=1, shape=21, colour="blue")
pp1 <- pp1+xlab("year")+ylab("female SSB")+ylim(c(0,max(sq[3,])))
pp1

# female SSB depletion

dsq <- apply(summres$SSBf/summres$SSBf[,1],2,quantile,c(0.025,0.5,0.975))
df2 <- data.frame(y=yrs,med=dsq[2,],lq=dsq[1,],uq=dsq[3,])
pp2 <- ggplot()+geom_hline(yintercept=1,lty=2)+geom_errorbar(data=df2,mapping=aes(x=ISOdate(y,1,1),ymin=lq,ymax=uq),width=0,col='blue')
pp2 <- pp2+geom_point(data=df2, mapping=aes(x=ISOdate(y,1,1), y=med), size=1, shape=21, colour="blue")
pp2 <- pp2+xlab("year")+ylab("female SSB depletion")+ylim(c(0,1.15))
pp2

# total recruitment

rq <- apply(summres$Rec,2,quantile,c(0.025,0.5,0.975))
df3 <- data.frame(y=yrs,med=rq[2,],lq=rq[1,],uq=rq[3,])
pp3 <- ggplot()+geom_errorbar(data=df3,mapping=aes(x=ISOdate(y,1,1),ymin=lq,ymax=uq),width=0,col='blue')
pp3 <- pp3+geom_point(data=df3, mapping=aes(x=ISOdate(y,1,1), y=med), size=1, shape=21, colour="blue")
pp3 <- pp3+xlab("year")+ylab("Recruitment")+ylim(c(0,max(rq[3,])))
pp3

# area specific female SSB

srq <- apply(summres$SSB[,,1,],c(2,3),quantile,c(0.025,0.5,0.975))
df4 <- expand.grid(y=yrs,r=c("North","South"),med=NA,lq=NA,uq=NA)
df4[df4$r == 'North','med'] <- srq[2,,1]
df4[df4$r == 'North','lq'] <- srq[1,,1]
df4[df4$r == 'North','uq'] <- srq[3,,1]
df4[df4$r == 'South','med'] <- srq[2,,2]
df4[df4$r == 'South','lq'] <- srq[1,,2]
df4[df4$r == 'South','uq'] <- srq[3,,2]
pp4 <- ggplot()+geom_errorbar(data=df4,mapping=aes(x=ISOdate(y,1,1),ymin=lq,ymax=uq),width=0,col='blue')
pp4 <- pp4+geom_point(data=df4, mapping=aes(x=ISOdate(y,1,1), y=med), size=1, shape=21, colour="blue")
pp4 <- pp4+facet_wrap(~r)+xlab("year")+ylab("female SSB")+ylim(c(0,max(srq[3,,])))
pp4

# area specific female SSB depletion

dsr <- summres$SSB[,,1,]
dsr[,,1] <- dsr[,,1]/dsr[,1,1]
dsr[,,2] <- dsr[,,2]/dsr[,1,2]
dsrq <- apply(dsr,c(2,3),quantile,c(0.025,0.5,0.975))
df5 <- expand.grid(y=yrs,r=c("North","South"),med=NA,lq=NA,uq=NA)
df5[df5$r == 'North','med'] <- dsrq[2,,1]
df5[df5$r == 'North','lq'] <- dsrq[1,,1]
df5[df5$r == 'North','uq'] <- dsrq[3,,1]
df5[df5$r == 'South','med'] <- dsrq[2,,2]
df5[df5$r == 'South','lq'] <- dsrq[1,,2]
df5[df5$r == 'South','uq'] <- dsrq[3,,2]
pp5 <- ggplot()+geom_errorbar(data=df5,mapping=aes(x=ISOdate(y,1,1),ymin=lq,ymax=uq),width=0,col='blue')
pp5 <- pp5+geom_point(data=df5, mapping=aes(x=ISOdate(y,1,1), y=med), size=1, shape=21, colour="blue")
pp5 <- pp5+facet_wrap(~r)+xlab("year")+ylab("female SSB depletion")+ylim(c(0,1.15))+geom_hline(yintercept=1,lty=2)
pp5

# spatial parameters

xdf <- expand.grid(var=c('eta','N-to-S','S-to-N'),val=NA,iter=1:nits)
xdf[xdf$var == 'eta','val'] <- ilogit(parmat[,33])
xdf[xdf$var == 'N-to-S','val'] <- 1-ilogit(parmat[,34])
xdf[xdf$var == 'S-to-N','val'] <- 1-ilogit(parmat[,35])
pp6 <- ggplot(xdf,aes(val))+geom_histogram(colour='blue')+facet_wrap(~var,scales = 'free_x')+xlab("Parameter value")+ylab("Posterior")
pp6

# spatial harvest rates (final year)

hfin <- apply(summres$hyafin,c(1,3,4),max)

################
# status stats #
################

# total SSB deleption

round(df2[dim(df2)[1],],2)

# spatial SSB depletion

subset(df5,y=='2020')

# recruitment fraction and migration

round(quantile(subset(xdf,var == 'eta')$val,c(0.025,0.5,0.975)),2)
round(quantile(subset(xdf,var == 'N-to-S')$val,c(0.025,0.5,0.975)),3)
round(quantile(subset(xdf,var == 'S-to-N')$val,c(0.025,0.5,0.975)),3)

###############
# projections #
###############

# check the most recent spatial catch distros

Cnow <- C[(ny-4):ny,]
t(apply(Cnow,1,function(x){x <- x/sum(x)}))

tac <- 638
tacatl <- 300
nssplit <- c(0.75,0.25)
tac.split <- get.tacsplit(tac,tacatl,nssplit,ffnm)
rectype <- 'nonparametric'
mat.type <- 'new'
set.seed(23)
prj <- proj2(nits,35,parmat,pars8,parnm1,parnm2,summres$SSBf[,1],summres$Nfin,summres$hyafin,unname(tau[ny,]),tac.split,rectype,mat.type)

# join historical with projected

ssbtot <- cbind(summres$SSBf,prj$SSBf)
dx <- ssbtot/ssbtot[,1]
pdx1 <- apply(dx,2,function(x){length(x[x>=0.5])/length(x)})
pdx1
round(tac.split)
pdx2 <- apply(dx,2,function(x){length(x[x<0.2])/length(x)})
pdx2
colnames(dx) <- c(yrs,c((yrs[length(yrs)]+1):(yrs[length(yrs)]+35)))
boxplot(dx,outline=FALSE,col='green',ylim=c(0,1.1),xlab='year',ylab='Female SSB depletion')
abline(h=0.5,lty=2)
abline(v=yrs[length(yrs)]-yrs[1]+1,lty=2)

save.image("proj_9tagev_newmat.rda")

