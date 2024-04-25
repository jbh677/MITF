##################################################################
# script that runs MSE loop (function-of-a-function nonsense...) #
##################################################################
# R. Hillary & P. Bessell-Browne CSIRO (2024) ####################
##################################################################

# initialisation of popn via Biol & Fleet

cat("Population initialisation...","\n")
prj <- proj.para(Biol,Fleet,msectrl$nyinit,FALSE,nitsx,ncore,TRUE) 
ycurr <- msectrl$nyinit

# generate initial tag data 

cat("Generate initial tag data...","\n")
prjlist <- listify.proj(prj)
tdflist <- mclapply(1:nits,wrap.get.tag.data2,mc.cores=ncore)
rm(prjlist)

# get initial estimates of harvest rates

cat("Get initial harvest rate estimates...","\n")
model <- msectrl$initmod
mctrl <- msectrl$mctrl
replist <- mclapply(1:nits,wrap.get.hrates,mc.cores=ncore)

if(model == 'tgm1') {
  
  hyhat <- array(dim=c(length(replist[[1]]$hy),nits))
  for(nn in 1:nits) hyhat[,nn] <- replist[[nn]]$hy

} 

if(model == 'tgm1re') {
  
  hyhat <- array(dim=c(length(replist[[1]]$hy),nits))
  sigmahat <- rep(NA,nits)
  for(nn in 1:nits) {

    hyhat[,nn] <- replist[[nn]]$hy
    sigmahat[nn] <- replist[[nn]]$sigmah

  }
} 

if(model == 'tmg2') {

  Phi.mc <- array(dim=c(nits,nr,nr)) 
  hyhat <- array(dim=c(dim(replist[[1]]$hy,nits)))
  for(nn in 1:nits) {

    Phi.mc[nn,,] <- replist[[nn]]$Phi 
    hyhat[,,nn] <- replist[[nn]]$hy 

  }
}

if(model == 'tmg2re') {

  Phi.mc <- array(dim=c(nits,nr,nr)) 
  hyhat <- array(dim=c(dim(replist[[1]]$hy,nits)))
  sigmahat <- rep(NA,nits) 
  for(nn in 1:nits) {

    Phi.mc[nn,,] <- replist[[nn]]$Phi 
    hyhat[,,nn] <- replist[[nn]]$hy 
    sigmahat[nn] <- replist[[nn]]$sigmah 

  }
} 
rm(replist)

##########################
# start of main MSE loop #
##########################

cat("Start of main MSE loop...","\n")
CatchControl <- TRUE

# calculate number and length of projection periods in between TAC decisions

ycurr <- msectrl$firstdec
ndec <- ceiling((max(msectrl$yrng)-msectrl$firstdec)/msectrl$tactau)
ydec <- ycurr+(0:(ndec-1))*msectrl$tactau

#################################
# begin loop over TAC decisions #
#################################

# first decision

yprj <- ydec[1]+msectrl$tactau 
Cterm <- prj$C[ycurr,,]
TACold <- apply(Cterm,2,sum)

# apply the MP for the "first time"

cat("First TAC decision...","\n")
TACnew <- get.TAC(tacctrl)

# prep Biol and Fleet objects

Biolx <- Biol
Fleetx <- Fleet
Biolx$Ninit <- prj$N[ycurr,,,,]
Fleetx$Cfix <- tacctrl$catchSplit %o% TACnew

# project forward

prjtmp <- get.popdyn(Biolx,Fleetx,msectrl$tactau,CatchControl,nits) 

# joining it all up

prj$N <- abind(prj$N,prjtmp$N,along=1)
prj$H <- abind(prj$H,prjtmp$H,along=1) 
prj$C <- abind(prj$C,prjtmp$C,along=1) 
prj$SSB <- abind(prj$SSB,prjtmp$SSB,along=1) 
prj$SSBtot <- abind(prj$SSBtot,prjtmp$SSBtot,along=1) 
rm(prjtmp) 

for(j in 2:ndec) { 
 
  # extend existing tag data into projection period

  prjlist <- listify.proj(prj)
  txdflist <- mclapply(1:nits,wrap.extend.tag.data,mc.cores=ncore) 

  # generate new tag data in the next projection period

  Tag$ytagon <- ycurr
  Tag$ytagoff <- yprj-1
  tnewdf <- mclapply(1:nits,wrap.get.tag.data2,mc.cores=ncore)

  # merge the continuation and newly released tag data frames

  tdflist <- list()
  for(nn in 1:nits) tdflist[[nn]] <- rbind(txdflist[[nn]],tnewdf[[nn]])
  rm(txdflist,tnewdf,prjlist)

  # re-estimate harvest rates

  model <- msectrl$prjmod
  mctrl <- msectrl$mctrl
  replist <- mclapply(1:nits,wrap.get.hrates,mc.cores=ncore)

  if(model == 'tgm1') {
  
    hyhat <- array(dim=c(length(replist[[1]]$hy),nits))
    for(nn in 1:nits) hyhat[,nn] <- replist[[nn]]$hy

  } 

  if(model == 'tgm1re') {
  
    hyhat <- array(dim=c(length(replist[[1]]$hy),nits))
    sigmahat <- rep(NA,nits)
    for(nn in 1:nits) {

      hyhat[,nn] <- replist[[nn]]$hy
      sigmahat[nn] <- replist[[nn]]$sigmah

    }
  } 

  if(model == 'tmg2') {

    Phi.mc <- array(dim=c(nits,nr,nr)) 
    hyhat <- array(dim=c(dim(replist[[1]]$hy,nits)))
    for(nn in 1:nits) {

      Phi.mc[nn,,] <- replist[[nn]]$Phi 
      hyhat[,,nn] <- replist[[nn]]$hy 

    }
  }

  if(model == 'tmg2re') {

    Phi.mc <- array(dim=c(nits,nr,nr)) 
    hyhat <- array(dim=c(dim(replist[[1]]$hy,nits)))
    sigmahat <- rep(NA,nits) 
    for(nn in 1:nits) {

      Phi.mc[nn,,] <- replist[[nn]]$Phi 
      hyhat[,,nn] <- replist[[nn]]$hy 
      sigmahat[nn] <- replist[[nn]]$sigmah 

    }
  } 
  rm(replist)

  # update current year etc

  ycurr <- ydec[j]
  yprj <- ydec[j]+msectrl$tactau 

  Cterm <- prj$C[ycurr,,]
  TACold <- apply(Cterm,2,sum)

  # apply the MP f

  cat("Current year: ",ycurr,"\n")
  cat("TAC decision ",j," of ",ndec,"\n")
  TACnew <- get.TAC(tacctrl)

  # prep Biol and Fleet objects

  Biolx <- Biol
  Fleetx <- Fleet
  Biolx$Ninit <- prj$N[ycurr,,,,]
  Fleetx$Cfix <- tacctrl$catchSplit %o% TACnew

  # project forward

  prjtmp <- get.popdyn(Biolx,Fleetx,msectrl$tactau,CatchControl,nits) 

  # joining it all up

  prj$N <- abind(prj$N,prjtmp$N,along=1)
  prj$H <- abind(prj$H,prjtmp$H,along=1) 
  prj$C <- abind(prj$C,prjtmp$C,along=1) 
  prj$SSB <- abind(prj$SSB,prjtmp$SSB,along=1) 
  prj$SSBtot <- abind(prj$SSBtot,prjtmp$SSBtot,along=1) 
  rm(prjtmp) 

}

cat("MSE loop completed...","\n")

