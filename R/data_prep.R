######################################################
# prep the data ######################################
######################################################
# R. Hillary CSIRO O & A 2019 ########################
######################################################

yrs <- 1985:2020 # back 10 years before catches begin
ny <- length(yrs)
iyrs <- 0:(ny-1) # TMB model index years
names(iyrs) <- as.character(yrs)
ages <- 1:52
na <- length(ages)
nr <- 2 # area 1 (N), area 2 (S)
ns <- 2
lbins <- c(0,10,20,30,seq(35,140,5),seq(150,190,10))
nbins <- length(lbins)-1
nl <- length(lbins)-1
mulbins <- 0.5*(lbins[-c(1)]+lbins[-c(nl+1)])
nf <- 5 # ATT, NVT, ATL, NMRL, SMRL longline
ffnm <- c("ATT","NVT","ATL","NMRL","SMRL")
rf <- c(2,1,2,1,2)-1 # areas fisheries operate in C++ offset to zero
names(rf) <- ffnm
nrecmax <- 9 # maxium years of recaptures post release year

######################################
# set up the data: ###################
# 1. catch biomass ###################
# 2. catch timing (fraction of year) #
# 3. length frequency data ###########
# 4. age-length frequency data #######
# 5. tagging data ####################
######################################

# catch biomass and timing by fishery

datdir <- "../../../data/toothfish2021/"
fnm <- paste(datdir,"CatchBiomass.RData",sep="")
load(fnm)
gnm <- paste(datdir,"CatchTiming.RData",sep="")
load(gnm)

C <- tau <- matrix(0,nrow=ny,ncol=nf)
rownames(C) <- rownames(tau) <- as.character(yrs)
colnames(C) <- colnames(tau) <- ffnm

for(i in 1:dim(cbdf)[1]) {

  yy <- as.character(cbdf$Year[i])
  ff <- as.character(cbdf$Fleet[i])
  cc <- cbdf$Catch_t[i]
  tt <- taudf$tau[i]
  C[yy,ff] <- cc
  tau[yy,ff] <- tt

}

# modify 2001 and 2010 missing catch to avoid zero probs in tag likelihood

C['2001',] <- rep(0.01,nf)
C['2010','NVT'] <- C['2010','NMRL'] <- 0.01  

# fisheries with logistic (0) or double-logistic (1) selectivty functions

self <- c(1,1,0,0,0)
names(self) <- ffnm
nlg <- length(self[self==0])
ndl <- length(self[self==1])
flg <- c(2,3,4)
fdl <- c(0,1)

# length frequency data (using nhauls as initial ESS)

fnm <- paste(datdir,"LF_nhaul.dat",sep="")
lfdat <- matrix(scan(fnm),ncol=nbins+3,byrow=T)

iy <- unname(iyrs[as.character(lfdat[,1])])
iff <- lfdat[,2]-1
LFcov <- unname(cbind(iy,iff))
LF <- lfdat[,-c(1,2)]

# age-given-length frequency data

fnm <- paste(datdir,"ALF.dat",sep="")
alfdat <- matrix(scan(fnm),ncol=na+5,byrow=T)
iy <- unname(iyrs[as.character(alfdat[,1])])
iff <- alfdat[,2]-1
isx <- alfdat[,3]-1
il <- alfdat[,4]-1
ALFcov <- unname(cbind(iy,iff,isx,il))
ALF <- alfdat[,-c(1:4)]

# tagging data

fnm <- paste(datdir,"tags.dat",sep="")
tdf <- read.table(fnm,header=T)
tmin <- 5 # minimum number of releases
tdf <- subset(tdf,T >= tmin)
relyrs <- sort(unique(tdf$relyr))

tlist <- list()
ii <- 1 
for(y in relyrs) {

  yrecmax <- min(max(yrs),y+nrecmax)
  xdf <- subset(tdf,relyr == y & recyr <= yrecmax)
  ll <- sort(unique(xdf$lrel))
  for(l in ll) {

    ydf <- subset(xdf,lrel == l)

    rrel <- unique(ydf$relarea)
    if(length(rrel) == 1) {

      rrel <- ydf$relarea[1] 
      Rtmp1 <- subset(ydf,recarea==1)$R
      Rtmp2 <- subset(ydf,recarea==2)$R 
      tlist[[ii]] <- list()
      tlist[[ii]]$T <- ydf$T[1]
      tlist[[ii]]$yrel <- y
      tlist[[ii]]$l <- l 
      tlist[[ii]]$rrel <- rrel
      tlist[[ii]]$nrec <- yrecmax-y
      tlist[[ii]]$R <- unname(cbind(Rtmp1,Rtmp2))
      ii <-ii+1

    } else {

      for(rr in rrel) {

        zdf <- subset(ydf,relarea == rr)
        Rtmp1 <- subset(zdf,recarea==1)$R
        Rtmp2 <- subset(zdf,recarea==2)$R 
        tlist[[ii]] <- list()
        tlist[[ii]]$T <- zdf$T[1]
        tlist[[ii]]$yrel <- y
        tlist[[ii]]$l <- l 
        tlist[[ii]]$rrel <- rr
        tlist[[ii]]$nrec <- yrecmax-y
        tlist[[ii]]$l <- l 
        tlist[[ii]]$R <- unname(cbind(Rtmp1,Rtmp2))
        ii <-ii+1 

      }
    }
  }
}

ntagev <- length(tlist)
tagrel <- tagnrec <- rep(NA,ntagev)
tagcov <- matrix(ncol=3,nrow=ntagev)
tagrec <- array(0,dim=c(ntagev,nrecmax,nr))
for(i in 1:ntagev) {

  tagrel[i] <- tlist[[i]]$T
  tagnrec[i] <- tlist[[i]]$nrec
  tagcov[i,] <- c(tlist[[i]]$l-1,tlist[[i]]$rrel-1,unname(iyrs[as.character(tlist[[i]]$yrel)]))
  nrectrue <- dim(tlist[[i]]$R)[1]
  tagrec[i,1:nrectrue,] <- tlist[[i]]$R

}

# reporting rates over time

rrdf <- read.table("tagrep.dat",header=T)

# save the data

fnm <- paste(paste("MIdata_maxrecev",nrecmax,sep="_"),"rda",sep=".")
save.image(fnm)


