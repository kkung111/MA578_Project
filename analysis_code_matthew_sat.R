#Kelly Kung

#metroplis Hastings Alg

#prior --> N()
#Sigma--> sample var/covariance matrix nu0 = 1

#we are gonna code up the metropolis algorithm and see what happens

#next steps--> maybe try making grid into every day?
#setwd("/Users/kkung/Documents/GitHub/MA578_Project")
east<-as.matrix(read.table("AllEastMinMax.tsv", sep = "\t"))
west<-as.matrix(read.table("AllWestMinMax.tsv", sep = "\t"))

eastCov<-cov(east)
westCov<-cov(west)

##west min: 149, 245
##west max: 224 303
##east min: 116, 249
##east max: 238, 304

#find the prior means
westMinMu0<-mean(c(149, 245))
westMaxMu0<-mean(c(224, 303))
eastMinMu0<-mean(c(116, 249))
eastMaxMu0<-mean(c(238, 304))

westMu0<-matrix(c(westMinMu0, westMaxMu0), nrow = 2)
eastMu0<-matrix(c(eastMinMu0, eastMaxMu0), nrow = 2)

#find the prior standard deviations
westMins20<-(245 - westMinMu0)/2
westMaxs20<-(303 - westMaxMu0)/2
eastMins20<-(249 - eastMinMu0)/2
eastMaxs20<-(304 - eastMaxMu0)/2

wests20<-matrix(c(westMins20, 0, 0, westMaxs20), nrow = 2)
easts20<-matrix(c(eastMins20, 0, 0, eastMaxs20), nrow = 2)

#do the Gibbs Sampler First just so we have it

#helper functions from class
rmv<-function(n,mu,Sigma){    # samples Y~MVN(mu,Sigma)
  cm<-chol(Sigma);d<-dim(Sigma)[1]
  Y0<-matrix(rnorm(n*d),nrow=d)
  t(cm)%*%Y0 + mu
}
riw<-function(n,nu0,Sm){ # Sigma~IW(nu0,Sigma^(-1)); requires rmv
  m<- solve(Sm)
  sapply(1:n,function(i) 
    solve(crossprod(t(rmv(nu0*n,0,m))[(i-1)*nu0+1:nu0,])), simplify = 'array')
}

library(MCMCpack)
library(mvtnorm)

#starting values
ybW<-apply(west, 2, mean)
ybE<-apply(east, 2, mean)
k0<-1
nW<-dim(west)[1]
nE<-dim(east)[1]
nu0<-2 #note: had to change this because was giving me errors with 1
nun<-nu0 + n

kn<-k0 + n

set.seed(1234)
SigmaE<-riwish(nu0, eastCov)
thetaE<-rmvnorm(1, eastMu0, easts20)

SigmaW<-riwish(nu0, westCov)
thetaW<-rmvnorm(1, westMu0, wests20)
nSim<-5000

THW<-S2W<-YtW<-NULL
THE<-S2E<-YtE<-NULL
for(i in 1:nSim){
  #sample Sigma
  LnW<-westCov + crossprod(west - outer(rep(1, nrow(west)), c(ybW))) + k0*n/kn*
    crossprod(t(westMu0 - ybW))
  SigmaW<-riw(1, nun, LnW)[,,1]
  
  LnE<-eastCov + crossprod(east - outer(rep(1, nrow(east)), c(ybE))) + k0*n/kn*
    crossprod(t(eastMu0 - ybE))
  SigmaE<-riw(1, nun, LnE)[,,1]
  
  #sample theta
  munW<-(k0*westMu0 + nW *ybW)/kn
  thetaW<-rmv(1, munW, SigmaW/kn)
  
  munE<-(k0*eastMu0 + nE *ybE)/kn
  thetaE<-rmv(1, munE, SigmaE/kn)
  
  #prediction
  ytW<-rmv(1, thetaW, SigmaW)
  THW<-cbind(THW, thetaW)
  S2W<-cbind(S2W, c(SigmaW))
  YtW<-cbind(YtW, ytW)
  
  ytE<-rmv(1, thetaE, SigmaE)
  THE<-cbind(THE, thetaE)
  S2E<-cbind(S2E, c(SigmaE))
  YtE<-cbind(YtE, ytE)
}

rowMeans(THW)
rowMeans(S2W)
rowMeans(YtW)

rowMeans(THE)
rowMeans(S2E)
rowMeans(YtE)

#How likely the the east side to be dry and west side wet
west_dry_season <- rep(0,365)
east_dry_season <- rep(0,365)
west_v_east <- rep(0,365)
for(i in seq(1,365)) {
  west_dry_season[i] <- sum(YtW[1,] <= i & YtW[2,] > i )/length(YtW[2,])
  east_dry_season[i] <- sum(YtE[1,] <= i & YtE[2,] > i )/length(YtE[2,])
  west_v_east[i] <- (sum(  (YtW[1,] > i & YtE[1,] <= i  ) | (YtW[2,] <= i & YtE[2,] > i ))) / length(YtW[1,])
}
plot(west_dry_season,type="l",main = "Prop of West Dry Season")
plot(east_dry_season,type="l",main = "Prop of East Dry Season")
plot(west_v_east,type="l",main = "posterior prob of going East for dry season",ylab = "prop",xlab="day of year")


#now try to do this with metropolis
#need to propose a variance --> want acceptance between 30 - 40
thetaW<-matrix(c(0,0), nrow = 2)
thetaE<-matrix(c(0,0), nrow = 2)
tuningparam<-seq(0, .1, .01)
#tuningparam<-.04
accTuningW<-c()
accTuningE<-c()
#or(j in 1:length(tuningparam)){
proposeVarW<-tuningparam[5]*westCov
proposeVarE<-tuningparam[2]*eastCov
THETAW<-NULL
THETAE<-NULL
accW<-0
accE<-0
set.seed(1234)
for(i in 1:nSim){
  theta.starW<-rmvnorm(1, ybW, proposeVarW)
  log.rW<-sum(dmvnorm(west, theta.starW, westCov, log = TRUE)) + 
    dmvnorm(theta.starW, westMu0, wests20, log = TRUE) - 
    sum(dmvnorm(west, thetaW, westCov, log = TRUE)) - 
    dmvnorm(c(thetaW), westMu0, wests20, log = TRUE)
  
  if(log(runif(1))<log.rW){
    thetaW<-theta.starW
    accW<-accW + 1
  }
  THETAW<-rbind(THETAW, thetaW)
  
  theta.starE<-rmvnorm(1, ybE, proposeVarE)
  log.rE<-sum(dmvnorm(east, theta.starE, eastCov, log = TRUE)) + 
    dmvnorm(theta.starE, eastMu0, easts20, log = TRUE) - 
    sum(dmvnorm(east, thetaE, eastCov, log = TRUE)) - 
    dmvnorm(c(thetaE), eastMu0, easts20, log = TRUE)
  
  if(log(runif(1))<log.rE){
    thetaE<-theta.starE
    accE<-accE + 1
  }
  THETAE<-rbind(THETAE, thetaE)
}
accTuningW<-c(accTuningW, accW)
accTuningE<-c(accTuningE, accE)
#}

accTuningW/nSim
accTuningE/nSim
#goal was to have between 30-40 so I choose index5 for W, 2 for E

#now got to do this for Sigma


#sample from the posterior 
#west
meanPostW<-colMeans(THETAW)
meanPostE<-colMeans(THETAE)


