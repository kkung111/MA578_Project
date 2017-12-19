#Kelly Kung

#metroplis Hastings Alg

#prior --> N()
#Sigma--> sample var/covariance matrix nu0 = 1

#we are gonna code up the metropolis algorithm and see what happens

#next steps--> maybe try making grid into every day?
#setwd("/Users/kkung/Documents/GitHub/MA578_Project")
library(MCMCpack)
library(mvtnorm)

output_th <- list(rep(0,9))
output_s2 <- list(rep(0,9))
output_yt <- list(rep(0,9))
for (k in seq(1,9)) {
west<-as.matrix(read.table(paste("Area",k,"MinMaxDataV2.tsv"), sep = "\t"))


westCov<-cov(west)




#find the prior means
westMinMu0<-170
westMaxMu0<-260


westMu0<-matrix(c(westMinMu0, westMaxMu0), nrow = 2)


#find the prior standard deviations
westMins20< 30 #-(245 - westMinMu0)/2
westMaxs20<-30 #(303 - westMaxMu0)/2

wests20<-matrix(c(westMins20, 0, 0, westMaxs20), nrow = 2)


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



#starting values
ybW<-apply(west, 2, mean)
k0<-1
nW<-dim(west)[1]
nu0<-2 #note: had to change this because was giving me errors with 1
nun<-nu0 + n

kn<-k0 + n

set.seed(1234)

SigmaW<-riwish(nu0, westCov)
thetaW<-rmvnorm(1, westMu0, wests20)
nSim<-20000

THW<-S2W<-YtW<-NULL
for(i in 1:nSim+1000){
  #sample Sigma
  LnW<-westCov + crossprod(west - outer(rep(1, nrow(west)), c(ybW))) + k0*n/kn*
    crossprod(t(westMu0 - ybW))
  SigmaW<-riw(1, nun, LnW)[,,1]
  
  
  #sample theta
  munW<-(k0*westMu0 + nW *ybW)/kn
  thetaW<-rmv(1, munW, SigmaW/kn)
  
 
  
  #prediction
  if(i > 1000) {
  ytW<-rmv(1, thetaW, SigmaW)
  THW<-cbind(THW, thetaW)
  S2W<-cbind(S2W, c(SigmaW))
  YtW<-cbind(YtW, ytW)
  }

}

rowMeans(THW)
rowMeans(S2W)
rowMeans(YtW)

output_th[[k]] <- t(THW)
output_s2[[k]] <- t(S2W)
output_yt[[k]] <- t(YtW)
}
#analysis

west_dry_season <- rep(0,365)
east_dry_season <- rep(0,365)
west_v_east <- rep(0,365)
for(i in seq(1,365)) {
  west_dry_season[i] <- sum(output_yt[[3]][,1] <= i & output_yt[[3]][,2] > i )/length(output_yt[[3]][,2])
  east_dry_season[i] <- sum(output_yt[[6]][,1] <= i & output_yt[[6]][,2] > i )/length(output_yt[[6]][,2])
  west_v_east[i] <- (sum(  (output_yt[[3]][,1] > i & output_yt[[6]][,1] <= i  ) | (output_yt[[3]][,2] <= i & output_yt[[6]][,2] > i ))) / length(output_yt[[3]][,1])
}
plot(west_dry_season,type="l",main = "Prop of Seattle Dry Season")
plot(east_dry_season,type="l",main = "Prop of Columbia Plateau Dry Season")
plot(west_v_east,type="l",main = "Predictive Distribution: \nProbability of going East from Seattle and going\n from rainy season to dry season",ylab = "Probability",xlab="Day of Year")

# Means
round(matrix(unlist(lapply(output_th,function(x) {apply(x,2,mean)})),ncol=2,byrow=T))
round(matrix(unlist(lapply(output_th,function(x) {apply(x,2,var)})),ncol=2,byrow=T))


west_dry_season <- rep(0,365)
east_dry_season <- rep(0,365)
west_v_east <- rep(0,365)
for(i in seq(1,365)) {
  west_dry_season[i] <- sum(output_yt[[4]][,1] <= i & output_yt[[4]][,2] > i )/length(output_yt[[4]][,2])
  east_dry_season[i] <- sum(output_yt[[9]][,1] <= i & output_yt[[9]][,2] > i )/length(output_yt[[9]][,2])
  west_v_east[i] <- (sum(  (output_yt[[4]][,1] > i & output_yt[[9]][,1] <= i  ) | (output_yt[[4]][,2] <= i & output_yt[[9]][,2] > i ))) / length(output_yt[[4]][,1])
}
plot(west_dry_season,type="l",main = "Prop of Portland Dry Season")
plot(east_dry_season,type="l",main = "Prop of Eastern Oregon Dry Season")
plot(west_v_east,type="l",main = "Predictive Distribution: \nProbability of going East from Portland and going\n from rainy season to dry season",ylab = "Probability",xlab="Day of Year")
# First/Last area
start_dates <- matrix( unlist(lapply(output_yt,function(x) {x[,1]} )),nrow = nSim,ncol=9,byrow=F)
table(apply(start_dates,1,which.min))/nSim
  
end_dates <- matrix( unlist(lapply(output_th,function(x) {x[,2]} )),nrow = nSim,ncol=9,byrow=F)
table(apply(start_dates,1,which.max))/nSim


# Length
summer_length <- matrix( unlist(lapply(output_yt,function(x) {x[,2]-x[,1]} )),nrow = nSim,ncol=9,byrow=F)
round(apply(summer_length,2,mean))
round(sqrt(apply(summer_length,2,var)))

#plot Seattle Data
plot(density(as.matrix(read.table(paste("Area",4,"MinMaxDataV2.tsv"), sep = "\t"))[,1],kernel = "epanechnikov"),lty = 1,xlim=c(1,365),main = "Area 4 Data Distribution",xlab = "Day of Year",ylim = c(0,1/50))
points(density(as.matrix(read.table(paste("Area",4,"MinMaxDataV2.tsv"), sep = "\t"))[,2],kernel = "epanechnikov"),lty = 2,type="l")
legend("topleft", c("Dry Start","Wet Start"),lty = c(1,2))