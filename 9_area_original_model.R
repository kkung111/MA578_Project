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
westMins20<-30 #-(245 - westMinMu0)/2
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
nW<-dim(west)[1]
nu0<-2 #note: had to change this because was giving me errors with 1
nun<-nu0 + nW


set.seed(1234+k)

SigmaW<-riwish(nu0, westCov)
thetaW<-rmvnorm(1, westMu0, wests20)
#nSim<-10000

nSim<-20000


THW<-S2W<-YtW<-NULL
for(i in 1:(nSim+1000)){
  #sample Sigma
  LnW<-westCov + crossprod(west - outer(rep(1, nrow(west)), c(thetaW)))
  SigmaW<-riw(1, nun, LnW)[,,1]
  
  
  #sample theta
  LN<-solve(solve(wests20) + nW *solve(SigmaW))
  munW<-LN%*%(solve(wests20)%*%westMu0 + nW*solve(SigmaW)%*%ybW)
  thetaW<-rmv(1, munW,LN)
  
 
  
  #prediction
  if(i > 1000) {
  yPred<-rmv(1, thetaW, SigmaW)
  while(yPred[1]<0 | yPred[2]>365){yPred<-rmv(1, thetaW, SigmaW)}
  ytW<-yPred
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


seqLength<-seq(1, dim(output_yt[[1]])[1], 10)
tempoutput_yt<-output_yt
tempoutput_th<-output_th
tempoutput_s2<-output_s2

output_yt<-lapply(output_yt, function(x){x[seqLength,]})
output_th<-lapply(output_th, function(x){x[seqLength,]})
output_s2<-lapply(output_s2, function(x){x[seqLength,]})
west_dry_season <- rep(0,365)
east_dry_season <- rep(0,365)
west_v_east <- rep(0,365)
for(j in seq(1,365)) {
  west_dry_season[j] <- sum(output_yt[[3]][,1] <= j & output_yt[[3]][,2] > j )/length(output_yt[[3]][,2])
  east_dry_season[j] <- sum(output_yt[[6]][,1] <= j & output_yt[[6]][,2] > j )/length(output_yt[[6]][,2])
  west_v_east[j] <- (sum( (!(output_yt[[3]][,1] <= j & output_yt[[3]][,2] > j))&(output_yt[[6]][,1] <= j & output_yt[[6]][,2] > j ))) / length(seqLength)#(sum(!(output_yt[[3]][,1] <= i & output_yt[[3]][,2] > i)))
  # west_v_east[i] <- (sum(  (output_yt[[3]][,1] > i & output_yt[[6]][,1] <= i  ) | (output_yt[[3]][,2] <= i & output_yt[[6]][,2] > i ))) / length(output_yt[[3]][,1])
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
for(j in seq(1,365)) {
  west_dry_season[j] <- sum(output_yt[[4]][,1] <= j & output_yt[[4]][,2] > j )/length(output_yt[[4]][,2])
  east_dry_season[j] <- sum(output_yt[[9]][,1] <= j & output_yt[[9]][,2] > j )/length(output_yt[[9]][,2])
  west_v_east[j] <- (sum( (!(output_yt[[4]][,1] <= j & output_yt[[4]][,2] > j))&(output_yt[[9]][,1] <= j & output_yt[[9]][,2] > j ))) / length(seqLength)
}
plot(west_dry_season,type="l",main = "Prop of Portland Dry Season")
plot(east_dry_season,type="l",main = "Prop of Eastern Oregon Dry Season")
plot(west_v_east,type="l",main = "Predictive Distribution: \nProbability of going East from Portland and going\n from rainy season to dry season",ylab = "Probability",xlab="Day of Year")
# First/Last area
start_dates <- matrix( unlist(lapply(output_yt,function(x) {x[,1]} )),nrow = length(seqLength),ncol=9,byrow=F)
table(apply(start_dates,1,which.min))/length(seqLength)
  
end_dates <- matrix( unlist(lapply(output_yt,function(x) {x[,2]} )),nrow = length(seqLength),ncol=9,byrow=F)
table(apply(start_dates,1,which.max))/length(seqLength)


# Length
summer_length <- matrix( unlist(lapply(output_yt,function(x) {x[,2]-x[,1]} )),nrow = length(seqLength),ncol=9,byrow=F)
round(apply(summer_length,2,mean))
round(sqrt(apply(summer_length,2,var)))


#plot Seattle Data
plot(density(as.matrix(read.table(paste("Area",4,"MinMaxDataV2.tsv"), sep = "\t"))[,1],kernel = "epanechnikov"),lty = 1,xlim=c(1,365),main = "Region 4 Data Distribution",xlab = "Day of Year",ylim = c(0,1/50))
points(density(as.matrix(read.table(paste("Area",4,"MinMaxDataV2.tsv"), sep = "\t"))[,2],kernel = "epanechnikov"),lty = 2,type="l")
legend("topleft", c("Dry Start","Wet Start"),lty = c(1,2))

plot(density(output_yt[[4]][,1]))
points(density(output_yt[[9]][,1]))

# 95% summer start date
round(unlist(lapply( output_yt,function(x) {quantile(x[,1],.95)})))

#plot some of the acf plots
par(mfrow=(c(2,2)))
acf(tempoutput_th[[3]][,1], main = expression(paste("ACF for ", theta, " for Seattle Area")))
acf(output_th[[3]][,1], main = expression(paste("ACF for Thinned ", theta, " for Seattle Area")))
acf(tempoutput_s2[[3]][,1], main = expression(paste("ACF for ", Sigma, " for Seattle Area")))
acf(output_s2[[3]][,1], main = expression(paste("ACF for Thinned ", Sigma, " for Seattle Area")))

#effective size calculation
effSizeth<-matrix(unlist(lapply(tempoutput_th, function(x){effectiveSize(x)})), ncol = 2, byrow = T)
effSizes2<-matrix(unlist(lapply(tempoutput_s2, function(x){effectiveSize(x)})), ncol = 2, byrow = T)

colMeans(effSizeth)
colMeans(effSizes2)

