# Metropolis for refined model
library(coda)
# pi-wet (beta), pi-dry (Beta), dry-start (normal mu_s), dry-end (normal)
#tf_weather_data2 <- weatherDat >= .5
# see code in create_9part_dataset.R to manipulate this into a list of lists.
# this uses an intermediate result - i.e. a list of 9 [ list of 46 [vector of 365 T/F]]
data_met <- tf_data2[[4]]
#priors
wet_a <- 3
wet_b <- 1
dry_a <- 1
dry_b <- 3

theta0 <- c(180,270)
Sigma0 <- matrix(c(30^2,0,0,30^2),ncol=2)

S <- 500
theta <- theta0
pi <- c(.7,.3)
THETASTAR <- THETA <- matrix(ncol=2,nrow = S)
PISTAR <- PI <-  matrix(ncol=2,nrow = S)
accept_probs <- rep(0,S)

delta <- .025
llike <- function(y,pi,Dry){
  # takes in the list of 365 data points
dry <- sapply(Dry,round)

 return( dbinom(sum(y[1:(dry[1]-1)]),dry[1]-1,pi[1],log=T) + 
    dbinom(sum(y[seq(dry[1],dry[2]-1)]) , dry[2]-dry[1] , pi[2],log=T) +
             dbinom(sum(y[seq(dry[2],365)]),366-dry[2], pi[1] ,log=T) )
}

accept <- 0
for( i in seq(1,S)){
  theta.star <- rmvnorm(1,theta,Sigma0*delta)
  if(theta.star[1] <= 1) {theta.star[1] <- 2+abs(theta.star[1])}
  if(theta.star[2] >= 364) {theta.star[2] <- 364-abs(364-theta.star[2])}
  while(theta.star[1] > theta.star[2]) {
    theta.star <- rmvnorm(1,theta,Sigma0*delta)
    if(theta.star[1] <= 1) {theta.star[1] <- 2+abs(theta.star[1])}
    if(dry[2] <= 364) {theta.star[2] <- 364-abs(364-theta.star[2])}
  }
  #pi should be alpha/(alpha+beta)
  #pi.star <- rbeta(2, pi*beta_dist_total, (rep((1-pi)*beta_dist_total,2)))
  pi.star <- 1-abs( 1-abs(rmvnorm(1,pi,delta * matrix(c(.05^2,0,0,.05^2),ncol=2))))
 # print(theta.star)
 # print(pi.star)
  accept_probs[i] <- exp(log.r)
  log.r <- sum(unlist(lapply(data_met,llike,pi=pi.star,Dry=theta.star))) - sum(unlist(lapply(data_met,llike,pi=pi,Dry=theta))) +
 dmvnorm(theta.star,theta0,Sigma0,log=T) + dbeta(pi.star[1],wet_a,wet_b,log=T) + dbeta(pi.star[2],dry_a,dry_b,log=T) -
    dmvnorm(theta,theta0,Sigma0,log=T) - dbeta(pi[1],wet_a,wet_b,log=T) -  dbeta(pi[2],dry_a,dry_b,log=T)
  accept_probs[i] <- exp(log.r)
  if(log(runif(1))<log.r) {
    
    #print("Accept")
    theta <- theta.star
    pi <- pi.star
    accept <- accept + 1
  }
  THETA[i,] <- theta
  PI[i,] <- pi
  THETASTAR[i,] <- theta.star
  PISTAR[i,] <- pi.star
}
sum(accept)/S
apply(THETA,2,effectiveSize)
apply(PI,2,effectiveSize)
apply(THETA,2,mean)
apply(PI,2,mean)
apply(THETA,2,var)
apply(PI,2,var)

plot(density(THETA[seq(100,S,100),1]))
plot(density(THETA[,2]))
plot(density(PI[,1]))
plot(density(PI[,2]))

#Analysis
# Rain freq by day
data_matrix <- matrix(unlist(data_met),ncol=365,byrow=T)
round(apply(data_matrix,2,sum)/46,2)
plot(round(apply(data_matrix,2,sum)/46,2),xlab = "Day of Year" , ylab = "Probability of Precipitation", main = "Area 4 Precipitation\nChance by day of Year")
