# Piecewise Linear Metropolis

# Metropolis for refined model
library(coda)
# pi-wet (beta), pi-dry (Beta), dry-start (normal mu_s), dry-end (normal)
#tf_weather_data2 <- weatherDat >= .5
# see code in create_9part_dataset.R to manipulate this into a list of lists.
# this uses an intermediate result - i.e. a list of 9 [ list of 46 [vector of 275 T/F]]
# Parameterize it by mindate,minval, a, b
drop_days <-90
data_met <- tf_data2[[6]]
year_aggregates <-  apply(matrix(as.numeric(unlist(data_met)),ncol=365,byrow=T),2,sum)[(drop_days+1):365]

#priors
mindate_mean <- mindate <- 210-drop_days
mindate_sd <- 30
minval_a  <- 5
minval_b <- 20
a_mean <-  -.6/150
a_sd <- abs(a_mean)/2
b_mean <- b <- .6/150
b_sd <- abs(b_mean)/2
minval <- minval_a / (minval_a + minval_b)

S <- 4000
MINDATE <- MINVAL <- A <- B <- rep(0,S)
OBSRAIN <- PRAIN <- matrix(ncol=365-drop_days,nrow=S)
accept_probs <- rep(0,S+500)

delta <- .095
llike <- function(y,theta_min,theta_val,a,b){
  # takes in the list of 365 data points
  after_min <- c(rep(0,floor(theta_min)),rep(1,365-drop_days-floor(theta_min)))
  j <- seq(1,(365-drop_days))
  prob <- theta_val + (j-theta_min)*(a*(1- after_min)+ b*after_min)
  prob <- sapply(prob,function(x) {max(min(x,1),0)})
  llike <- sum(dbinom(y,46,prob,log=T))
  return(llike)
}

prob_by_day <- function(theta_min,theta_val,a,b)
{
  after_min <- c(rep(0,floor(theta_min)),rep(1,365-drop_days-floor(theta_min)))
  j <- seq(1,(365-drop_days))
  prob <- theta_val + (j-theta_min)*(a*(1- after_min)+ b*after_min)
  prob <- sapply(prob,function(x) {max(min(x,1),0)})
  return(prob)
}

accept <- 0
for( i in seq(1,(S+500))){
 mindate.star <- rnorm(1,mindate,mindate_sd*delta)
 minval.star <- rnorm(1,minval,sqrt(minval_a*minval_b/((minval_a+minval_b)^2* (minval_a+minval_b+1)))*delta)
 a.star <- rnorm(1,a,delta*a_sd)
 b.star <- rnorm(1,b,delta*b_sd)
 
  log.r <- llike(year_aggregates,mindate.star,minval.star,a.star,b.star) -llike(year_aggregates,mindate,minval,a,b) +
    dnorm(a.star,a_mean,a_sd,log=T) + dbeta(minval.star,minval_a,minval_b,log=T) + dnorm(a.star,a_mean,a_sd,log=T) + dnorm(b.star,b_mean,b_sd,log=T) -
    dnorm(a,a_mean,a_sd,log=T) - dbeta(minval,minval_a,minval_b,log=T) - dnorm(a,a_mean,a_sd,log=T) - dnorm(b,b_mean,b_sd,log=T) 
  accept_probs[i] <- exp(log.r)
  if(log(runif(1))<log.r){
    if(i > 500) accept <- accept+1
    #print("Accept")
  a <- a.star
    b <- b.star
    minval <- minval.star
    mindate <- mindate.star
  }
  output_filter <- 50
  if((i > 500) ) {
A[i-500] <- a
B[i-500] <- b
MINDATE[i-500] <- mindate
MINVAL[i-500] <- minval
PRAIN[i-500,] <- prob_by_day(mindate,minval,a,b)
OBSRAIN[i-500,] <- rbinom(365-drop_days,46,prob_by_day(mindate,minval,a,b))/46
}
}
MINDATE <- MINDATE[seq(1,S,by=50)]
MINVAL <- MINVAL[seq(1,S,by=50)]
A <- A[seq(1,S,by=50)]
B <- B[seq(1,S,by=50)]
OBSRAIN <- OBSRAIN[seq(1,S,by=50),]
PRAIN <- PRAIN[seq(1,S,by=50),]

sum(accept)/S
c(effectiveSize(MINDATE),effectiveSize(MINVAL),effectiveSize(A),effectiveSize(B))
c(mean(MINDATE),mean(MINVAL),mean(A),mean(B))
c(var(MINDATE),var(MINVAL),var(A),var(B))
plot(density(MINDATE+drop_days,adjust = 1.5),xlab = "Day of Year", main = "Distribution of Nicest day of the Year")
plot(seq(drop_days+1,365),year_aggregates/46,ylab = "Probability of Rain",xlab = "Day of Year",main = "Probability of Rain by Day\nWith 95% Credible Interval")
points(seq(drop_days+1,365),prob_by_day(mean(MINDATE),mean(MINVAL),mean(A),mean(B)),type="l")
points(seq(drop_days+1,365),apply(OBSRAIN,2,quantile,.975),type="l",col="red")
points(seq(drop_days+1,365),apply(OBSRAIN,2,quantile,.025),type="l",col="red")


#Ad hoc stuff
#prain3 <- PRAIN
#prain6 <- PRAIN

     

#Analysis Ideas:
# Biggest Change

#P escape rain by day
p_escape_rain <- rep(0,365-drop_days)
for(k in 1:(365-drop_days)){
  p_escape_rain[k] <- mean(sapply(prain3[,k],rbinom,n=1,size=1)&sapply(1-prain6[,k],rbinom,n=1,size=1))
}
plot(seq(drop_days+1,365) , p_escape_rain , xlab = "day of year",ylab = "Probability of Escaping Rain",main = "Daily Chance that Area 3 rains and Area 6 is Dry") 

#7 constant days of rain
week_of_rain = rep(0,365-drop_days-7)
for(k in 1:(365-drop_days-7)){
  week_of_rain[k] <- mean(sapply(prain3[,k],rbinom,n=1,size=1)*sapply(prain3[,k+1],rbinom,n=1,size=1)
  *sapply(prain3[,k+2],rbinom,n=1,size=1)
  *sapply(prain3[,k+3],rbinom,n=1,size=1)
  *sapply(prain3[,k+4],rbinom,n=1,size=1)
  *sapply(prain3[,k+5],rbinom,n=1,size=1)
  *sapply(prain3[,k+6],rbinom,n=1,size=1))
}
