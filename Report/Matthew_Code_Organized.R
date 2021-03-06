# Organized Code

#-----------------------------#
# Load in the data            #
#-----------------------------#
library(ncdf4)
weather<-nc_open("pnwrain.50km.daily.4994.nc")
#details about the dataset
print(weather)
weatherDat<-ncvar_get(weather, attributes(weather$var)$names[1])
weatherDat
dim(weatherDat) #latitude (17), longitude (16), time (16801 in days since 1949)
#each cross section shows the mm/day amount of rainfall 
#32767 missing data points

# Define the 9 regions
finer_grids <- 
  c( 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     0,0,0,0,0,3,5,5,5,6,6,6,6,7,7,7,7,
     0,0,0,0,0,3,5,5,6,6,6,6,6,7,7,7,7,
     0,0,1,1,3,3,5,5,6,6,6,6,6,6,6,6,6,
     0,0,0,3,3,3,5,5,6,6,6,6,6,6,6,6,6,
     0,0,0,3,3,3,5,5,6,6,6,6,6,6,6,6,6,
     0,0,0,3,3,3,5,5,6,6,6,6,6,6,6,6,6,
     0,0,0,2,3,3,5,6,6,6,6,6,8,8,8,8,8,
     0,0,0,2,3,3,5,9,9,9,9,9,8,8,8,8,8,
     0,0,0,2,4,4,5,9,9,9,9,9,8,8,8,8,8,
     0,0,0,2,4,4,5,9,9,9,9,9,9,9,9,9,9,
     0,0,0,2,4,4,5,9,9,9,9,9,9,9,9,9,9,
     0,0,0,2,4,4,9,9,9,9,9,9,9,9,9,9,9,
     0,0,2,2,4,4,9,9,9,9,9,9,9,9,9,9,9,
     0,0,2,2,4,4,9,9,9,9,9,9,9,9,0,0,0,
     0,0,0,2,4,4,9,9,9,9,9,9,0,0,0,0,0
  )
grid2 <- t(matrix(finer_grids,ncol=17,byrow=T)[,-17])

#
# Kelly model 1 code here
# I'm skipping the algorithm that computes the min/max dates


#-----------------------------------------------------------#
# Model 3.2                                                 #
# Implemented using Metropolis for each of the 9 areas      #
#-----------------------------------------------------------#

# Piecewise Linear Metropolis

# Metropolis for refined model
library(coda)
#tf_weather_data2 <- weatherDat >= .5
# see code in create_9part_dataset.R to manipulate this into a list of lists.
# this uses an intermediate result - i.e. a list of 9 [ list of 46 [vector of 275 T/F]]
# Parameterize it by mindate,minval, a, b
drop_days <-90
mindate_all <- list(rep(0,9))
minval_all<-list(rep(0,9))
a_all <- list(rep(0,9))
b_all <- list(rep(0,9))
obsrain_all <- list(rep(0,9))
prain_all <- list(rep(0,9))
posterior_means <- matrix(nrow=9,ncol=4)
for (l in seq(1,9)){
  data_met <- tf_data2[[l]]
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
  
  S <- 10000
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
  for( i in seq(1,(S+500))){ # 500 warm up iterations
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
  
  mindate_all[[l]] <- MINDATE
  minval_all[[l]] <- MINVAL
  a_all[[l]] <- A
  b_all[[l]] <- B
  obsrain_all[[l]] <- OBSRAIN
  prain_all[[l]] <- PRAIN
  posterior_means[l,] <- c(mean(MINDATE),mean(MINVAL),mean(A),mean(B))
}

# Model 3.2 Analysis

# Plot the distributions of the means
library(ggplot2)
plt_df <- data.frame(region = rep("1",512),DayOfYear = density(mindate_all[[1]]+drop_days,adjust = 1.5)$x,Density = density(mindate_all[[1]]+drop_days,adjust = 1.5)$y)
for(l in seq(2,9)){
  plt_df <- rbind(plt_df,data.frame(region = rep(paste(l),512),DayOfYear = density(mindate_all[[l]]+drop_days,adjust = 1.5)$x,Density = density(mindate_all[[l]]+drop_days,adjust = 1.5)$y))
}

ggplot(data=plt_df, aes(x=DayOfYear, y=Density, group=region) ) +geom_line(aes(color=region)) + ggtitle("Posterior Distributions of the Nicest Day of the Year")

# Credible Intervals
par(mfrow = c(3,3))
for (l in 1:9) {
  plot(seq(drop_days+1,365),apply(matrix(as.numeric(unlist(tf_data2[[l]])),ncol=365,byrow=T),2,sum)[(drop_days+1):365]/46,ylab = "Probability of Rain",xlab = "Day of Year",main = paste("Region",l))
  points(seq(drop_days+1,365),prob_by_day(mean(mindate_all[[l]]),mean(minval_all[[l]]),mean(a_all[[l]]),mean(b_all[[l]])),type="l")
  points(seq(drop_days+1,365),apply(obsrain_all[[l]],2,quantile,.975),type="l",col="red")
  points(seq(drop_days+1,365),apply(obsrain_all[[l]],2,quantile,.025),type="l",col="red")
}
par(mfrow = c(1,1))

#Escape Seattle Rain
escape_p <- seq(1,365-drop_days)
num_escape <- seq(1,365-drop_days)
for( i in seq(1,365-drop_days)){
  searain <- rbinom(length(prain_all[[3]][,i]),1,prain_all[[3]][,i])
  eastdry <- rbinom(length(prain_all[[3]][,i]),1,1-prain_all[[6]][,i])
  escape_p[i] <-  sum(searain*eastdry)/sum(searain)
  num_escape[i] <- sum(searain*eastdry)/length(searain*eastdry)
}
plot(seq(drop_days+1,365),escape_p,ylim = c(0,1),col="blue",cex=.5,xlab = "Day Of Year",ylab = "Probability",main = "Escaping Seattle Rain",pch=16)
points(seq(drop_days+1,365),num_escape,col="red",cex=.5,pch=16)
legend("topright",legend = c("Given Rain in Seattle","Unconditional Probability"),col = c("blue","red"),cex = .5,pch=c(16,16))
