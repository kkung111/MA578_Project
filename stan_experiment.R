#STAN Code
#library(devtools)
#install.packages("devtools")
#
#library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_code <- '
data {
int<lower=0, upper=46> y[365];       // data  
}
parameters {
real<lower=0, upper=365> drystart;
real<lower=drystart, upper=365> wetstart;
real<lower=0, upper=1> wet_prob;
real<lower=0, upper=1> dry_prob;

}


model {

drystart ~ normal(180,30);
wetstart ~ normal(270,30);
wet_prob ~ beta(3,1);
dry_prob ~ beta(1,3);

for (j in 1:365) {
real p;
if (j < drystart) 
p = wet_prob;
else if (j >= wetstart) 
p = wet_prob;
else 
p = dry_prob;
y ~ binomial(46,p);

}
}
'
standat <- apply(matrix(as.numeric(unlist(data_met)),ncol=365,byrow=T),2,sum)
fit <- stan(model_code = stan_code, data = list(y = standat), iter = 1000, chains = 4,control = list(max_treedepth = 15))
# If nothing seems to kickoff
#.rs.restartR()

#sum(sapply(standat[1:201],dbinom,46,.49,log=T) + sum(sapply(standat[202:203],dbinom,46,.47,log=T)) +sum(sapply(standat[204:365],dbinom,46,.49,log=T)))



