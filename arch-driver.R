# The file "arch-driver-R". The purpose is to simulate ARCH(1) stochastic process with 
# user defined parameters, and then simulate the "latent" models' parameters using Stran 
# high-order programming language ARCH(1) model. 

library(rstan)   # we load in the Stan package
library(tseries)
setwd("C:/Users/computer/Desktop/thesis")  

# In first step, we generate ARCH(1) process
#
alpha0 <- 0.06321  # user defined subjective constant
alpha1 <- 0.03548  # user defined subjective constant
T <- 10000         # number of replications 

r_t <- rep(0,T)

sigma2_t <- 0.1    # initial guess; user defined parameter
r_t[1] <- sqrt(sigma2_t)*rnorm(1,0,1)

for (i in 2:T) {
    eps <- rnorm(1,0,1)    
    sigma2_t <- alpha0 + alpha1*r_t[i-1]^2
    r_t[i] <- sqrt(sigma2_t)*eps
}

plot.ts(r_t)
hist(r_t)

## Stan modelling in R 
#
modelString = "

    data {
        int<lower=0> T; // number of time points
        real r_t[T]; // return at time t
    }

    parameters {
        real mu;
        real<lower=0> alpha0; // noise intercept
        real<lower=0,upper=1> alpha1; // noise slope
    }

    model {

        for (t in 2:T)
            r_t[t] ~ normal(mu, sqrt(alpha0 + alpha1 * pow(r_t[t-1] - mu, 2)));
    }
"

sm <- stan_model(model_code = modelString) # Stan model 

fit <- sampling(sm, data=list(r_t=r_t, T=T), iter=1000, chains=4, warmup=250, thin=5)

print(fit, digits_summary = 5)
plot(fit)

fit_sim <- extract(fit)
hist(fit_sim$alpha0)
hist(fit_sim$alpha1)
# Done!..
