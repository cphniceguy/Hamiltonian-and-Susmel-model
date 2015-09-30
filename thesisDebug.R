rm(list=ls())

library(tseries)
alpha0 <- 0.025 # user defined subjective constant
alpha1 <- 0.057 # user defined subjective constant
beta0 <- beta1 <- 0.005
T <- 10000       # number of replications 
r_t <- rep(0,T)
sigma2_t <- 1

for (i in 2:T) {
    eps <- rnorm(1,0,1)    
    r_t[i] <- sqrt(sigma2_t)*eps
    sigma2_t <- alpha0 + alpha1*r_t[i-1]^2
}

r_t <- r_t^2; r_t
hist(r_t)
r2_t <- r_t^2; r2_t
hist(r2_t)

plot.ts(r_t)

log_posterior <- function(alpha0,beta0,alpha1,beta1,r2_t) {
    sigma2_t <- alpha0+alpha1*r2_t[-length(r2_t)]
    return(sum(-0.5*log(sigma2_t)-r2_t[-1]/(2*sigma2_t))-alpha0/beta0-alpha1/beta1)
} 

# input: theta parameter is 2 dim vector, contains dynamic alpha hyperparameters
#
log_target <- function(theta, r2_t) { # theta, r2_t
    beta0 <- 0.0001  # a constant
    beta1 <- 0.0002 # a constant
    return(log_posterior(theta[1], beta0, theta[2], beta1, r2_t)) 
}

# our proposal function(s)
#
proposal_function <- function(lastTheta, nextTheta) {
    r1 <- dexp(nextTheta,1/lastTheta)
    r2 <- dexp(lastTheta,1/nextTheta)
    if (is.na(r1) || is.na(r2)) {
        browser()
    }
    if (is.na(r1/r2)) {
        browser()
    }
    return(r1/r2)
}

oldTheta <- c(alpha0,alpha1)
newTheta <- oldTheta

numSteps <- 2000  # replications
burnin <- 150    # phase 
mu_alpha0_GE <- 0.00025
mu_alpha1_GE <- 0.5
oldTheta <- c(mu_alpha0_GE,mu_alpha1_GE) # a position in a previous state
tryTheta <- oldTheta

# main loop GE, begin to collect info about draws from the log_target sample
#
accept_rate <- c(0,0) # rate of acceptance of successive moves
targetSample_alpha0 <- numeric(numSteps)
targetSample_alpha1 <- numeric(numSteps)
#
for (i in burnin+1:numSteps){
    
    accept_draw <- runif(1,0,1) # Generate an innovation
    tryTheta[1] <- oldTheta[1]*rexp(1,1) # generate a draw for alpha0
    p_accept <- min(1,proposal_function(oldTheta[1],tryTheta[1])*exp(log_target(tryTheta,r2_t)-log_target(oldTheta,r2_t)))
    #p_accept <- min(1,1/proposal_function(oldTheta[1],tryTheta[1])/exp(log_target(tryTheta,r2_t)-log_target(oldTheta,r2_t)))
    
    #print(p_accept)
    if (accept_draw <= p_accept) { # accept new move
        oldTheta[1] <- tryTheta[1] # update new theta value
        accept_rate[1] <- accept_rate[1] + 1 # count good ones    
    }
    
    targetSample_alpha0[i] <- oldTheta[1] #assign new value of accepted move
    
    #tryTheta <- oldTheta
    
    ####################################################################################
    
accept_draw <- runif(1,0,1)
    tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1
    pfr = proposal_function(oldTheta[2], tryTheta[2])
    if (is.infinite(pfr))
        pfr = proposal_function(oldTheta[2], tryTheta[2])
    lt1 = log_target(oldTheta,r2_t)
    if (is.na(lt1)) {
        browser()
    }
    lt2 = log_target(tryTheta,r2_t)
    if (is.na(lt2))
        browser()
    e1 = exp(lt1-lt2)
    if (is.na(lt1))
        browser()
    p_accept <- min(1,pfr*exp(lt2-lt1))
    if (is.na(accept_draw))
        browser()
        if (is.na(p_accept))
            browser()
    if (accept_draw <= p_accept) {
        oldTheta[2] <- tryTheta[2]
        #accept_rate[2] <- accept_rate[2] + 1        
    }
    
    targetSample_alpha1[i] <- oldTheta[2] #assign values
    
    #tryTheta <- oldTheta
    
}#end main loop 



par(mfrow=c(2,1))
par(mfrow=c(1,1))
plot(targetSample_alpha0,xlab="Iterations",ylab=expression(alpha[0]),type='l', main="Fubini")
plot(targetSample_alpha1,xlab="Iterations",ylab=expression(alpha[1]),type='l', main="Fibonacci")

par(mfrow=c(2,1))
hist(targetSample_alpha0,40,xlab="Iterations",ylab=expression(alpha[0]), main="Fubini")
hist(targetSample_alpha1,40,xlab="Iterations",ylab=expression(alpha[1]), main="Fibonnaci")

summary(targetSample_alpha0)
targetSample_alpha1[1850]
targetSample_alpha1[min(which(targetSample_alpha1 != 0))]


a <- c(1,2,3,4)
a[-1]

## matrix_function_to_do_some_hairy_algebra goes here!..

# done