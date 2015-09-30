# Hacked by Serge Boris Nesterenko

rm(list=ls()) 

setwd("C:/Users/computer/Desktop/thesis")
set.seed(55555)

library(tseries)
library(coda)

load(file="BAC.RData")
load(file="GE.RData")

BAC$Return  <- diff(log(BAC$AdjClose))
GE$Return   <- diff(log(GE$AdjClose))
BAC$atReturn <- BAC$Return - mean(BAC$Return, na.rm=TRUE)
r2_t_BAC <-  na.omit(GE$atReturn^2)
GE$atReturn <- GE$Return - mean(GE$Return, na.rm=TRUE)      
r2_t_GE <-  na.omit(GE$atReturn^2)

MSFT   <- get.hist.quote('MSFT', start='2005-01-01', end='2014-11-04', quote='AdjClose')
plot(MSFT)
MSFT$Return  <- diff(log(MSFT$AdjClose))
MSFT$atReturn <- MSFT$Return - mean(MSFT$Return, na.rm=TRUE)
r2_t_MSFT <-  na.omit(MSFT$atReturn^2)

# The Posteriori Distibution (i.e., the loglikelihood function; the unscaled target density)
#
log_posterior <- function(alpha0,beta0,alpha1,beta1,r2_t) {
    sigma2_t <- alpha0+alpha1*r2_t[-length(r2_t)]
    return(sum(-0.5*log(sigma2_t)-r2_t[-1]/(2*sigma2_t))-alpha0/beta0-alpha1/beta1)
} 
# input: theta parameter is 2 dim vector, contains dynamic alpha hyperparameters
log_target <- function(theta, r2_t) { # theta, r2_t
    beta0 <- 0.0001  # a constant
    beta1 <- 0.0002 # a constant
    return(log_posterior(theta[1], beta0, theta[2], beta1, r2_t)) 
}

proposal_function <- function(lastTheta, nextTheta) {
    dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)    
}


proposal_fun_next <- function(lastTheta, nextTheta) {
    dexp(nextTheta,1/lastTheta)
}

proposal_fun_old <- function(lastTheta, nextTheta) {
    dexp(lastTheta,1/nextTheta)
}


proposal_function_1 <- function(lastTheta, nextTheta) {
    dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
}

proposal_function_2 <- function(lastTheta, nextTheta) {
    dexp(nextTheta,lastTheta)/dexp(lastTheta,nextTheta) 
}

proposal_function_3 <- function(lastTheta, nextTheta) {
    dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta)    
}

proposal_function_4 <- function(lastTheta, nextTheta) {
    dexp(lastTheta,nextTheta)/dexp(nextTheta,lastTheta)    
}

proposal_function_5 <- function(lastTheta, nextTheta) {
    1/(dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta))
}

proposal_function_6 <- function(lastTheta, nextTheta) {
    1/(dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta))    
}

proposal_function_7 <- function(lastTheta, nextTheta) {
    1/(dexp(lastTheta,nextTheta)/dexp(nextTheta,lastTheta)) 
}

proposal_function_8 <- function(lastTheta, nextTheta) {
    1/dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta)    
}


########################################################################################################
########################################################################################################
########################################################################################################
    
    numSteps <- 10000  # replications
    burnin <- 1000    # phase 
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
        p_accept <- min(1,1/proposal_function(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_MSFT)-log_target(tryTheta,r2_t_MSFT)))
        
        #p_accept <- min(1,1/proposal_fun_old(oldTheta[1],tryTheta[1])/proposal_fun_next(oldTheta[1],tryTheta[1])
        #                /exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        #p_accept <- min(1,(1/proposal_function_1(oldTheta[1],tryTheta[1]))*exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        #p_accept <- min(1,(1/proposal_function_1(oldTheta[1],tryTheta[1]))/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))        
        #p_accept <- min(1,proposal_function_2(oldTheta[1],tryTheta[1])*exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_3(oldTheta[1],tryTheta[1])*exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_4(oldTheta[1],tryTheta[1])/exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        
        #p_accept <- min(1,proposal_function_5(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_6(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_7(oldTheta[1],1/tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_8(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))

        #print(p_accept)
        if (accept_draw <= p_accept) { # accept new move
            oldTheta[1] <- tryTheta[1] # update new theta value
            accept_rate[1] <- accept_rate[1] + 1 # count good ones    
        }
        
        targetSample_alpha0[i] <- oldTheta[1] #assign new value of accepted move
        
        ####################################################################################
        
        accept_draw <- runif(1,0,1)
        tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1
        p_accept <- min(1,1/proposal_fununction(oldTheta[2],tryTheta[2])/exp(log_target(oldTheta,r2_t_MSFT)-log_target(tryTheta,r2_t_MSFT)))
        #p_accept <- min(1,1/proposal_fun_old(oldTheta[2],tryTheta[2])/proposal_fun_next(oldTheta[2],tryTheta[1])
         #               /exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        
        #p_accept <- min(1,1/proposal_function_1(oldTheta[1],tryTheta[1])*exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        #p_accept <- min(1,(1/proposal_function_1(oldTheta[1],tryTheta[1]))/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))       
        #p_accept <- min(1,proposal_function_2(oldTheta[1],tryTheta[1])*exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_3(oldTheta[1],tryTheta[1])*exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_4(oldTheta[1],tryTheta[1])/exp(log_target(tryTheta,r2_t_GE)-log_target(oldTheta,r2_t_GE)))
        
        #p_accept <- min(1,proposal_function_5(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_6(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_7(oldTheta[1],1/tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        #p_accept <- min(1,proposal_function_8(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)-log_target(tryTheta,r2_t_GE)))
        
        if (accept_draw <= p_accept) {
            oldTheta[2] <- tryTheta[2]
            accept_rate[2] <- accept_rate[2] + 1        
        }
    
        targetSample_alpha1[i] <- oldTheta[2] #assign values
    
    }#end main loop GE

#########################################################################################
#########################################################################################
#########################################################################################
## Visual Inspection
#########################################################################################
#########################################################################################
#########################################################################################
##
## The plot() generates the values that the parameter took during the runtime of the chain

par(mfrow=c(2,1))
plot(targetSample_alpha0,xlab="Iterations",ylab=expression(alpha[0]),type='l', 
     main="Trace values of the chain (GE)")
plot(targetSample_alpha1,xlab="Iterations",ylab=expression(alpha[1]),type='l', 
     main="Trace values of the chain (GE)")

hist(targetSample_alpha0,40,xlab="Iterations",ylab=expression(alpha[0]), main=" ")
hist(targetSample_alpha1,40,xlab="Iterations",ylab=expression(alpha[1]), main=" ")
#par(mfrow=c(1,1))

#############################################################################
## basic summary statistics GE
#############################################################################
#targetSample_alpha0
#mean(targetSample_alpha0
cat("mean_alpha0 = ",  mean(targetSample_alpha0))
cat("var_alpha0 = ", var(targetSample_alpha0))
summary(targetSample_alpha0)
#targetSample_alpha1
#mean(targetSample_alpha1)
cat("mean_alpha1 = ", mean(targetSample_alpha1))
cat("var_alpha1 = ", var(targetSample_alpha1))
summary(targetSample_alpha1)
#
###############################################################
## End Gibbs-within-Metropolis sampler BAC 
###############################################################



##########################################################################################
##########################################################################################
##########################################################################################
## Thesis Testing
## Function Testing 
## Simulation af ARCH process, with known parameters
##########################################################################################
##########################################################################################
##########################################################################################
#
# We simulate a timeseries for an ARCH process in the following way
# Definition of ARCH process: 
# \sigma_t^2 = \alpha_0 + \alpha_1*\r^2_{t-1}
# r_t = \sigma_t*\epsilon_t
#
########################################################################################################
# Generate time series, part 1
########################################################################################################
#
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
plot.ts(r_t)
hist(r_t)
r2_t <- r_t^2
hist(r2_t)
mean(r2_t)
str(r2_t)
garch(r_t, 0:1)
garch(r2_t, 0:1)

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
    dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
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
    p_accept <- min(1,1/proposal_function(oldTheta[1],tryTheta[1])
                    /exp(log_target(tryTheta,r2_t)-log_target(oldTheta,r2_t)))
    
    #print(p_accept)
    if (accept_draw <= p_accept) { # accept new move
        oldTheta[1] <- tryTheta[1] # update new theta value
        accept_rate[1] <- accept_rate[1] + 1 # count good ones    
    }
    
    targetSample_alpha0[i] <- oldTheta[1] #assign new value of accepted move
    
    tryTheta <- oldTheta
    
    ####################################################################################
    
    accept_draw <- runif(1,0,1)
    tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1
    p_accept <- min(1,1/proposal_function(oldTheta[2],tryTheta[2])/exp(log_target(oldTheta,r2_t)-log_target(tryTheta,r2_t)))
    
    
    if (accept_draw <= p_accept) {
        oldTheta[2] <- tryTheta[2]
        #accept_rate[2] <- accept_rate[2] + 1        
    }
    
    targetSample_alpha1[i] <- oldTheta[2] #assign values

    tryTheta <- oldTheta
    
}#end main loop 

par(mfrow=c(2,1))
plot(targetSample_alpha0,xlab="Iterations",ylab=expression(alpha[0]),type='l', main="Fubini")
plot(targetSample_alpha1,xlab="Iterations",ylab=expression(alpha[1]),type='l', main="Fibonacci")

par(mfrow=c(2,1))
hist(targetSample_alpha0,40,xlab="Iterations",ylab=expression(alpha[0]), main="Fubini")
hist(targetSample_alpha1,40,xlab="Iterations",ylab=expression(alpha[1]), main="Fibonnaci")

summary(targetSample_alpha0)

# done

