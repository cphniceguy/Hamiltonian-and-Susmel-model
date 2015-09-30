# Thesis @ Copenhagen Business School 2015     
# The MCMC_MAGIC_HACK Metropolis-within-Gibbs Algorithm for ARCH processes
# Purpose of this program is to simulate a Markov chain from a time series, 
# by using Metropolis-within-Gibbs Algorithm.
# DISCLAMER OF WARRANTY: The author of this software provide "as is" with 
# no warranties of any kind. Use at own risk!..
# Everyone is permitted to copy, change, modify and distribute this code 
# in any possible and impossible way. 

/*Thesis in progrress,.. 
Copenhagen Business School @ 2015     
Cand.merc.(mat.) Serge Boris Nesterenko
Subject "Financial Engineering of a stock volatility for Bank of Amerika (BAC) by using Markov 
    Regime Switching Models (MS-ARCH(1)) or Hidden Markov Models (HMM) with ARCH(1) Component Regime     
    Detection"
*/

rm(list=ls()) # clear global environment
setwd("C:/Users/computer/Desktop/thesis") # set the working directory 

################################################################################
## Data manipulation 
################################################################################
# time series manipulation
library(tseries); #library(fGarch)
str(BAC)
################################################################################
## Download: Financial time series, from finance.yahoo.com
################################################################################
#
BAC  <- get.hist.quote('BAC', start='2000-01-01', end='2014-11-04', quote='AdjClose')
#GE   <- get.hist.quote('GE', start='2000-01-01', end='2014-11-04', quote='AdjClose')
#par(mfrow=c(1,2))
#plot(BAC)
#plot(GE)
#length(BAC)/2 # Note, we divide by 2, because we have 2 columns (i.e., R specific)
#length(GE)/2
#save(BAC, file="BAC.RData")       
#save(GE, file="GE.RData")
load(file="BAC.RData")
load(file="GE.RData")

# Log_Excess returns for our time series (i.e., continues excess returns)
#
BAC$Return  <- diff(log(BAC$AdjClose))
BAC$Return*100
length(BAC$Return)
BAC$Return # equivalent to {r_t}_{t=1}^T
GE$Return   <- diff(log(GE$AdjClose))
GE$Return  # equivalent to {r_t}_{t=1}^T
summary(BAC); summary(GE)
head(BAC)
class(BAC)
str(BAC)
head(BAC)
# Fubini
#
################################################################################
## Visual Inspections
################################################################################
# par(mfrow=c(2,2))
# the histograms
hist(BAC$Return, breaks=50, main="BAC")
hist(GE$Return, breaks=50, main="GE")
# Q-Q plots
qqnorm(BAC$Return)
qqline(BAC$Return)
qqnorm(GE$Return)
qqline(GE$Return)
plot(BAC, main="BAC")
legend(x="topleft", legend="BAC", lty=1:2) 
plot(GE, main="GE")
legend(x="topleft", legend="GE", lty=1:2) 
BACGE <- cbind(BAC, GE)
plot(BACGE, plot.type="single", ylab="Adjusted close price", col=c("blue", "red"), lty=1:2) 
legend(x="topleft", legend=c("BAC","GE"), col=c("blue","red"), lty=1:2) 
#
################################################################################
## Test for autocorrelation, using two tests, Ljung-Box and Box-Ljung test 
################################################################################
#
# Note that ARCH(m) processes per definition, have no autocorrelation, but posseses 
# variability in variance (cond_var). Hence when testing two null hypotehsises, which 
# (i)  test for no autocorrelation, H_0: ACF == 0
# (ii) test for no ARCH effect,     H_0: ACF^2 == 0
#acf(as.numeric(BAC$Return[c(2:1867)]))   
acf(as.numeric(BAC$Return[c(2:1867)], lag=8))   # serial correlation test  
#acf(as.numeric(BAC$Return[c(2:1867)], lag=12))  
#acf(as.numeric(BAC$Return[c(2:1867)])^2) 
acf(as.numeric(BAC$Return[c(2:1867)], lag=8)^2) # ARCH Effect Test
#acf(as.numeric(GE$Return[c(2:1867)]))
acf(as.numeric(GE$Return[c(2:1867)], lag=8))
#acf(as.numeric(GE$Return[c(2:1867)])^2)  
acf(as.numeric(GE$Return[c(2:1867)], lag=8)^2)


#################################################################################################
## Box Test for ARCH Effect (i.e., conditional heteroscedasticity)
#################################################################################################
#
#Box.test(as.numeric(BAC$Return[c(2:1867)]), type="Ljung")         
Box.test(as.numeric(BAC$Return[c(2:1867)]), lag=8, type="Ljung")# testing for a serial correlation
#Box.test(as.numeric(BAC$Return[c(2:1867)]), lag=12, type="Ljung") 
BAC$atReturn <- BAC$Return - mean(BAC$Return, na.rm=T)
mean(BAC$atReturn)

#Box.test(as.numeric(BAC$atReturn[c(2:1365)])^2, type="Ljung") 
Box.test(as.numeric(BAC$atReturn[c(2:1365)])^2, lag=8, type="Ljung") # testing for a conditional heteroscedasticity
#Box.test(as.numeric(BAC$atReturn[c(2:1365)])^2, lag=12, type="Ljung")
summary(BAC)
#Box.test(as.numeric(GE$Return[c(2:1867)]), type="Ljung")  
Box.test(as.numeric(GE$Return[c(2:1867)]), lag=8, type="Ljung") 
#Box.test(as.numeric(GE$Return[c(2:1867)]), lag=12, type="Ljung")  
GE$atReturn <- GE$Return - mean(GE$Return, na.rm=T)             
Box.test(as.numeric(GE$atReturn[c(2:1867)])^2, lag=8, type="Ljung")
summary(GE)
#Box.test(as.numeric(GE$atReturn[c(2:1867)])^2, type="Ljung")
#Box.test(as.numeric(GE$atReturn[c(2:1867)])^2, lag=5, type="Ljung")
Box.test(as.numeric(GE$atReturn[c(2:1867)])^2, lag=8, type="Ljung")


##################################################################################
##################################################################################
##################################################################################
## 
## Markov Chain Monte Carlo (MCMC) 
##  
##################################################################################
##################################################################################
##################################################################################
##################################################################################
#
# library(coda); library(lattice); library(MASS); library(MCMCpack); library(mvtnorm) 
#set.seed(55555)  # initialize Random Number Generator
#library(gtools) # debugging tool, invalide() command
#list1 <- c(1:10); list1
#list1[-length(list1)] # hack_trick to remove last element in a list
#list1[-5:-length(list1)]
#list1[-1:-3]
#
##########################################################################################
## Metropolis-within-Gibbs sampler (MWH)
##########################################################################################
# 
# The Posteriori Distibution (i.e., the loglikelihood function; the unscaled target density)
#
log_posterior <- function(alpha0,beta0,alpha1,beta1,r2_t) {
    sigma2_t <- alpha0+alpha1*r2_t[-length(r2_t)]
    sum(-0.5*log(sigma2_t)-r2_t[-1]/(2*sigma2_t))-alpha0/beta0-alpha1/beta1
} 
# input: theta parameter is 2 dim vector, contains dynamic alpha hyperparameters
log_target <- function(theta, r2_t) { # theta, r2_t
    beta0 <- 0.1  # a constant
    beta1 <- 0.2 # a constant
    log_posterior(theta[1], beta0, theta[2], beta1, r2_t) 
}

# our proposal function
#
q_prob <- function(lastTheta, nextTheta) {
    1/dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
}

q_prob1 <- function(lastTheta, nextTheta) {
    #version_1
    dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta) 
}

q_prob2 <- function(lastTheta, nextTheta) {
    #version_1
    #dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
    #version_2
    #dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta)
    #version_3
    #1/(dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta))
    1/dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
}

q_prob3 <- function(lastTheta, nextTheta) {
    #version_1
    #dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
    #version_2
    #dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta)
    #version_3
    #1/(dexp(lastTheta,1/nextTheta)/dexp(nextTheta,1/lastTheta))
    1/dexp(nextTheta,1/lastTheta)/dexp(lastTheta,1/nextTheta)
}



# Calc r^2_t vector from above mentioned returns, result is the squared time series vector
r2_t_BAC <- na.omit(BAC$atReturn^2) # na.omit() removes NA obs in a timeSeries

# Prior Values definition
#mu_alpha0_BAC <- mean(r2_t_BAC); mu_alpha0_BAC 
#mu_alpha1_BAC <- mean(r2_t_BAC); mu_alpha1_BAC
mu_alpha0_BAC <- 0.25
mu_alpha1_BAC <- 0.5
# allocation of memory space to save MWH sampling draws
# 
numSteps <- 1500  # replications
burnin <- 250


########################################################
########################################################
## Begin Metropolis-within-Gibbs (MwG) Algorithm
########################################################
########################################################
#
##########################################################################################
## 1) MwG - BAC Bank of America 
##########################################################################################
# 
oldTheta <- c(mu_alpha0_BAC,mu_alpha1_BAC) # a position in a previous state
oldTheta
tryTheta <- oldTheta; tryTheta

for (i in 2:burnin) { # the burn-in phase

    accept_draw <- runif(1,0,1) # Generate an innovation
    tryTheta[1] <- oldTheta[1]*rexp(1,1) # generate a draw for alpha0
    p_accept <- min(1,q_prob(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_BAC)-log_target(tryTheta,r2_t_BAC)))
#    p_accept <- min(1,q_prob(oldTheta[1],tryTheta[1])*exp(log_target(tryTheta,r2_t_BAC)-log_target(oldTheta,r2_t_BAC)))
#    print(p_accept)
    if (accept_draw <= p_accept) { # accept new move
        oldTheta[1] <- tryTheta[1] # update new theta value
    }
    
    accept_draw <- runif(1,0,1)
    tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1
    
    p_accept <- min(1,q_prob(oldTheta[2],tryTheta[2])/exp(log_target(oldTheta,r2_t_BAC)-log_target(tryTheta,r2_t_BAC)))
#    p_accept <- min(1,q_prob(oldTheta[2],tryTheta[2])*exp(log_target(tryTheta,r2_t_BAC)-log_target(oldTheta,r2_t_BAC)))
    
    if (accept_draw <= p_accept) {
        oldTheta[2] <- tryTheta[2]
    }
        
} # end burnin loop BAC
#
###############################################################################
# main loop BAC, begin to collect info about draws from the log_target sample
###############################################################################
#
accept_rate <- c(0,0) # rate of acceptance of successive moves
targetSample_alpha0 <- numeric(numSteps); str(targetSample_alpha0)
targetSample_alpha1 <- numeric(numSteps)
#
for (i in burnin+1:numSteps){

    accept_draw <- runif(1,0,1) # Generate an innovation
    tryTheta[1] <- oldTheta[1]*rexp(1,1) # generate a draw for alpha0
    p_accept <- min(1,q_prob(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_BAC)-log_target(tryTheta,r2_t_BAC)))
#   p_accept <- min(1,q_prob(oldThetaa[1],tryTheta[1])*exp(log_target(tryTheta,r2_t_BAC)-log_target(oldTheta,r2_t_BAC)))
#    print(p_accept)
    if (accept_draw <= p_accept) { # accept new move
        oldTheta[1] <- tryTheta[1] # update new theta value
        accept_rate[1] <- accept_rate[1] + 1 # count good ones    
    }
    
    targetSample_alpha0[i-burnin] <- oldTheta[1] #assign new value of accepted move
            
    accept_draw <- runif(1,0,1)
    tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1

    p_accept <- min(1,q_prob(oldTheta[2],tryTheta[2])/exp(log_target(oldTheta,r2_t_BAC)-log_target(tryTheta,r2_t_BAC)))
    
    if (accept_draw <= p_accept) {
        oldTheta[2] <- tryTheta[2]
        accept_rate[2] <- accept_rate[2] + 1        
    }
        
    targetSample_alpha1[i-burnin] <- oldTheta[2] #assign values

}
#
## end main loop BAC
#


###############################################################################
## "BAC" Visual inspections 
###############################################################################

par(mfrow=c(2,1))
#par(mfrow=c(2,2))
plot(targetSample_alpha0,xlab="Iterations",ylab=expression(alpha[0]),type='l', 
     main="Trace values of the c hain (BAC)")
plot(targetSample_alpha1,xlab="Iterations",ylab=expression(alpha[1]),type='l', 
     main="Trace values of the chain (BAC)")

hist(targetSample_alpha0,40,xlab="Iterations",ylab=expression(alpha[0]),main=" ")
hist(targetSample_alpha1,40,xlab="Iterations",ylab=expression(alpha[1]),main=" ")
#par(mfrow=c(1,1))

#############################################################################
## basic summary statistics BAC
#############################################################################
#
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


##########################################################################################
##########################################################################################
##########################################################################################
## 2) MwG - GE General Electric 
##########################################################################################
##########################################################################################
##########################################################################################
#
numSteps <- 10000  # replications
burnin <- 1000

# Prior Values definition
r2_t_GE <-  na.omit(GE$atReturn^2); r2_t_GE 
#mu_alpha0_GE <- mean(r2_t_GE); mu_alpha0_GE 
#mu_alpha1_GE <- mean(r2_t_GE); mu_alpha1_GE
mu_alpha0_GE <- 0.00025
mu_alpha1_GE <- 0.5

oldTheta <- c(mu_alpha0_GE,mu_alpha1_GE) # a position in a previous state
oldTheta
tryTheta <- oldTheta; tryTheta

# burn-in phase GE
for (i in 2:burnin) { 
    
    accept_draw <- runif(1,0,1) # Generate an innovation
    tryTheta[1] <- oldTheta[1]*rexp(1,1) # generate a draw for alpha0
    p_accept <- min(1,1/q_prob(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)
                                                            -log_target(tryTheta,r2_t_GE)))
    print(p_accept)
    if (accept_draw <= p_accept) { # accept new move
        oldTheta[1] <- tryTheta[1] # update new theta value
    }
    
    accept_draw <- runif(1,0,1)
    tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1
    p_accept <- min(1,1/q_prob(oldTheta[2],tryTheta[2])/exp(log_target(oldTheta,r2_t_GE)
                                                            -log_target(tryTheta,r2_t_GE)))
    
    if (accept_draw <= p_accept) {
        oldTheta[2] <- tryTheta[2]
    }
    
}#end burnin loop GE

# main loop GE, begin to collect info about draws from the log_target sample
#
accept_rate <- c(0,0) # rate of acceptance of successive moves
targetSample_alpha0 <- numeric(numSteps); str(targetSample_alpha0)
targetSample_alpha1 <- numeric(numSteps)
#
for (i in burnin+1:numSteps){
    
    accept_draw <- runif(1,0,1) # Generate an innovation
    tryTheta[1] <- oldTheta[1]*rexp(1,1) # generate a draw for alpha0
    #   p_accept <- min(1,q_prob(oldTheta[1],tryTheta[1])*exp(log_target(tryTheta,r2_t_BAC)-log_target(oldTheta,r2_t_BAC)))
    p_accept <- min(1,1/q_prob(oldTheta[1],tryTheta[1])/exp(log_target(oldTheta,r2_t_GE)
                                                            -log_target(tryTheta,r2_t_GE)))
    print(p_accept)
    if (accept_draw <= p_accept) { # accept new move
        oldTheta[1] <- tryTheta[1] # update new theta value
        accept_rate[1] <- accept_rate[1] + 1 # count good ones    
    }
    
    targetSample_alpha0[i-burnin] <- oldTheta[1] #assign new value of accepted move
    
    accept_draw <- runif(1,0,1)
    tryTheta[2] <- oldTheta[2]*rexp(1,1) # generate a draw for alpha1
    p_accept <- min(1,1/q_prob(oldTheta[2],tryTheta[2])/exp(log_target(oldTheta,r2_t_GE)
                                                            -log_target(tryTheta,r2_t_GE)))
    
    if (accept_draw <= p_accept) {
        oldTheta[2] <- tryTheta[2]
        accept_rate[2] <- accept_rate[2] + 1        
    }
    
    targetSample_alpha1[i-burnin] <- oldTheta[2] #assign values
    
}#end main loop GE

# The plot() generates the values that the parameter took during the runtime of the chain

par(mfrow=c(3,2))
#par(mfrow=c(2,2))

plot(targetSample_alpha0,xlab="Iterations",ylab=expression(alpha[0]),type='l', 
     main="Trace values of the chain (GE)")
plot(targetSample_alpha1,xlab="Iterations",ylab=expression(alpha[1]),type='l', 
     main="Trace values of the chain (GE)")

hist(targetSample_alpha0,40,xlab="Iterations",ylab=expression(alpha[0]), main=" ")
hist(targetSample_alpha1,40,xlab="Iterations",ylab=expression(alpha[1]), main=" ")

plot(targetSample_alpha0,col=1:100000)
plot(targetSample_alpha1,col=1:100000)

#par(mfrow=c(1,1))

#############################################################################
## basic summary statistics GE
#############################################################################
#
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






###############################################################
###############################################################
###############################################################
## 
## Markov-Switching (MS) 
## 
###############################################################
###############################################################
###############################################################
#
#library(MSwM)

#garch(as.numeric(na.omit(GE$atReturn)))































































# Serge Boris was here!..