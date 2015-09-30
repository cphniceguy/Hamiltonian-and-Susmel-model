library(tseries)
library(fGarch)

setwd("C:/Users/computer/Desktop/thesis") # set the working directory 

load(file="BAC.RData")

BAC$Return  <- diff(log(BAC$AdjClose))
data <- na.omit(as.vector(BAC$Return))
gf<-garchFit(~garch(1, 1), data)
gf
plot(gf)
1
2
3
