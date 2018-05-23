getwd()
setwd("/Users/denniswang/Desktop")
emissiondata<-read.csv("/Users/denniswang/Desktop/MER_T12_06.csv",header=TRUE)

TXEIEUS<-read.csv("/Users/denniswang/Desktop/Pstat 174 files/TXEIEUS.csv",header=TRUE)
Total_Energy_Electric_Power <- ts(TXEIEUS[,3],start=c(1996,1),frequency=12) #Time series data

#Plotting original data
par(mfrow=c(1,1))
plot(Total_Energy_Electric_Power,main="Total Energy Electric Power CO2 Emission per Million Metric Ton") #Plotting time series
#Nonconstant variance
par(mfrow = c(1,2))
acf(Total_Energy_Electric_Power,lag.max = 60,main = "")
#Clearly cyclicial trends
pacf(Total_Energy_Electric_Power,lag.max = 60,main = "")
title("Original Time Series", line = -1, outer=TRUE)
var(Total_Energy_Electric_Power) #Very big variance
max(Total_Energy_Electric_Power)
min(Total_Energy_Electric_Power)
#Detrending
par(mfrow=c(1,1))
t = 1:length(Total_Energy_Electric_Power)
fit = lm(Total_Energy_Electric_Power ~ t)
bcTransform = boxcox(Total_Energy_Electric_Power ~ t,plotit = TRUE) #package qpcR
#Lab03 since the CI includes 0, then BC transformation is lamda or log(data)
lambda = bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
Total_Energy_Electric_Power.tr <- (1/lambda)*(Total_Energy_Electric_Power^lambda - 1)
Total_Energy_Electric_Power.bc<-log(Total_Energy_Electric_Power)
op <- par(mfrow = c(1,3))
ts.plot(Total_Energy_Electric_Power,main = "Original data",ylab = expression(X[t]))
ts.plot(Total_Energy_Electric_Power.bc,main = "Box-Cox tranformed data-log", ylab = expression(Y[t]))
ts.plot(Total_Energy_Electric_Power.tr,main = "Box-Cox tranformed data-lambda", ylab = expression(Y[t]))
#As seen, the plots do not differ much from each other
#B.C transformed ACF and PACF
par(mfrow=c(1,2))
acf(Total_Energy_Electric_Power.bc,lag.max = 60,main = "")
pacf(Total_Energy_Electric_Power.bc,lag.max = 60,main = "")
title("Box-Cox Transformed Time Series", line = -1, outer=TRUE)
#It looks like significant correlation at every 12 lag and we can conclude seasonal component d=12
#Comparing variance between the original and transformed
var(Total_Energy_Electric_Power)
var(Total_Energy_Electric_Power.bc) #Way much lower - log trans
var(Total_Energy_Electric_Power.tr) #The lowest variance - lambda transformation

#Differencing once to remove trend
par(mfrow=c(1,1))
y1 = diff(Total_Energy_Electric_Power.tr, 2)
plot(y1,main = "De-trended Time Series at lag 2",ylab = expression(nabla~Y[t])) #quad
abline(h = 0,lty = 2)
acf(y1, main="")
pacf(y1, main="")
var(y1) #lower than B.C. no transformation
var(diff(y1,1)) #Differencing more than once actually increase the variance, indicating over differencing 
#Differencing more than once at lag of 1 increased variance
#Differencinng more than once at lag of 2 actually lowered variance

# Diference at lag = 12 (cycle determined by the ACF) to remove seasonal component
y12 = diff(y1, 12)
ts.plot(y12,main = "De-trended/seasonalized Time Series of 12",ylab = expression(nabla^{12}~nabla~Y[t]))
abline(h = 0,lty = 2)
var(y12)
# Re-calculate the sample variance and examine the ACF and PACF for first difference
op = par(mfrow = c(1,2))
acf(y1,lag.max = 60,main = "")
pacf(y1,lag.max = 60,main = "") #P is possibly 1
title("De-trended Time Series", line = -1, outer=TRUE)

# Re-calculate the sample variance and examine the ACF and PACF for 12th
op = par(mfrow = c(1,2))
acf(y12,lag.max = 60,main = "")
pacf(y12,lag.max = 60,main = "")
title("De-trended/seasonalized Time Series",line = -1, outer=TRUE)
#Possibly MA(1) from ACf as the rest are in 95 CI but PACF does not look like it
ar(y12, method="yule-walker")
ar(y12, aic = TRUE, order.max = NULL, method = c("ols")) #using OLS method
ar(y12, aic = TRUE, order.max = NULL, method = c("mle")) #using MLE method
#All three method suggest order of 12, so AR(12)
