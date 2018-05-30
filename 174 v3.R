# packages
library(MASS)
library(qpcR)
library(stats)
library(survMisc)
library(ggplot2)
library(ggpubr)

# read in data
getwd()
emissiondata<-read.csv("MER_T12_06.csv",header=TRUE)
TXEIEUS<-read.csv("TXEIEUS.csv",header=TRUE)
# time series data
Total_Energy_Electric_Power <- ts(TXEIEUS[,3],start=c(1996,1),frequency=12) 

# time series
ts=ts(TXEIEUS$Value,start=c(1996,1),frequency=12)
ts=tsclean(ts)
t = 1:length(ts)
fit = lm(ts ~ t)

# transformation
bcTransform = boxcox(ts ~ t,plotit = TRUE)
lambda <- bcTransform$x[which(bcTransform$y == max(bcTransform$y))]
# close to 0 so we'll use log of the data
ts.log<-log(Total_Energy_Electric_Power)
plot(ts, main = "Log-Transformed Time Series")

# decomposing
ggtsdisplay(ts,lag.max=120, main = "Time Series with ACF and PACF")
Decomp=decompose(ts) # survMisc package
autoplot(Decomp,main="Decomposed Time Series") + theme(axis.text.y=element_text(size=6),text=element_text(size=10))+ xlab("Year") #forecast package

# detrending the data
dif2 = diff (ts , lag = 2)
ggtsdisplay ( dif2 , lag.max = 60, main = "De-Trended Data")
autoplot(decompose(dif2),main="De-Trended Data") +
  theme(axis.text.y=element_text(size=6),text=element_text(size=10))+
  xlab("Year")

# de-seasoning the data
dif12 = diff ( dif2 , lag = 12)
ggtsdisplay ( dif12 , lag.max = 60, main = "De-trended and De-seasonalized Data")
autoplot(decompose(dif12),main= "De-trended and De-seasonalized Data") +
  theme(axis.text.y=element_text(size=6),text=element_text(size=10))+
  xlab("Year")

# generate possible models to use with a matrix of the AICc values to compare
Fit = list ()
AICc=matrix (, nrow =200 , ncol =7)
colnames ( AICc ) =c (" p " ," q " ," P " ," Q " ," AICc " ," Shapiro " ," LjungBox ")
i =0
for (P in c (0:3) ) {
  for (Q in c (0:3) ) {
    for (p in c (0:3) ){
      for (q in c (0:3) ){
        Fit[[i+1]]=sarima (ts ,p ,2 ,q ,P ,1 ,Q ,12 , Model = TRUE , details = FALSE )$fit
        plot.new()
        AICc[i+1,1]=p
        AICc[i+1,2]=q
        AICc[i+1,3]=P
        AICc[i+1,4]=Q
        AICc[i+1,5]=sarima(Total_Energy_Electric_Power,p,1,q,P,1,Q,12,Model=TRUE,details=FALSE)$AICc
        plot.new()
        AICc[i+1,6]=shapiro.test(resid(Fit[[i+1]]))$p.value
        AICc[i+1,7]=Box.test(resid(Fit[[i+1]]),type=c("Ljung-Box"),lag=12)$p.value
        i=i+1
      }
    }
  }
}

# gives error and 12-13 warnings
 AICc<-data.frame(AICc)
# sort by AICc ( Increasing order)
AICc_sorted<-AICc[order(AICc$X.AICc.),]
# Filter the models by Shapiro test
AICc_cmp<-subset(AICc_sorted,AICc_sorted$X.Shapiro.>0.05)

# Filter the models by LjungBox test (for checking independence)
AICc_cmp<-subset(AICc_cmp,AICc_cmp$X.LjungBox.>0.05)

# final models to consider 
AICc_cmp

model1 = sarima(ts ,3 ,2 ,3 ,0 ,1 ,1 ,12) #astsa package
model2 = sarima(ts ,3 ,2 ,3 ,0 ,1 ,2 ,12)

# diagnostic checking

# 1. normality
par ( mfrow = c (1 ,2) )
resid1 = residuals(model1$fit )
resid2 = residuals(model2$fit )
qqnorm (resid1,main ="Model 1 Normal Q - Q Plot")
qqline(resid1,col ="blue")
qqnorm (resid2,main ="Model 2 Normal Q - Q Plot")
qqline(resid2,col ="blue")
par(mfrow=c(1,2))
hist(resid1, main = "Histogram of Residuals from Model 1")
hist(resid2, main = "Histogram of Residuals from Model 2")

# 2. residual plots (when they fall within the confidence intervals
# then we know that they have constant variance. They almost all fall
# within the boundaries so we can assume this.)
ggtsdisplay ( resid1 , main ="Residual Plots for Model 1" , xlab =" Year ")
p1 = ggAcf ( resid1 ^2 , lag.max = 60 , main ="")
p2 = ggPacf ( resid1 ^2 , lag.max = 60 , main ="")
ggarrange ( p1 ,p2 , nrow = 2, ncol = 2) #ggpubr
ggtsdisplay ( resid2 , main ="Residual Plots for Model 2" , xlab =" Year ")
p11 = ggAcf ( resid2 ^2 , lag.max = 60 , main ="")
p22 = ggPacf ( resid2 ^2 , lag.max = 60 , main ="")
ggarrange ( p11 ,p22 , nrow = 2, ncol = 2) #ggpubr

# choosing our fit
fit=model1 #sarima (3,2,3)x(0,1,1)_12

# Forecast and transfer to original scale
library(forecast)

TXEIEUS2 <- read.csv("TXEIEUS2.csv",header=TRUE)
ts2 <- ts(TXEIEUS2$Value,start=c(1996,1),frequency=12) # we include the actual values for 2016

predict<-forecast(ts, level=c(95))
autoplot(ts2,main="Original Data") + ylab("Carbon Emission Per metric Ton actual values")+xlab("Month")+theme(legend.position="None")
autoplot(predict, main="Forecasted Values") + ylab("Carbon Emission Per metric Ton predicted values")+xlab("Month")+theme(legend.position="None")