knitr :: opts_chunk$set ( echo = TRUE )
install.packages(c("forecast","astsa" ,"TSA" ," GeneCycle " ," zoo " ," MASS " ," astsa " ," ggfortify " ," tseries " ," forecast " ,"
tidyverse " ," ggplot2 ", " knitr ", " readr " ," ggpubr " ," ggridges " ," tibble " ," stringr " ," tidyr " ,"
                    anytime " ," pracma " ," hydrostats ")) 
#only do it once
Packages = c ("forecast","astsa" ,"TSA" ," GeneCycle " ," zoo " ," MASS " ," astsa " ," ggfortify " ," tseries " ," forecast " ,"
tidyverse " ," ggplot2 ", " knitr ", " readr " ," ggpubr " ," ggridges " ," tibble " ," stringr " ," tidyr " ,"
anytime " ," pracma " ," hydrostats ")
library(MASS)
library(qpcR)
library(stats)
library(survMisc)
library(ggplot2)
library(ggpubr)
invisible(lapply(Packages,library,character.only = TRUE ))
knitr :: opts_chunk$set ( fig.width =7 , fig.height =5)
options ( digits = 4)
opts_chunk$set (tidy.opts = list ( width.cutoff =50) , tidy = TRUE )


getwd()
emissiondata<-read.csv("MER_T12_06.csv",header=TRUE)
TXEIEUS<-read.csv("TXEIEUS.csv",header=TRUE)
Total_Energy_Electric_Power <- ts(TXEIEUS[,3],start=c(1996,1),frequency=12) #Time series data
ts=ts(TXEIEUS$Value,start=c(1996,1),frequency=12)
ts=tsclean(ts) #forecast package
t = 1:length(ts)
fit = lm(ts ~ t)
#transform
bcTransform = boxcox(ts ~ t,plotit = TRUE) #package qpcR
#Lab03 since the CI includes 0, then BC transformation is lamda or log(data)
ts.log<-log(Total_Energy_Electric_Power)
plot(ts)

ggtsdisplay(ts,lag.max=120)
Decomp=decompose(ts) #survMisc package
autoplot(Decomp,main="") + theme(axis.text.y=element_text(size=6),text=element_text(size=10))+ xlab("Year") #forecast package

#Detrend
dif2 = diff (ts , lag = 2)
ggtsdisplay ( dif2 , lag.max = 60)

autoplot(decompose(dif2),main="") +
  theme(axis.text.y=element_text(size=6),text=element_text(size=10))+
  xlab("Year")

#de - seasonality
dif12 = diff ( dif2 , lag = 12)
ggtsdisplay ( dif12 , lag.max = 60)
autoplot(decompose(dif12),main="") +
  theme(axis.text.y=element_text(size=6),text=element_text(size=10))+
  xlab("Year")

# Generate all interesting models
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
AICc<-data.frame(AICc)
# sort by AICc ( Increasing order)
AICc_sorted<-AICc[order(AICc$X.AICc.),]
# Filter the models by Shapiro test
AICc_cmp<-subset(AICc_sorted,AICc_sorted$X.Shapiro.>0.05)
# Filter the models by LjungBox test
AICc_cmp<-subset(AICc_cmp,AICc_cmp$X.LjungBox.>0.05)
# Final possible models
AICc_cmp

model1 = sarima(ts ,3 ,2 ,3 ,0 ,1 ,1 ,12) #astsa package
model2 = sarima(ts ,3 ,2 ,3 ,0 ,1 ,2 ,12)

# check invertible / stationary


# normality
par ( mfrow = c (1 ,2) )
resid1 = residuals(model1$fit )
resid2 = residuals(model2$fit )
qqnorm (resid1,main =" Normal Q - Q Plot for Model I ")
qqline(resid1,col ="blue")
qqnorm (resid2,main =" Normal Q - Q Plot for Model II ")
qqline(resid2,col ="blue")
par(mfrow=c(1,2))
hist(resid1)
hist(resid2)
# residual plots
ggtsdisplay ( resid1 , main ="" , xlab =" Year ")
p1 = ggAcf ( resid1 ^2 , lag.max = 36 , main ="")
p2 = ggPacf ( resid1 ^2 , lag.max = 36 , main ="")
ggarrange ( p1 ,p2 , nrow = 2, ncol = 2) #ggpubr
ggtsdisplay ( resid2 , main ="" , xlab =" Year ")
p11 = ggAcf ( resid2 ^2 , lag.max = 36 , main ="")
p22 = ggPacf ( resid2 ^2 , lag.max = 36 , main ="")
ggarrange ( p11 ,p22 , nrow = 2, ncol = 2) #ggpubr
# choose fit
fit=model1 #sarima (3,2,3)x(0,1,1)_12

# Forecast and transfer to original scale
library(forecast)
forecast(ts, level=c(95))

Box.test(resid1, type="Ljung")
shapiro.test(resid1)
Box.test(resid2, type="Ljung")
shapiro.test(resid2)

