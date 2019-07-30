library(TSA)

##DATA SETUP##
indiaData <- read.table("C:\\Users\\Nick\\Desktop\\Urbanization Research_CHINA_TS\\India_Urbanization_Data(noYears).txt" , skip=1)
names(indiaData) <- c("urbanRate_India" , "totalPopulation_India")

totalPopulation_India <- indiaData[,1]
urbanRate_India <- indiaData[,2]
length(urbanRate_India)
##DATA INTERPOLATION##
library(imputeTS)
urbanRate_India.I <- na.interpolation(urbanRate_India)
indiaData.I <- na.interpolation(indiaData)

##TIME SERIES CREATION##
urbanRate_India.ts <- ts(urbanRate_India.I , start = 1961 , end = 2011)

##GRAPHICAL REPRESENTATION##
plot(urbanRate_India.ts)

##Stationarity Check && Model Selection##
acf(urbanRate_India.ts)
pacf(urbanRate_India.ts)

eacf(urbanRate_India.ts)

library(fUnitRoots)

adfTest(urbanRate_India.ts , lags = 1 , type = c("ct"))
urbanRate_India.ts.diff <- diff(urbanRate_India.ts)
pacf(urbanRate_India.ts.diff)

adfTest(urbanRate_India.ts.diff , lags = 2 , type = c("c"))
plot(urbanRate_India.ts.diff)

acf(urbanRate_India.ts.diff)
eacf(urbanRate_India.ts , ar.max = 11 , ma.max = 11)
plot(armasubsets(urbanRate_India.ts.diff , nar = 10 , nma = 10 , y.name='test',ar.method='ols'))


##TS MODEL CLEANING -- [FAILURE]## --Note: This code is only hear for experimental purposes. At present
library(forecast)
urbanRate_India.ts.clean <- tsclean(urbanRate_India.ts , replace.missing = TRUE)
par(mfrow=c(1,2))
plot(urbanRate_India.ts.clean , ylab = 'Cleaned Urbanization Rate')
plot(urbanRate_India.ts , ylab = 'Standard Urbanization Rate')
library(forecast)
auto.arima(urbanRate_India.ts.clean)
urbanRate_India.ari.clean <- arima(urbanRate_India.ts.clean , order = c(1,1,0) , method = "ML")

urbanRate_India.ts.clean.diff <- diff(diff(urbanRate_India.ts.clean))
par(mfrow=c(1,1))
plot(urbanRate_India.ts.clean.diff)
adfTest(diff(urbanRate_India.ts.clean.diff) , lags = 5, type = c("c"))

bc <- BoxCox.ar(urbanRate_India.ts , lambda = seq(-5,5) , method = "ols")
testModel <- auto.arima(1/(urbanRate_India.ts)^4)

##Model Specification -- [UNCLEAN]##
acf(urbanRate_India.ts.diff)
pacf(urbanRate_India.ts.diff)
eacf(urbanRate_India.ts.diff , ar.max = 11)

par(mfrow=c(1,1))
plot(armasubsets(urbanRate_India.ts.diff , nar = 12 , nma = 12))

par(mfrow=c(2,1))
acf(urbanRate_India.ts.diff)
pacf(urbanRate_India.ts.diff)

##MODEL SPECIFICATION -- [CLEAN]##
acf(urbanRate_India.ts.clean.diff)
pacf(urbanRate_India.ts.clean.diff)
eacf(urbanRate_India.ts.clean.diff)

plot(armasubsets(urbanRate_India.ts.diff , nar = 10 , nma = 10))
auto.arima(urbanRate_India.ts.clean.diff)

urbanRate_India.ari.clean.log <- arima(log(urbanRate_India.ts.clean) , order=c(10,2,0))

#Will create two models to analyze. ARIMA(10,1,10) and IMA(1,10) for [UNCLEAN]#
urbanRate_India.arima <- arima(urbanRate_India.ts , order = c(10,1,10))
urbanRate_India.ima <- arima(urbanRate_India.ts , order = c(0,1,10))

##ANALYSIS OF RESIDUALS -- [UNCLEAN]##
plot(residuals(urbanRate_India.arima))
plot(residuals(urbanRate_India.ima))
tsdiag(urbanRate_India.arima)
tsdiag(urbanRate_India.ima)

runs(residuals(urbanRate_India.ima))

#Creating models to analyze for [CLEAN]#
urbanRate_India.ari.clean <- arima(urbanRate_India.ts.clean.diff , order = c(10,2,0))
urbanRate_India.arima.clean.auto <- auto.arima(urbanRate_India.ts.clean.diff , trace = TRUE)

##ANALYSIS OF RESIDUALS -- [CLEAN]##
plot(residuals(urbanRate_India.ari.clean))
plot(residuals(urbanRate_India.arima.clean.auto))

tsdiag(urbanRate_India.ari.clean , gof.lag=20)
tsdiag(urbanRate_India.arima.clean.auto)

##OUTLIER DETECTION -- [UNCLEAN]##
#31 corresponds to 1991 observation. 21 to 1981
library(tsoutliers)
detectAO(urbanRate_India.ima)
detectIO(urbanRate_India.ima)

locate.outliers(residuals(urbanRate_India.ima), coefs2poly(urbanRate_India.ima, add = TRUE)
                , cval = 3.5, types = c("AO" , "IO"),
                delta = 0.7)

#Marking Most significant Outliers
urbanRate_India.ima.detOut <- arima(urbanRate_India.ts,order=c(0,1,10),io=c(21,31))
detectIO(urbanRate_India.ima.detOut)


##Residual Plots##
plot(residuals(urbanRate_India.ima))

##NORMALITY CHECK -- [UNCLEAN]##
urbanRate_India.ima.log <- arima(log(urbanRate_India.ts), order = c(0,1,10))
par(mfrow=c(1,2))
qqnorm(residuals(urbanRate_India.ima))
qqline(residuals(urbanRate_India.ima.log))

library(normtest)
jb.norm.test(residuals(urbanRate_India.ima.log), nrepl=2000)

shapiro.test(residuals(urbanRate_India.ima.log))

##NORMALITY CHECK -- [CLEAN]##
qqnorm(residuals(urbanRate_India.ari.clean))
qqline(residuals(urbanRate_India.ari.clean))

qqnorm(residuals(urbanRate_India.arima.clean.auto))
qqline(residuals(urbanRate_India.arima.clean.auto))

qqnorm(residuals(urbanRate_India.ari.clean.log))
qqline(residuals(urbanRate_India.ari.clean.log))

##TRANSFORMATION##
bc <- BoxCox.ar(urbanRate_India.ts.clean , lambda = seq(-5,5) , method = "ols")
normalityTest <- 1/(urbanRate_India.ts^2)
urbanRate_India.ima.normal <- arima(normalityTest , order = c(0,1,10))

shapiro.test(residuals(urbanRate_India.ima.normal))

#Comparing pre transformation and post-transformation QQ-Plots#
par(mfrow=c(2,1))
qqnorm(residuals(urbanRate_India.ima) , ylab = "Untransformed Quantiles")
qqline(residuals(urbanRate_India.ima))
qqnorm(residuals(urbanRate_India.ima.normal) , ylab = "Transformed Quantiles")
qqline(residuals(urbanRate_India.ima.normal))

##COEF SIGNIFICANCE TEST##
library(lmtest)
urbanRate_India.ima.detOut <- arima(urbanRate_India.ts,order=c(0,1,10) , method='ML' , io=c(21 , 31))
urbanRate_India.ima2 <- arima(urbanRate_India.ts , order = c(0,1,1) , method = 'ML')

coeftest(urbanRate_India.ima.detOut)
coeftest(urbanRate_India.ima)

##PREDICT -- [UNCLEAN]##
rate<-predict(urbanRate_India.ima2,n.ahead=1)

par(mfrow=c(1,1))

plot(urbanRate_India.ima.detOut , n.head=1 , type='b' , xlab='Year' , ylab='Urban Rate')
plot(urbanRate_India.ima2,n.ahead=1,type='b', col="red",xlab='Time')

library(forecast)
library(tseries)
future <- forecast.ar(urbanRate_India.arima, h = 10)
par(mfrow=c(1,2))
plot(future2)
plot(forecast(auto.arima(urbanRate_India.ts)))

plot(urbanRate_India.ima.detOut , type = "ma")

##PREDICT -- [CLEAN]##

#1. Predicting second difference.
pred <- predict(urbanRate_India.ari.clean, n.ahead = 10)
pred
plot(urbanRate_India.ts.clean.diff, type= 'l' , xlim=c(1961,2021),ylim=c(-.01,.01),xlab = 'Year' ,ylab = 'Urban Rate')
lines(pred$pred , col="blue")
lines(pred$pred+2*pred$se,col="orange")
lines(pred$pred-2*pred$se,col="orange")

future2 <- forecast(auto.arima(urbanRate_India.ts.clean))
plot(future2)

##DATA SPLIT##
urbanRate_India.ts.train <- ts(urbanRate_India.I , start = 1961 , end = 2001)
urbanRate_India.ts.test <- ts(urbanRate_India.I , start= 2002 , end = 2011)

pred2002 <- predict(urbanRate_India.ari.clean , n.ahead = 1)


##ONE-STEP FORECASTING##
pred1<-rep(NA,10)

for(i in 1:10){
  urbanRate_India.ts.train<-urbanRate_India.ari.clean[1:(2002+i-1)]
  
  arimaModel.forecasted<-arima(urbanRate_India.ts, order=c(0,0,10),method='ML')
  
  pred1<-predict(arimaModel.forecasted,n.ahead=10)$pred
  
  pred1_1<-pred1[1]
  pred1_3<-pred1[3]
  
}

y_pred<-pred1
