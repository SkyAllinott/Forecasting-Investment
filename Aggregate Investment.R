library(cansim)
library(forecast)
library(stargazer)
library(MTS)
library(vars)
library(ggplot2)
library(AER)
options(scipen=999)

buscapform <- get_cansim_vector("v62144081", start_time=as.Date("1961-01-01"), 
                                end_time=as.Date("2004-12-31"))
govtcapform <- get_cansim_vector("v62144105", start_time=as.Date("1961-01-01"), 
                                 end_time=as.Date("2004-12-31"))
gdp <- get_cansim_vector("v62305814", start_time=as.Date("1961-01-01"), 
                         end_time=as.Date("2004-12-31"))

CPI <- get_cansim_vector("v62307283", start_time = as.Date("1961-01-01"), 
                         end_time=as.Date("2004-12-31"))
baseline <- mean(CPI$VALUE[125:128])
CPI <- (CPI$VALUE/baseline)*100

# Capital portion of CCPI92
CPI_buscapital <- get_cansim_vector("v62307268", start_time = as.Date("1961-01-01"), 
                                    end_time=as.Date("2004-12-31"))
basebus1992 <- mean(CPI_buscapital$VALUE[125:128])
CPI_buscapital$real <- (CPI_buscapital$VALUE/basebus1992)*100

CPI_govtcapital <- get_cansim_vector("v62307275", start_time = as.Date("1961-01-01"), 
                                     end_time=as.Date("2004-12-31"))
basegovt1992 <- mean(CPI_govtcapital$VALUE[125:128])
CPI_govtcapital$real <- (CPI_govtcapital$VALUE/basegovt1992)*100


data.all <- data.frame(Date=seq(as.Date("1961-01-01"), as.Date("2004-12-31"), "quarter"))
data.all$SCPI1992 <- CPI
data.all$Capital <- ((buscapform$VALUE/CPI_buscapital$real)*100)+
  (govtcapform$VALUE/CPI_govtcapital$real)*100
data.all$Investment <- as.numeric(rbind(NA, as.matrix(diff(data.all$Capital))))
data.all$GDP2 <- (gdp$VALUE/data.all$SCPI1992)*100



data.all.trimmed <- data.all[1:156,]
testdata <- cbind(data.all$Investment[157:176], data.all$Capital[157:176], data.all$GDP[157:176])
colnames(testdata) <- c("Investment", "Capital", "GDP")
xregtest <- cbind(testdata[,"Capital"], testdata[,"GDP"])
xreg <- cbind(lag(data.all.trimmed$Capital), data.all.trimmed$GDP)

forecasts <- data.frame(TrueValues=testdata[,"Investment"])

# Univariate case:
d <- 1                                # Specified non-seasonal difference
D <- 1 
n.lag.max <- 5                        # Maximum AR
n.err.max <- 5                        # Max MA
n.seaslag.max <- 3                    # Max seasonal AR
n.seaserr.max <- 3                    # Max seasonal MA


seasonal.mse <- matrix(NA, nrow=n.seaserr.max+1, ncol=n.seaslag.max+1)
ARMA.mse <- matrix(NA, nrow=n.err.max+1, ncol=n.lag.max+1)
ARMA.arr <- matrix(NA, nrow=n.err.max+1, ncol=n.lag.max+1)
final.mse <- matrix(NA, nrow=1, ncol=1)
final.ARMA <- matrix(NA, nrow=1, ncol=1)
final.seasonal <- matrix(NA, nrow=1, ncol=1)
seasonal <- matrix(NA, nrow=1, ncol=2)
for(i in 1:n.lag.max){
  for(j in 1:n.err.max){
    for(l in 1:n.seaslag.max){
      for(m in 1:n.seaserr.max){
        model <- Arima(ts(data.all.trimmed$Investment, start=c(1961,1), frequency=4), 
                       order=c((i-1),d,(j-1)), seasonal=c((l-1),D, (m-1)))
        predict <- forecast(model, h=20)
        seasonal.mse[l,m] <- mean((forecasts[,"TrueValues"] - predict$mean)^2)
        ARMA.mse[i,j] <- min(seasonal.mse, na.rm=TRUE)
        ARMA.arr[i,j] <- paste(which(seasonal.mse == min(seasonal.mse, na.rm=TRUE), 
                                     arr.ind = TRUE), collapse=', ')
        seasonal <- which(ARMA.mse == min(ARMA.mse, na.rm=TRUE), arr.ind=TRUE)
        final.mse[1,] <- round(min(ARMA.mse, na.rm=TRUE), digits = 2)
        final.ARMA[1,] <- paste(which(ARMA.mse == min(ARMA.mse, na.rm=TRUE), 
                                      arr.ind=TRUE), collapse = ', ')
        final.seasonal[1,] <- paste(ARMA.arr[seasonal[1,1],seasonal[1,2]], collapse = ', ')
        final.df <- cbind(final.mse, final.ARMA, final.seasonal)
        colnames(final.df) <- c("Mean Squared Error", "ARMA Order (p-1, d, q-1)", 
                                "Seasonal Order (P-1, D, Q-1)")
      }
    }
  }
}

final.df

# Best model is SARIMA(0,1,1)(0,1,1)
mlmodel2 <- Arima(ts(data.all.trimmed$Investment, start=c(1961,1), frequency=4), 
                  order=c(0,1,1), seasonal=c(0,1,1))
forecastsmarima <- forecast(mlmodel2, h=20)

forecasts$SARIMA <- forecastsmarima$mean


# VAR MODEL:
data.var <- cbind(data.all.trimmed$Investment, data.all.trimmed$GDP)
data.var.trim <- data.var[2:156,]
colnames(data.var.trim) <- c("Investment", "GDP")
data.var.ts <- ts(data.var.trim, start=c(1961,2), frequency=4)

# Selecting optimal VAR
mse.vars <- matrix(NA, nrow=10,ncol=1)
forecastmatrix <- matrix(NA, nrow=20, ncol=10)

for(i in 1:10){
  var1 <- vars::VAR(data.var.ts, p=i, type='const')
  varforecast <- predict(var1, n.ahead=20)
  varforecastvalue <- varforecast$fcst[[1]]
  varforecastvalue <- varforecastvalue[,1]
  forecastmatrix[,i] <- varforecastvalue
  mse.vars[i,] <- mean((data.all$Investment[157:176]-forecastmatrix[,i])^2)
}
mse.vars
# Out of sample chose p=6

# makes png file, only used to make paper.
#png(file="C:/Users/Bret/OneDrive/U of A/Masters/ECON 509 - Time Series/Term Paper/
#LaTeX/Final Paper/Fig1.png")
plot(mse.vars, type='b', pch=19, xlab='VAR Lags', ylab='Mean Squared Error') +
  abline(v=which.min(mse.vars), col='red')
#dev.off()

var1 <- VAR(data.var.ts, p=6, type='const')
varforecast <- predict(var1, n.ahead=20)
varforecastvalue <- varforecast$fcst[[1]]
varforecastvalue <- varforecastvalue[,1]
forecasts$VAR <- varforecastvalue


# Plot:
ggplot(forecasts, aes(seq(as.Date("2000-01-01"), as.Date("2004-12-31"), "quarter"))) +
  geom_line(aes(y=TrueValues, colour="Actual Value"), size=1) +
  #geom_line(aes(y=MARIMA, color="MARIMA"), size=0.5) +
  geom_line(aes(y=SARIMA, color='SARIMA (MSE=908,002.7)'), size=1) +
  #geom_line(aes(y=SMARIMA, color="SMARIMA (MSE=2,306,940)"), size=1) +
  geom_line(aes(y=VAR, color='VAR (MSE=739,321.8)'), size=1) + 
  ylab("Investment (Millions)") + 
  xlab("Year") +
  theme(legend.title=element_blank()) +
  geom_hline(yintercept=0) + 
  geom_vline(xintercept=40000) + 
  theme(axis.line.x=element_line(color='black', size=0.75, linetype='solid'),
        axis.line.y=element_line(color='black', size=0.75, linetype='solid'),
        panel.background=element_rect('white', 'white', size=0.5),
        legend.key=element_rect(fill='transparent', color='transparent'))



errors <- data.frame(VAR=(forecasts$TrueValues-forecasts$VAR))
errors$SARIMA <- forecasts$TrueValues-forecasts$SARIMA


# Predictive power
dm.test(errors$VAR, errors$SARIMA, alternative='two.sided')
