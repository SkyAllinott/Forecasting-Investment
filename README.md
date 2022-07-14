# Forecasting Investment
## Overview
I construct the aggregate investment from aggregate capital (K) according to the formula:

I = K-L(K)

Where L(*) is the lag operator. I then construct two forecasts
1. Univariate ARIMA model
2. Bivariate VAR model (second variable is GDP)

## Results
![Figure](https://user-images.githubusercontent.com/52394699/179077453-5b3d3d9f-f0fb-4c43-9981-7d5aafbe3489.png)

The bivariate VAR outperforms the univariate ARIMA model, and the Diebold-Mariano test for forecast accuracy states the VAR model is statistically significantly better.
