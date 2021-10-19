#' ---
#' title: "Assignment1"
#' author: "Devosmita Chatterjee"
#' date: "2 April 2021"
#' output:
#'   html_document: default
#' ---
#' 
## ----setup, include<-FALSE-----------------------------------------------
knitr::opts_chunk$set(echo <- TRUE)

#' 
#' 
## ------------------------------------------------------------------------
rm(list = ls())# Clear workspace environment
cat("\014")# Clear console

#' 
#' ##### Practical part 2. #####
## ------------------------------------------------------------------------
# Read the given data intel.csv.
data <- read.csv("C:\\Users\\devcha\\Desktop\\TMS088\\intel.csv",header<-TRUE)
data

#' 
## ------------------------------------------------------------------------
# Number of datapoints
length(data[,1])
# Number of missing values in VolumeMissing
sum(is.na(data[,2]))

#' 
## ------------------------------------------------------------------------
# Compute and plot the time series X where X[t] are the data points in VolumeMissing.
#png("TimeSeries.png")
X <- data[,2]# VolumeMissing
plot(X, type = "l", main = "Time Series Plot", sub =  expression("Daily trading volume data for Intel Corporation stock at Nasdaq."), xlab = "day", ylab = "daily trading volume data X[t]",
col.main = "red", col.sub = "red", col.lab = "blue")
#dev.off()

#' 
## ------------------------------------------------------------------------
# Compute and plot the differenced log time series Y[t] = X[t+1] - X[t] where X[t] are the data points in VolumeMissing.
#png("logTimeSeries.png")
X <- data[,2]# VolumeMissing
Y <- matrix(, nrow <- length(X)-1, ncol <- 1)
for (i in 1:length(X)-1) 
{
  Y[i] <- log(X[i+1]) - log(X[i])
}
plot(Y, type = "l", main = "Time Series Plot", sub = expression("Log-returns of daily trading volume data for Intel Corporation stock at Nasdaq."), xlab = "day", ylab = "log-returns Y[t] <- log(X[t+1]) - log(X[t])", col.main = "red", col.sub = "red", col.lab = "blue")
#dev.off()

#' 
## ------------------------------------------------------------------------
# Compute the sample autocorrelation function (ACF) of the log-returns.
L <- acf(Y, lag = 20, type = "correlation", na.action = na.pass, plot = FALSE)
ACF <- L$acf

# Compute whether autocorrelations are statistically significant.
for (i in 1:21)
{
if (ACF[i] > 1.96/sqrt(4958))
{
print(i)
}
else if (ACF[i] < -1.96/sqrt(4958))
{
print(i)
}
end
}
end

#' 
#' 
## ------------------------------------------------------------------------
# Plot the sample autocorrelation function (ACF) of the log-returns.
#png("ACF.png")
plot(L, xaxt = 'no', main = "", sub  = expression("Sample ACF of the log-returns h = 0,...,20."), ylab = "ACF of log-returns", xlab = "Lag", col.sub = "red", col.lab = "blue", axes = FALSE, xaxt = 'n', yaxt = 'n')
par(new = TRUE)
plot(L$acf, main = "ACF Plot", xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', col.main = "red", pch = 20)
axis(1, at = seq(1,21,1), labels = seq(0,20,1), cex.axis = 0.6)
axis(2, at = seq(-1,1,0.1), labels = seq(-1,1,0.1), cex.axis = 0.6)
par(las = 2)
axis(2, at = 1.96/sqrt(4958), labels = expression(+1.96/sqrt(N)), cex.axis = 0.5, col.axis = "green")
axis(2, at = -1.96/sqrt(4958), labels = expression(-1.96/sqrt(N)), cex.axis = 0.5, col.axis = "green")
#dev.off()

#' 
#' ##### Practical part 3. #####
## ------------------------------------------------------------------------
# Set of indices of the missing values of Y.
M <- which(is.na(Y))

# Initialize best linear predictors for all missing values.
bestlinearpredictor <- matrix(, length(M), 1)

# Loop starts
for (j in 1:length(M))
{
  # Index of missing value
  t <- M[j]
  
  # Lag
  q <- 4
  
  # Total number of data points
  N <- 4958-1
  
  # Maximum
  max <- max(1,t-q)
  
  # Minimum
  min <- min(N,t+q)
  
  # Sequence of integers
  int <- seq.int(max,min)
  
  # Sequence of integers which are missing and can not be used for prediction
  drop <- intersect(int,M)

  # Set of indices without missing values
  s <- int[!int %in% drop]

  # Length of set of indices without missing values
  n <- length(Y[s])
  
  # Compute autocorrelation function for Y[s].
  L <- acf(Y[s], lag = n-1, type = "correlation", plot = FALSE)

  # Construct left hand side of linear system of equations from corollary 2.4.6.
  # Construct covariance matrix
  gamma <- data.matrix(L$acf)
  s1 <- dim(gamma)
  Gamma <- matrix(1, s1, s1)
  d <- row(Gamma) - col(Gamma)
  for (i in 1:(s1 - 1))
  {
    Gamma[d == i | d == (-i)] <- gamma[i + 1]
  }
  end
  
  # Construct right hand side of linear system of equations from corollary 2.4.6.
  ACF <- data.matrix(L$acf)
  b <- matrix(, length(s), 1)
  for (k in length(s):-1:1)
  {
    c <- abs(t-s[k])
    b[k] <- ACF[c]
  }
  
  # Solve the linear system of equations from corollary 2.4.6 to find the coefficients.
  A <- Gamma
  a1 <- solve(A, b)

  # Find mean
  result.mean <-  mean(Y[-M])

  # Compute first coefficient
  a0 <- result.mean*(1-colSums(a1))
  
  # Compute best linear predictors for each missing value
  bestlinearpredictor[j] <- a0+t(Y[s])%*%a1 

}
bestlinearpredictor

#' 
#' ##### Practical part 4. #####
## ------------------------------------------------------------------------
# Compute the differenced log time series Y2[t] = X2[t+1] - X2[t] where X2[t] are the data points in Volume.
X2 <- data[,1]# Volume
Y2 <- matrix(, nrow = length(X2)-1, ncol = 1)
for (i in 1:length(X2)-1) 
{
  Y2[i] <- log(X2[i+1]) - log(X2[i])
}

# Compute modified Y where each missing value is replaced by best linear predictor.
Ycap <- Y
Ycap[M] <- bestlinearpredictor
RMSE <- sqrt((1/length(M))*colSums(data.matrix((Y2[M]-Ycap[M])^2)))

# Compute modified Y where each missing value is replaced by a value using simple linear interpolation.
#install.packages("baytrends")
#install.packages("gapfill")
Ybar <- fillMissing(Y)
RMSE2 <- sqrt((1/length(M))*colSums(data.matrix((Y2[M]-Ybar[M])^2)))

