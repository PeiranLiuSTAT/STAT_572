library(glmnet)
library(ggplot2)
library(xtable)
library(plyr)
library(reshape)
library(mvtnorm)
library(gridExtra)
library(truncnorm)
library(grid)
library(rstiefel)
library(Matrix)
setwd("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Codes")


canonical.svd <- function(X) 
{
  X.svd <- svd(X)
  for (j in 1:min(dim(X))) {
    i = which.max(abs(X.svd$u[,j]))
    if (X.svd$u[i,j] < 0) {
      X.svd$u[,j] <- -X.svd$u[,j]
      X.svd$v[,j] <- -X.svd$v[,j]
    }
  }
  return(X.svd)
}

normc <- function(X) 
{
  X.scaled = scale(X, center=FALSE, scale=sqrt(colSums(X^2)))
  X.scaled[,] # No attributes
}
decompose.X <- function(X, randomize) 
{
  n = nrow(X); p = ncol(X)
  stopifnot(n >= 2*p)
  
  result = canonical.svd(X)
  Q = qr.Q(qr(cbind(result$u, matrix(0,n,p))))
  u_perp = Q[,(p+1):(2*p)]
  if (randomize) {
    Q = qr.Q(qr(rnorm_matrix(p,p)))
    u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}

thres_compute <- function(W, q=0.2, type="knockoff")
{
  temp <- sort(unique(abs(W)))[-1]
  k <- length(temp)
  thres <- Inf
  for(i in 1:k)
  {
    numerator <- sum(W <= -temp[i])
    denumerator <- max(sum(W >= temp[i]),1)
    #
    if(type == "knockoff")
    {
      if((numerator)/denumerator <= q)
      {
        thres <- temp[i]
        break
      }
    }
    else if(type == "knockoff+")
    {
      if((1+numerator)/denumerator <= q)
      {
        thres <- temp[i]
        break
      }
    }
    else if(type == "BHq")
    {
      numerator <- length(W) * 2 * pnorm(-temp[i])
      denumerator <- sum(abs(W) >= temp[i])
      if(numerator / denumerator <= q)
      {
        thres <- temp[i]
        break
      }
    }
    else if(type == "BHq_corrected")
    {
      numerator <- length(W) * 2 * pnorm(-temp[i])
      denumerator <- sum(abs(W) >= temp[i])
      if(numerator / denumerator <= q/(log(length(W))-digamma(1)))
      {
        thres <- temp[i]
        break
      }
    }
  }
  return(thres)
}

Trials_figure6 <- function(A = 3.5, nobs = 2000, npar = 1000, nsupp = 200)
{
  ### Randomized X 
  x <- rnorm(nobs * npar, 0, 1)
  X <- matrix(x, ncol = npar)
  X.svd <- svd(X)
  X <- X.svd$u
  
  ### Randomized Y
  eps <- rnorm(nobs, 0, 1)
  signs <- sample(c(-1,1), nsupp, replace=T)*A
  Y <- X[,1:nsupp] %*% signs + eps
  
  X <- normc(X)
  npar <- ncol(X)
  Sigma <- t(as.matrix(X)) %*% as.matrix(X)
  
  X.svd <- decompose.X(X,FALSE)
  ## Compute s
  ## equi-correlated knockoffs
  lambda0 <- lambda_min <- min(X.svd$d)^2
  s <- min(2*lambda_min, 1)
  s_equiv <- rep(s, npar)
  s_diff = pmax(0, 2*s - (s/X.svd$d)^2) # can be negative due to numerical error
  invSigma <- (X.svd$v %*% diag(1/X.svd$d)^2) %*% t(X.svd$v)
  C_equiv <- diag(sqrt(s_diff)) %*% t(X.svd$v)
  #tildeU <- NullC(as.matrix(X))[,1:npar]
  
  tildeU <- X.svd$u_perp
  tildeX_equiv <- as.matrix(X) %*% (diag(npar) - invSigma %*% diag(s_equiv)) + tildeU %*% C_equiv
  #tildeX_equiv_package <- knockoff.create(X, method = "equicorrelated")
  ####### Fit Augmented Lasso
  model <- glmnet(cbind(as.matrix(X), tildeX_equiv), Y, intercept=F, standardize=F,
                  standardize.response=F)
  lambda_max <- max(model$lambda)
  lambda_min <- lambda_max/2e3
  ks <- (0:(5*npar-1))/(5*npar)
  lambdas <- lambda_max * (lambda_min/lambda_max)^ks
  model_equiv <- glmnet(cbind(as.matrix(X), tildeX_equiv), Y, lambda = lambdas, intercept=F, standardize=F,
                        standardize.response=F)
  Z_equiv <- numeric(2*npar)
  for(i in 1:(2*npar))
  {
    temp <- which(abs(model_equiv$beta[i,])>=1e-7)
    if(length(temp) == 0) {Z_equiv[i] <- 0}
    else
    {
      temp <- min(temp)
      Z_equiv[i] <- model_equiv$lambda[temp]
    }
  }
  Z_equiv <- matrix(Z_equiv, ncol=2)
  W_equiv <- apply(Z_equiv, 1, max) * sign(Z_equiv[,1]- Z_equiv[,2])

  Thres <- thres_compute(W_equiv, type="knockoff+", q=q)
  index_equiv_plus <- which(W_equiv>= Thres)
  
  ## Just use package
  fdr_trial_plus<- sum(!(index_equiv_plus %in% (1:nsupp)))/max(length(index_equiv_plus),1)
  power_knockoff_plus <- sum(index_equiv_plus %in% (1:nsupp))/nsupp
  #browser()
  ############### BHq
  
  model_ls <- lm(Y~X)
  Z_ls <- summary(model_ls)$coef[-1, "t value"]
  Thres <- thres_compute(Z_ls, type="BHq")
  index.BHq <- which(abs(Z_ls) >= Thres)
  fdr_trial_BHq <- sum(!(index.BHq %in% (1:nsupp)))/max(length(index.BHq),1)
  power_BHq <- sum(index.BHq %in% (1:nsupp))/nsupp
  
  return(c(fdr_trial_plus, fdr_trial_BHq, power_knockoff_plus, power_BHq))
}

  result.data <- matrix(nrow =230, ncol=5)
  result.data <- as.data.frame(result.data)
  names(result.data) <- c("Magnitude","fdr_trial_plus", "fdr_trial_BHq", "power_knockoff_plus", 
                          "power_BHq")
  
  As <- rep(seq(1,6.5,0.25), each=100)
  
  
  for(i in 1:230)
  {
    result.data[i,1] <- As[i]
    result.data[i,2:5] <- Trials_figure6(A=As[i])
    gc()
    if(j%%10==0){cat("Trial ",j," finished. ",as.character(Sys.time()), "\n", sep="")}
  }
  
  write.csv(result.data, "figure6.csv")

result_fig6 <- read.csv("figure6_1.csv", header=T)[,-1]
for(i in 2:10)
{
  filename <- paste("figure6_",i,".csv",sep="")
  result_fig6 <- rbind(result_fig6, read.csv(filename, header=T)[,-1])
}
result_fig6[,1] <- As
result_mean_fig6 <- matrix(nrow =23, ncol=5)
result_mean_fig6 <- as.data.frame(result_mean_fig6)
names(result_mean_fig6) <- c("Magnitude","fdr_trial_plus", "fdr_trial_BHq", "power_knockoff_plus", 
                        "power_BHq")
for(i in 1:23)
{
  result_mean_fig6[i,] <- apply(result_fig6[(i*250-249):(i*250),], 2, mean)
}

pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "compare_BHq", ".pdf", sep=""), width=10,height=5)
par(mfrow= c(1,2))
plot(result_mean_fig6$Magnitude, result_mean_fig6$fdr_trial_plus*100, type="l", col="red",xlab="Signal Magnitude", ylab="FDR(%)", ylim=c(0,25))
lines(result_mean_fig6$Magnitude, result_mean_fig6$fdr_trial_BHq*100, col="blue")
abline(h=20, col="black", lty=2)
points(result_mean_fig6$Magnitude, result_mean_fig6$fdr_trial_plus*100, col="red", pch=17)
points(result_mean_fig6$Magnitude, result_mean_fig6$fdr_trial_BHq*100, col="blue", pch=17)
legend("bottomright", c("Nominal", "Knockoff+", "BHq"), col=c("black","red", "blue"), lwd = c(2.5, 2.5,2.5),lty=c(2,1,1))

plot(result_mean_fig6$Magnitude, result_mean_fig6$power_knockoff_plus*100, type="l", col="red",xlab="Signal Magnitude", ylab="Power(%)", ylim=c(0,100))
lines(result_mean_fig6$Magnitude, result_mean_fig6$power_BHq*100, col="blue")
points(result_mean_fig6$Magnitude, result_mean_fig6$power_knockoff_plus*100, col="red", pch=17)
points(result_mean_fig6$Magnitude, result_mean_fig6$power_BHq*100, col="blue", pch=17)
legend("bottomright", c("Nominal", "Knockoff+", "BHq"), col=c("black","red", "blue"), lwd = c(2.5, 2.5,2.5),lty=c(2,1,1))

dev.off()