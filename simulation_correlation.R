library(glmnet)
library(plyr)
library(reshape)
library(mvtnorm)
library(truncnorm)
library(rstiefel)
library(knockoff)
library(Matrix)

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

Trials_comparison <- function(nobs = 300, npar = 100, rho = 0.3, nsupp = 30, A = 3.5)
{
  ### Randomized X 
  covX <- abs(matrix(rep(1:npar,npar), ncol=npar) -  matrix(rep(1:npar,each = npar), ncol=npar))
  covX <- rho^covX
  X <- rmvnorm(nobs, mean = rep(0,npar), sigma = covX)
  
  X <- apply(X, 2, scale)/sqrt(nobs-1)
  
  ### Randomized Y
  eps <- rnorm(nobs, 0, 1)
  signs <- sample(c(-1,1), nsupp, replace=T)*A
  truesupp <- sample(npar, nuspp)
  Y <- X[,truesupp] %*% signs + eps
  
  ### Creating Knockoffs
  Sigma <- t(X) %*% X
  
  ## Compute s
  ## equi-correlated knockoffs
  lambda0 <- s <- min(eigen(Sigma)$values)
  s <- min(2*s-(1e-10),1-(1e-10))
  s_equiv <- rep(s, npar)
  invSigma <- chol2inv(chol(Sigma))
  #browser()
  C_equiv <- chol(2*diag(s_equiv) - diag(s_equiv) %*% invSigma %*% diag(s_equiv) )
  tildeU <- NullC(X)[,1:npar]
  tildeX_equiv <- X %*% (diag(npar) - invSigma %*% diag(s_equiv)) + tildeU %*% C_equiv
  
  #### Knockoff
  model <- glmnet(cbind(X, tildeX_equiv), Y)
  lambdas <- max(model$lambda)
  lambdas <- seq(0.001*lambdas, lambdas, 0.001*lambdas)
  model_equiv <- glmnet(cbind(X, tildeX_equiv), Y, lambda = lambdas)
  
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
  Thres <- thres_compute(W_equiv, type="knockoff")
  index_equiv <- which(W_equiv>= Thres)
  Thres <- thres_compute(W_equiv, type="knockoff+")
  index_equiv_plus <- which(W_equiv>= Thres)
  fdr_trial <- sum(!(index_equiv %in% truesupp))/max(length(index_equiv),1)
  power_knockoff <- sum(index_equiv %in% truesupp)/nsupp
  fdr_trial_plus<- sum(!(index_equiv_plus %in% truesupp))/max(length(index_equiv_plus),1)
  power_knockoff_plus <- sum(index_equiv_plus %in% truesupp)/nsupp
  
  #browser()
  ############### BHq
  
  model_ls <- lm(Y~X)
  Z_ls <- summary(model_ls)$coef[-1, "t value"]
  Thres <- thres_compute(Z_ls, type="BHq")
  index.BHq <- which(abs(Z_ls) >= Thres)
  fdr_trial_BHq <- sum(!(index.BHq %in% truesupp))/max(length(index.BHq),1)
  power_BHq <- sum(index.BHq %in% truesupp)/nsupp
  return(c(fdr_trial, fdr_trial_plus, fdr_trial_BHq, 
           power_knockoff, power_knockoff_plus, power_BHq))
}

result_correlation <- matrix(nrow =2000, ncol=7)
result_correlation <- as.data.frame(result_correlation)
names(result_correlation) <- c("correlation","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                               "power_knockoff","power_knockoff_plus", "power_BHq")
for(i in 1:10)
{
  result_correlation[(i*200-199):(i*200),1] <- i/10-0.1
  for(j in 1:200)
  {
    result<- Trials_comparison(nobs = 3000, npar = 1000, rho = i/10-0.1, nsupp = 10*i)
    result_correlation[i*200-200+j,2:7] <- result
    gc()
    if(j%%10==0){cat("Trial ",j," finished. ",as.character(Sys.time()), "\n", sep="")}
  }
}

write.csv(result_correlation, "result_correlation.csv")



