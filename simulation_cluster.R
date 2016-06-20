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
library(knockoff)
library(CVXfromR)
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
Trials <- function(nobs = 300, npar = 100, rho = 0.3, nsupp = 30, A = 3.5)
{
  ### Randomized X 
  v <- rnorm(nobs, 0, 1)
  V <- matrix(rep(v, npar), ncol=npar)
  x <- rnorm(nobs * npar, 0, 1)
  X <- matrix(x, ncol = npar)
  
  X <- rho*V + sqrt(1-rho^2) * X
  X <- apply(X, 2, scale)/sqrt(nobs-1)
  
  ### Randomized Y
  eps <- rnorm(nobs, 0, 1)
  signs <- sample(c(-1,1), nsupp, replace=T)*A
  Y <- X[,1:nsupp] %*% signs + eps
  
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
  
  cvxcode <- paste("variable s(npar)",
                   "minimize sum(-s)",
                   "subject to",
                   "0 <= s <= 1",
                   "2*Sigma - diag(s)  == semidefinite(npar)", sep=";")
  s_sdp <- CallCVX(cvxcode, const.vars = list(npar = npar, Sigma=Sigma), opt.var.names = "s")
  s_sdp <- s_sdp$s
  #browser()
  s_sdp <- s_sdp
  #browser()
  C_sdp <- chol(2*diag(s_sdp) - diag(s_sdp) %*% invSigma %*% diag(s_sdp) + 1e-5* diag(npar))
  tildeX_sdp <- X %*% (diag(npar) - invSigma %*% diag(s_sdp)) + tildeU %*% C_sdp
  
  
  ### Create Permutations
  index <- sample(nobs)
  permuteX <- X[index, ]
  
  
  #### Knockoff
  model <- glmnet(cbind(X, tildeX_equiv), Y)
  lambdas <- max(model$lambda)
  lambdas <- seq(0.001*lambdas, lambdas, 0.001*lambdas)
  model_equiv <- glmnet(cbind(X, tildeX_equiv), Y, lambda = lambdas)
  model <- glmnet(cbind(X, tildeX_sdp), Y)
  lambdas <- max(model$lambda)
  lambdas <- seq(0.001*lambdas, lambdas, 0.001*lambdas)
  model_sdp <- glmnet(cbind(X, tildeX_sdp), Y, lambda = lambdas)
  
  model.permute <- glmnet(cbind(X, permuteX), Y, lambda = lambdas)
  Z_equiv <- numeric(2*npar)
  Z_sdp <- numeric(2*npar)
  Z.permute <- numeric(2*npar)
  for(i in 1:(2*npar))
  {
    temp <- which(abs(model_equiv$beta[i,])>=1e-7)
    if(length(temp) == 0) {Z[i] <- 0}
    else
    {
      temp <- min(temp)
      Z_equiv[i] <- model_equiv$lambda[temp]
    }
    temp <- which(abs(model_sdp$beta[i,])>=1e-7)
    if(length(temp) == 0) {Z[i] <- 0}
    else
    {
      temp <- min(temp)
      Z_sdp[i] <- model_sdp$lambda[temp]
    }
    temp <- which(abs(model.permute$beta[i,])>=1e-7)
    if(length(temp) == 0) {Z.permute[i] <- 0}
    else
    {
      temp <- min(temp)
      Z.permute[i] <- model.permute$lambda[temp]
    }
  }
  Z_equiv <- matrix(Z_equiv, ncol=2)
  Z_sdp <- matrix(Z_sdp, ncol=2)
  Z.permute <- matrix(Z.permute, ncol=2)
  W_equiv <- apply(Z_equiv, 1, max) * sign(Z_equiv[,1]- Z_equiv[,2])
  W_sdp <- apply(Z_sdp, 1, max) * sign(Z_sdp[,1]- Z_sdp[,2])
  W.permute <- apply(Z.permute, 1, max) * sign(Z.permute[,1]- Z.permute[,2])
  Thres <- thres_compute(W_equiv, type="knockoff")
  index_equiv <- which(W_equiv>= Thres)
  Thres <- thres_compute(W_equiv, type="knockoff+")
  index_equiv_plus <- which(W_equiv>= Thres)
  Thres <- thres_compute(W_sdp, type="knockoff")
  index_sdp <- which(W_sdp>= Thres)
  Thres <- thres_compute(W_sdp, type="knockoff+")
  index_sdp_plus <- which(W_sdp>= Thres)
  Thres <- thres_compute(W.permute, type="knockoff+")
  index.permute <- which(W.permute>= Thres)
  fdr_trial <- sum(!(index_equiv %in% (1:nsupp)))/max(length(index_equiv),1)
  power_knockoff <- sum(index_equiv %in% (1:nsupp))/nsupp
  fdr_trial_plus<- sum(!(index_equiv_plus %in% (1:nsupp)))/max(length(index_equiv_plus),1)
  power_knockoff_plus <- sum(index_equiv_plus %in% (1:nsupp))/nsupp
  fdr_trial_sdp <- sum(!(index_sdp %in% (1:nsupp)))/max(length(index_sdp),1)
  power_knockoff_sdp <- sum(index_sdp %in% (1:nsupp))/nsupp
  fdr_trial_sdp_plus <- sum(!(index_sdp_plus %in% (1:nsupp)))/max(length(index_sdp_plus),1)
  power_knockoff_sdp_plus <- sum(index_sdp_plus %in% (1:nsupp))/nsupp
  fdr_trial_permute <- sum(!(index.permute %in% (1:nsupp)))/max(length(index.permute),1)
  power_permute <- sum(index.permute %in% (1:nsupp))/nsupp
  
  #browser()
  ############### BHq
  
  model_ls <- lm(Y~X)
  Z_ls <- summary(model_ls)$coef[-1, "t value"]
  Thres <- thres_compute(Z_ls, type="BHq")
  index.BHq <- which(abs(Z_ls) >= Thres)
  Thres <- thres_compute(Z_ls, type="BHq_corrected")
  index.BHq.corrected <- which(abs(Z_ls) >= Thres)
  invSigma <- chol2inv(chol(Sigma))
  Z_ls_whitened <- rmvnorm(1,sigma = 1.01/(lambda0)*diag(npar)- invSigma )
  Z_whitened <- (model_ls$coef[-1] + Z_ls_whitened)*sqrt(lambda0)/sqrt(1.01)
  #browser()
  Thres <- thres_compute(Z_whitened, type="BHq")
  index.BHq.whitened <- which(abs(Z_whitened) >= Thres)
  fdr_trial_BHq <- sum(!(index.BHq %in% (1:nsupp)))/max(length(index.BHq),1)
  power_BHq <- sum(index.BHq %in% (1:nsupp))/nsupp
  fdr_trial_BHq_Corrected <- sum(!(index.BHq.corrected %in% (1:nsupp)))/max(length(index.BHq.corrected),1)
  power_BHq_Corrected <- sum(index.BHq.corrected %in% (1:nsupp))/nsupp
  fdr_trial_BHq_whitened <- sum(!(index.BHq.whitened %in% (1:nsupp)))/max(length(index.BHq.whitened),1)
  power_BHq_Whitened <- sum(index.BHq.whitened %in% (1:nsupp))/nsupp
  return(c(fdr_trial, fdr_trial_plus, fdr_trial_permute, fdr_trial_sdp,fdr_trial_sdp_plus, fdr_trial_BHq, fdr_trial_BHq_Corrected, fdr_trial_BHq_whitened,
           power_knockoff, power_knockoff_plus, power_permute, power_knockoff_sdp, power_knockoff_sdp_plus,  power_BHq, power_BHq_Corrected, power_BHq_Whitened))
}

result.data <- matrix(nrow =600, ncol=16)
result.data <- as.data.frame(result.data)
names(result.data) <- c("fdr_trial", "fdr_trial_plus", "fdr_trial_permute","fdr_trial_sdp", "fdr_trial_sdp_plus", 
                        "fdr_trial_BHq", "fdr_trial_BHq_Corrected", "fdr_trial_BHq_whitened",
                        "power_knockoff","power_knockoff_plus", "power_permute", "power_knockoff_sdp", "power_knockoff_sdp_plus",
                        "power_BHq", "power_BHq_Corrected", "power_BHq_Whitened")
for(i in 1:140)
{
  result<- Trials(nobs = 3000, npar = 1000, rho = 0, nsupp = 30)
  result.data[i,] <- result
  if(i%%10==0){cat("Trial ",i," finished. ",as.character(Sys.time()), "\n", sep="")}
}
save.image("~/result.RData")

