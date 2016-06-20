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

install_url("http://www.bscb.cornell.edu/~bien/CVXfromR/CVXfromR_1.6.tar.gz")

setwd("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Codes")

{v <- rnorm(3000, 0, 1)
V <- matrix(rep(v, 1000), ncol=1000)
x <- rnorm(3000000, 0, 1)
X <- matrix(x, ncol = 1000)

X <- 0*V + sqrt(1) * X
X <- apply(X, 2, scale)/sqrt(2999)
eps <- rnorm(3000, 0, 1)
Y <- X[,1:30] %*% rep(3.5, 30) + eps


### Creating Knockoffs
Sigma <- t(X) %*% X
 
## Compute s
## equi-correlated knockoffs
s <- eigen(Sigma)$values[1000]
s <- min(2*s-(1e-10),1-(1e-10))
s <- rep(s, 1000)
C <- chol(2*diag(s) - diag(s) %*% solve(Sigma) %*% diag(s))
tildeU <- NullC(X)[,1:1000]
tildeX <- X %*% (diag(1000) - solve(Sigma, diag(s))) + tildeU %*% C}
{
model <- glmnet(cbind(X, tildeX), Y)
lambdas <- max(model$lambda)
lambdas <- seq(0.001*lambdas, lambdas, 0.001*lambdas)
model <- glmnet(cbind(X, tildeX), Y, lambda = lambdas)
Z <- numeric(2000)
setup.dir = "C:/Users/Peiran Liu/Softwares/cvx"
cvxcode <- paste("variable s(npar)",
                 "minimize sum(-s)",
                 "subject to",
                 "0 <= s <= 1",
                 "2*Sigma - diag(s)  == semidefinite(npar)", sep=";")
s_sdp <- CallCVX(cvxcode, const.vars = list(npar = npar, Sigma=Sigma), opt.var.names = "s")

for(i in 1:2000)
{
  temp <- which(abs(model$beta[i,])>=1e-7)
  if(length(temp) == 0) {Z[i] <- 0}
  else
  {
    temp <- min(temp)
    Z[i] <- model$lambda[temp]
  }
}
Z <- matrix(Z, ncol=2)
W <- apply(Z, 1, max) * sign(Z[,1]- Z[,2])}
Thres <- thres_compute(W, type="knockoff+")
index <- which(W>=Thres)

m1 <- knockoff.filter(X, Y, normalize = F, threshold = "knockoff+")
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
  X <- apply(X, 2, scale, center=F)/sqrt(nobs-1)
  
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
    if(length(temp) == 0) {Z_equiv[i] <- 0}
    else
    {
      temp <- min(temp)
      Z_equiv[i] <- model_equiv$lambda[temp]
    }
    temp <- which(abs(model_sdp$beta[i,])>=1e-7)
    if(length(temp) == 0) {Z_sdp[i] <- 0}
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

Trials.permute <- function(nobs = 300, npar = 100, rho = 0.3, nsupp = 30, A = 3.5)
{
  ### Randomized X 
  CovX <- rho * matrix(1, nrow=npar, ncol=npar) + (1-rho)*diag(npar)
  X <- rmvnorm(nobs, mean = rep(0, npar), sigma = CovX)
  X <- apply(X, 2, scale)/sqrt(nobs-1)
  
  ### Randomized Y
  eps <- rnorm(nobs, 0, 1)
  Y <- X[,1:nsupp] %*% rep(A, nsupp) + eps
  
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
  ### Create Permutations
  index <- sample(nobs)
  permuteX <- X[index, ]
  
  
  #### Knockoff
  model <- glmnet(cbind(X, tildeX_equiv), Y)
  lambdas <- max(model$lambda)
  lambdas <- seq(0.001*lambdas, lambdas, 0.001*lambdas)
  model_equiv <- glmnet(cbind(X, tildeX_equiv), Y, lambda = lambdas)
  
  model.permute <- glmnet(cbind(X, permuteX), Y, lambda = lambdas)
  Z_equiv <- numeric(2*npar)
  Z.permute <- numeric(2*npar)
  for(i in 1:(2*npar))
  {
    temp <- which(abs(model_equiv$beta[i,])>=1e-7)
    if(length(temp) == 0) {Z_equiv[i] <- 0}
    else
    {
      temp <- min(temp)
      Z_equiv[i] <- model_equiv$lambda[temp]
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
  Z.permute <- matrix(Z.permute, ncol=2)
  W_equiv <- apply(Z_equiv, 1, max) * sign(Z_equiv[,1]- Z_equiv[,2])
  W.permute <- apply(Z.permute, 1, max) * sign(Z.permute[,1]- Z.permute[,2])
  Thres <- thres_compute(W_equiv, type="knockoff")
  index_equiv <- which(W_equiv>= Thres)
  Thres <- thres_compute(W.permute, type="knockoff")
  index.permute <- which(W.permute>= Thres)
  fdr_trial <- sum(!(index_equiv %in% (1:nsupp)))/max(length(index_equiv),1)
  fdr_trial_permute <- sum(!(index.permute %in% (1:nsupp)))/max(length(index.permute),1)
  
  #browser()
  ############### BHq
  return(c(fdr_trial, fdr_trial_permute))
}
fdr_trial <- numeric(600)
power_knockoff_plus <- numeric(600)
fdr_trial_permute <- numeric(600)
power_permute <- numeric(600)
fdr_trial_BHq <- numeric(600)
power_BHq <- numeric(600)
fdr_trial_BHq_Corrected <- numeric(600)
power_BHq_Corrected
fdr_trial_BHq_whitened <- numeric(600)
result.data <- matrix(nrow =600, ncol=16)
result.data <- as.data.frame(result.data)
names(result.data) <- c("fdr_trial", "fdr_trial_plus", "fdr_trial_permute","fdr_trial_sdp", "fdr_trial_sdp_plus", 
                        "fdr_trial_BHq", "fdr_trial_BHq_Corrected", "fdr_trial_BHq_whitened",
                        "power_knockoff","power_knockoff_plus", "power_permute", "power_knockoff_sdp", "power_knockoff_sdp_plus",
                        "power_BHq", "power_BHq_Corrected", "power_BHq_Whitened")
for(i in 541:600)
{
  result<- Trials(nobs = 3000, npar = 1000, rho = 0, nsupp = 30)
  result.data[i,] <- result
  if(i%%10==0){cat("Trial ",i," finished. ",as.character(Sys.time()), "\n", sep="")}
}

result.permute <- matrix(nrow =1000, ncol=2)
result.permute <- as.data.frame(result.permute)
names(result.permute) <- c("fdr_trial", "fdr_trial_permute")
for(i in 1:1000)
{
  result<- Trials.permute(nobs = 300, npar = 100, rho = 0.3, nsupp = 30)
  result.permute[i,] <- result
  gc()
  if(i%%10==0){cat("Trial ",i," finished. ",as.character(Sys.time()), "\n", sep="")}
}

############33 Comparison
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
  index <- sample()
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
  fdr_trial <- sum(!(index_equiv %in% (1:nsupp)))/max(length(index_equiv),1)
  power_knockoff <- sum(index_equiv %in% (1:nsupp))/nsupp
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
  return(c(fdr_trial, fdr_trial_plus, fdr_trial_BHq, 
           power_knockoff, power_knockoff_plus, power_BHq))
}

result_sparsity <- matrix(nrow =4000, ncol=7)
result_sparsity <- as.data.frame(result_sparsity)
names(result_sparsity) <- c("sparsity","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                        "power_knockoff","power_knockoff_plus", "power_BHq")
for(i in 1:20)
{
  result_sparsity[(i*200-199):(i*200),1] <- i*10
  for(j in 1:200)
  {
    result<- Trials_comparison(nobs = 3000, npar = 1000, rho = 0, nsupp = 10*i)
    result_sparsity[i*200-200+j,2:7] <- result
    if(j%%10==0){cat("Trial ",j," finished. ",as.character(Sys.time()), "\n", sep="")}
  }
}

result_sparsity <- read.csv("result_sparsity4.csv")[,-1]
result_sparsity <- rbind(result_sparsity, read.csv("result_sparsity5.csv")[,-1])


result_amplitude <- matrix(nrow =3000, ncol=7)
result_amplitude <- as.data.frame(result_amplitude)
names(result_amplitude) <- c("amplitude","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                             "power_knockoff","power_knockoff_plus", "power_BHq")
for(i in 1:15)
{
  result_amplitude[(i*200-199):(i*200),1] <- 2.7+0.1*i
  for(j in 1:10)
  {
    result<- Trials_comparison(nobs = 3000, npar = 1000, rho = 0, nsupp = 30, A = 2.7+0.1*i)
    result_amplitude[i*200-200+j,2:7] <- result
    if(j%%10==0){cat("Trial ",j," finished. ",as.character(Sys.time()), "\n", sep="")}
  }
}
result_amplitude <- read.csv("result_amplitude.csv", header=T)
result_amplitude <- result_amplitude[,-1]


result_correlation <- matrix(nrow =2000, ncol=7)
result_correlation <- as.data.frame(result_correlation)
names(result_correlation) <- c("correlation","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                            "power_knockoff","power_knockoff_plus", "power_BHq")
for(i in 1:10)
{
  result_correlation[(i*200-199):(i*200),1] <- i/10-0.1
  for(j in 1:10)
  {
    result<- Trials_comparison(nobs = 3000, npar = 1000, rho = i/10-0.1, nsupp = 10*i)
    result_correlation[i*200-200+j,2:7] <- result
    if(j%%10==0){cat("Trial ",j," finished. ",as.character(Sys.time()), "\n", sep="")}
  }
}

result_correlation <- read.csv("result_correlation.csv", header=T)
result_correlation <- result_correlation[,-1]

result_mean_sparsity <- matrix(nrow =20, ncol=7)
result_mean_sparsity <- as.data.frame(result_mean_sparsity)
names(result_mean_sparsity) <- c("sparsity","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                                 "power_knockoff","power_knockoff_plus", "power_BHq")

for(i in 1:20)
{
  result_mean_sparsity[i,1] <- i*10
  result_mean_sparsity[i,2:7] <- apply(result_sparsity[(i*200-199):(i*200),2:7],2,mean)
}

result_mean_correlation <- matrix(nrow =10, ncol=7)
result_mean_correlation <- as.data.frame(result_mean_correlation)
names(result_mean_correlation) <- c("correlation","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                                 "power_knockoff","power_knockoff_plus", "power_BHq")

for(i in 1:10)
{
  result_mean_correlation[i,1] <- 0.1*i-0.1
  result_mean_correlation[i,2:7] <- apply(result_correlation[(i*200-199):(i*200),2:7],2,mean)
}



result_mean_amplitude <- matrix(nrow =10, ncol=7)
result_mean_amplitude <- as.data.frame(result_mean_amplitude)
names(result_mean_amplitude) <- c("amplitude","fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                                    "power_knockoff","power_knockoff_plus", "power_BHq")

for(i in 1:15)
{
  result_mean_amplitude[i,1] <- 2.7+0.1*i
  result_mean_amplitude[i,2:7] <- apply(result_amplitude[(i*200-199):(i*200),2:7],2,mean)
}


pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Slides Tex/Img/", 
                 "amplitude", ".pdf", sep=""), width=10,height=5)
par(mfrow= c(1,2))
plot(result_mean_amplitude$amplitude, result_mean_amplitude$fdr_trial*100, type="l", col="blue",xlab="amplitude level", ylab="FDR(%)", ylim=c(0,30))
lines(result_mean_amplitude$amplitude, result_mean_amplitude$fdr_trial_plus*100, col="red")
lines(result_mean_amplitude$amplitude, result_mean_amplitude$fdr_trial_BHq*100, col="black")
abline(h=20, col="black", lty=2)
points(result_mean_amplitude$amplitude, result_mean_amplitude$fdr_trial*100, col="blue", pch=16)
points(result_mean_amplitude$amplitude, result_mean_amplitude$fdr_trial_plus*100, col="red", pch=16)
points(result_mean_amplitude$amplitude, result_mean_amplitude$fdr_trial_BHq*100, col="black", pch=16)
legend("bottomright", c("Nominal","Knockoff", "Knockoff+", "BHq"), col=c("black","blue", "red", "black"), lwd = c(2.5, 2.5,2.5,2.5),lty=c(2,1,1,1))

plot(result_mean_amplitude$amplitude, result_mean_amplitude$power_knockoff*100, type="l", col="blue",xlab="amplitude level", ylab="Power(%)", ylim=c(0,100))
lines(result_mean_amplitude$amplitude, result_mean_amplitude$power_knockoff_plus*100, col="red")
lines(result_mean_amplitude$amplitude, result_mean_amplitude$power_BHq*100, col="black")
points(result_mean_amplitude$amplitude, result_mean_amplitude$power_knockoff*100, col="blue", pch=16)
points(result_mean_amplitude$amplitude, result_mean_amplitude$power_knockoff_plus*100, col="red", pch=16)
points(result_mean_amplitude$amplitude, result_mean_amplitude$power_BHq*100, col="black", pch=16)
legend("bottomright", c("Knockoff", "Knockoff+", "BHq"), col=c("blue", "red", "black"), lwd = c(2.5,2.5,2.5))

dev.off()

pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Slides Tex/Img/", 
                 "correlation", ".pdf", sep=""), width=10,height=5)
par(mfrow= c(1,2))
plot(result_mean_correlation$correlation, result_mean_correlation$fdr_trial*100, type="l", col="blue",xlab="correlation level", ylab="FDR(%)", ylim=c(0,30))
lines(result_mean_correlation$correlation, result_mean_correlation$fdr_trial_plus*100, col="red")
lines(result_mean_correlation$correlation, result_mean_correlation$fdr_trial_BHq*100, col="black")
abline(h=20, col="black", lty=2)
points(result_mean_correlation$correlation, result_mean_correlation$fdr_trial*100, col="blue", pch=16)
points(result_mean_correlation$correlation, result_mean_correlation$fdr_trial_plus*100, col="red", pch=16)
points(result_mean_correlation$correlation, result_mean_correlation$fdr_trial_BHq*100, col="black", pch=16)
legend("topleft", c("Nominal","Knockoff", "Knockoff+", "BHq"), col=c("black","blue", "red", "black"), lwd = c(2.5, 2.5,2.5,2.5),lty=c(2,1,1,1))

plot(result_mean_correlation$correlation, result_mean_correlation$power_knockoff*100, type="l", col="blue",xlab="correlation level", ylab="Power(%)", ylim=c(0,100))
lines(result_mean_correlation$correlation, result_mean_correlation$power_knockoff_plus*100, col="red")
lines(result_mean_correlation$correlation, result_mean_correlation$power_BHq*100, col="black")
points(result_mean_correlation$correlation, result_mean_correlation$power_knockoff*100, col="blue", pch=16)
points(result_mean_correlation$correlation, result_mean_correlation$power_knockoff_plus*100, col="red", pch=16)
points(result_mean_correlation$correlation, result_mean_correlation$power_BHq*100, col="black", pch=16)
legend("topleft", c("Knockoff", "Knockoff+", "BHq"), col=c("blue", "red", "black"), lwd = c(2.5,2.5,2.5))

dev.off()



pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Slides Tex/Img/", 
                 "sparsity", ".pdf", sep=""), width=10,height=5)
par(mfrow= c(1,2))
plot(result_mean_sparsity$sparsity, result_mean_sparsity$fdr_trial*100, type="l", col="blue",xlab="sparsity level", ylab="FDR(%)", ylim=c(0,30))
lines(result_mean_sparsity$sparsity, result_mean_sparsity$fdr_trial_plus*100, col="red")
lines(result_mean_sparsity$sparsity, result_mean_sparsity$fdr_trial_BHq*100, col="black")
abline(h=20, col="black", lty=2)
points(result_mean_sparsity$sparsity, result_mean_sparsity$fdr_trial*100, col="blue", pch=16)
points(result_mean_sparsity$sparsity, result_mean_sparsity$fdr_trial_plus*100, col="red", pch=16)
points(result_mean_sparsity$sparsity, result_mean_sparsity$fdr_trial_BHq*100, col="black", pch=16)
legend("bottomright", c("Nominal","Knockoff", "Knockoff+", "BHq"), col=c("black","blue", "red", "black"), lwd = c(2.5, 2.5,2.5,2.5),lty=c(2,1,1,1))

plot(result_mean_sparsity$sparsity, result_mean_sparsity$power_knockoff*100, type="l", col="blue",xlab="sparsity level", ylab="Power(%)", ylim=c(0,100))
lines(result_mean_sparsity$sparsity, result_mean_sparsity$power_knockoff_plus*100, col="red")
lines(result_mean_sparsity$sparsity, result_mean_sparsity$power_BHq*100, col="black")
points(result_mean_sparsity$sparsity, result_mean_sparsity$power_knockoff*100, col="blue", pch=16)
points(result_mean_sparsity$sparsity, result_mean_sparsity$power_knockoff_plus*100, col="red", pch=16)
points(result_mean_sparsity$sparsity, result_mean_sparsity$power_BHq*100, col="black", pch=16)
legend("bottomright", c("Knockoff", "Knockoff+", "BHq"), col=c("blue", "red", "black"), lwd = c(2.5,2.5,2.5))

dev.off()

temp <- which(is.na(result_sparsity$fdr_trial))
for(i in temp)
{
  result<- Trials_comparison(nobs = 3000, npar = 1000, rho = 0, nsupp = result_sparsity[i,1])
  result_sparsity[i,2:7] <- result
  if(i%%10==0){cat("Trial ",i," finished. ",as.character(Sys.time()), "\n", sep="")}
}
