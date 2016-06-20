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
library(PerformanceAnalytics)
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

Tableconstruction <- function(filename= "PI")
{
  fileread <- paste(filename, "_DATA.csv", sep="")
  data <- read.csv(fileread, header=T, stringsAsFactors = F)
  x_start <- which(names(data)=="P1")
  pos_index <- x_start:ncol(data)
  if(filename == "PI") 
  {
    data <- data[-c(600,621),]
  }
  else if(filename == "NRTI")
  {
    data <- data[-c(56,148,237,472,480),]
  }
  else 
  {
    data <- data[-c(558, 565),]
  }
  
  
  allmutations <- NULL
  for(i in pos_index)
  {
    maxmut <- max(sapply(data[,i], nchar))
    for(j in 1:maxmut)
    {
      allmutations <- c(allmutations, sapply(data[,i], substr, start= j, stop=j))
    }
    allmutations <- unique(allmutations)
  }
  
  
  ### Construct X
  X <- data.frame(id=1:nrow(data))
  
  columns <- 2
  for(i in pos_index)
  {
    temp <- data[,i]
    maxmut <- max(nchar(temp))
    allmutations <- NULL
    for(j in 1:maxmut)
    {
      allmutations <- c(allmutations, sapply(data[,i], substr, start= j, stop=j))
    }
    allmutations <- unique(allmutations)
    allmutations <- allmutations[which(allmutations!="-")]
    allmutations <- allmutations[which(allmutations!="*")]
    allmutations <- allmutations[which(allmutations!=".")]
    allmutations <- allmutations[which(allmutations!="")]
    if(length(allmutations)==0){next}
    for(j in 1:length(allmutations))
    {
      Xj <- numeric(nrow(data))
      for(l in 1:maxmut)
      {
        Xj <- Xj + (substr(data[,i],start=l,stop=l)==allmutations[j])
      }
      X <- cbind(X,Xj)
      names(X)[columns] <- paste("P", i-x_start+1,"_", allmutations[j], sep="")
      columns <- columns + 1
    }
  }
  X <- X[,-1]
  indicesremove <- which(colSums(X)<3)
  X <- X[,-indicesremove]
   
  ### Construct Y ####
  Y <- data[,4:(x_start-1)]
  return(list(namedata = filename, X=X, Y=Y))
}

PI_DATA <- Tableconstruction()
NRTI_DATA <- Tableconstruction(filename = "NRTI")
NNRTI_DATA <- Tableconstruction(filename = "NNRTI")

thres_compute <- function(W, q=0.2, type="knockoff")
{
  temp <- sort(unique(abs(W)))[-1]
  k <- length(temp)
  thres <- Inf
  for(i in 1:k)
  {
    numerator <- sum(W <= -temp[i])
    denumerator <- max(sum(W >= temp[i]),1)
    #browser()
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


analysis <- function(dataset = PI_DATA, q =0.2)
{
  responses <-dataset$Y
  features <- dataset$X
  numofdrugs <- ncol(responses)
  result_list <- list()
  for(k in 1:numofdrugs)
  {
    ## Begin Analysis
    Yk <- log10(responses[,k])
    ## NA Remove
    featuresk <- features[!is.na(Yk),]
    Yk <- Yk[!is.na(Yk)]
    featuresk <- featuresk[!duplicated(as.list(featuresk))]
    X <- featuresk[,colSums(featuresk)>=3]
    n <- nrow(X)
    if(n < 2*ncol(X))
    {
      #browser()
      npar <- ncol(X)
      X.svd <- svd(X, nu=n, nv=0)
      u2 <- X.svd$u[,(npar+1):n]
      sigma <- sqrt(mean((t(u2)%*%Yk)^2))
      set.seed(0)
      y.extra <- rnorm(2*npar-n, sd=sigma)
      X <- rbind(as.matrix(X), matrix(0,2*npar-n,npar))
      X <- as.data.frame(X)
      Yk <- c(Yk, y.extra)
    }
    
    #X <- apply(X, 2, scale, center=F)/sqrt(nrow(X)-1)
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
    model <- glmnet(cbind(as.matrix(X), tildeX_equiv), Yk, intercept=F, standardize=F,
                    standardize.response=F)
    lambda_max <- max(model$lambda)
    lambda_min <- lambda_max/2e3
    ks <- (0:(5*npar-1))/(5*npar)
    lambdas <- lambda_max * (lambda_min/lambda_max)^ks
    model_equiv <- glmnet(cbind(as.matrix(X), tildeX_equiv), Yk, lambda = lambdas, intercept=F, standardize=F,
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
    Thres <- thres_compute(W_equiv, type="knockoff", q = q)
    index_equiv <- which(W_equiv>= Thres)
    Thres <- thres_compute(W_equiv, type="knockoff+", q=q)
    index_equiv_plus <- which(W_equiv>= Thres)
    
    model_ls <- lm(Yk~X+0)
    Z_ls <- summary(model_ls)$coef[, "t value"]
    Thres <- thres_compute(Z_ls, type="BHq")
    index.BHq <- which(abs(Z_ls) >= Thres)
    
    #packageresult <- knockoff.filter(X1, Yk)
    #browser()
    
    result_selected <- list(index_knockoff = index_equiv, index_knockoff_plus = index_equiv_plus, index_BHq = index.BHq,
                            names_knockoff = colnames(X)[index_equiv], names_BHq = colnames(X)[index.BHq], ncolumn = ncol(X), nrows=n)
    result_list[[k]] <- result_selected
  }
  names(result_list) <- names(responses)
  return(result_list)
}

analysis_paper <- function(dataset = PI_DATA, q =0.2)
{
  responses <-dataset$Y
  features <- dataset$X
  numofdrugs <- ncol(responses)
  result_list <- list()
  for(k in 1:numofdrugs)
  {
    ## Begin Analysis
    Yk <- log10(responses[,k])
    ## NA Remove
    featuresk <- features[!is.na(Yk),]
    Yk <- Yk[!is.na(Yk)]
    X <- featuresk[,colSums(featuresk)>=3]
    duplicates <- which(colSums((cor(X)-1)^2 <1e-8)>1)
    if(length(duplicates) > 0)
    {
      X <- X[,-duplicates]
    }
    
    n <- nrow(X)
    if(n < 2*ncol(X))
    {
      #browser()
      npar <- ncol(X)
      X.svd <- svd(X, nu=n, nv=0)
      u2 <- X.svd$u[,(npar+1):n]
      sigma <- sqrt(mean((t(u2)%*%Yk)^2))
      set.seed(0)
      y.extra <- rnorm(2*npar-n, sd=sigma)
      X <- rbind(as.matrix(X), matrix(0,2*npar-n,npar))
      X <- as.data.frame(X)
      Yk <- c(Yk, y.extra)
    }
    
    #X <- apply(X, 2, scale, center=F)/sqrt(nrow(X)-1)
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
    model <- glmnet(cbind(as.matrix(X), tildeX_equiv), Yk, intercept=F, standardize=F,
                    standardize.response=F)
    lambda_max <- max(model$lambda)
    lambda_min <- lambda_max/2e3
    ks <- (0:(5*npar-1))/(5*npar)
    lambdas <- lambda_max * (lambda_min/lambda_max)^ks
    model_equiv <- glmnet(cbind(as.matrix(X), tildeX_equiv), Yk, lambda = lambdas, intercept=F, standardize=F,
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
    Thres <- thres_compute(W_equiv, type="knockoff", q = q)
    index_equiv <- which(W_equiv>= Thres)
    Thres <- thres_compute(W_equiv, type="knockoff+", q=q)
    index_equiv_plus <- which(W_equiv>= Thres)
    
    model_ls <- lm(Yk~X+0)
    Z_ls <- summary(model_ls)$coef[, "t value"]
    Thres <- thres_compute(Z_ls, type="BHq")
    index.BHq <- which(abs(Z_ls) >= Thres)
    
    #packageresult <- knockoff.filter(X1, Yk)
    #browser()
    
    result_selected <- list(index_knockoff = index_equiv, index_knockoff_plus = index_equiv_plus, index_BHq = index.BHq,
                            names_knockoff = colnames(X)[index_equiv], names_BHq = colnames(X)[index.BHq], ncolumn = ncol(X), nrows = n)
    result_list[[k]] <- result_selected
  }
  names(result_list) <- names(responses)
  return(result_list)
}

result_PI <- analysis()
result_NRTI <- analysis(dataset = NRTI_DATA)
result_NNRTI <- analysis(dataset = NNRTI_DATA)

result_PI_paper <- analysis_paper()
result_NRTI_paper <- analysis_paper(dataset = NRTI_DATA)
result_NNRTI_paper <- analysis_paper(dataset = NNRTI_DATA)


tsm <- function(filename = "PI")
{
  fileread <- paste(filename, ".txt", sep="")
  tsm_uncleaned <- readLines(fileread)
  tsm <- numeric(length(tsm_uncleaned))
  for(i in 1:length(tsm_uncleaned))
  {
    temp <- strsplit(tsm_uncleaned[i],"")
    temp_end <- which(temp[[1]] =="\t")-1
    temp <- paste(temp[[1]][1:(temp_end)], collapse="")
    tsm[i] <- as.numeric(temp)
  }
  return(tsm)
}

tsm_PI <- tsm("PI")
tsm_NRTI <- tsm("NRTI")
tsm_NNRTI <- tsm("NNRTI")

pickselected <- function(result = result_PI_paper, tsm = tsm_PI)
{
  selected_ko <- list(length(result))
  selected_BHq <- list(length(result))
  discoveries_ko <- numeric(length(result))
  discoveries_BHq <- numeric(length(result))
  falsediscoveries_ko <- numeric(length(result))
  falsediscoveries_BHq <- numeric(length(result))
  for(i in 1:length(result))
  {
    starts <- 2
    #browser()
    stops <- length(result[[i]]$names_knockoff)
    for(j in 1:length(result[[i]]$names_knockoff))
    {
      stops[j] <- which(strsplit(result[[i]]$names_knockoff[j], split="")[[1]] == "_")
    }
    selected_pos <- as.numeric(substr(result[[i]]$names_knockoff,start = 2, stop = stops-1))
    selected_ko[[i]] <- unique(selected_pos)
    discoveries_ko[i] <- length(selected_ko[[i]])
    falsediscoveries_ko[i] <- length(setdiff(selected_ko[[i]], tsm))
    
    stops <- length(result[[i]]$names_BHq)
    for(j in 1:length(result[[i]]$names_BHq))
    {
      stops[j] <- which(strsplit(result[[i]]$names_BHq[j], split="")[[1]] == "_")
    }
    selected_pos <- as.numeric(substr(result[[i]]$names_BHq,start = 2, stop = stops-1))
    selected_BHq[[i]] <- unique(selected_pos)
    discoveries_BHq[i] <- length(selected_BHq[[i]])
    falsediscoveries_BHq[i] <- length(setdiff(selected_BHq[[i]], tsm))
  }
  result_data <- data.frame(discoveries_ko = discoveries_ko, discoveries_BHq = discoveries_BHq, 
                            falsediscoveries_ko = falsediscoveries_ko, falsediscoveries_BHq = falsediscoveries_BHq)
  return(list(result_data = result_data, selected_ko = selected_ko, selected_BHq = selected_BHq))
}

selection_PI <- pickselected(result = result_PI)$result_data
selection_PI_paper <- pickselected()$result_data
selection_NRTI <- pickselected(result = result_NRTI, tsm=tsm_NRTI)$result_data
selection_NRTI_paper <- pickselected(result = result_NRTI_paper, tsm=tsm_NRTI)$result_data
selection_NNRTI <- pickselected(result = result_NNRTI, tsm=tsm_NNRTI)$result_data
selection_NNRTI_paper <- pickselected(result = result_NNRTI_paper, tsm=tsm_NNRTI)$result_data

simulation_nonGaussian <- function(q=0.2)
{
  responses <-NNRTI_DATA$Y
  features <- NNRTI_DATA$X
  z <- NULL
  for(i in 1:3)
  {
    Yk <- log(responses[,i])
    ## NA Remove
    featuresk <- features[!is.na(Yk),]
    Yk <- Yk[!is.na(Yk)]
    duplicates <- which(colSums((cor(featuresk)-1)^2 <1e-8)>1)
    if(length(duplicates) > 0)
    {
      Xfeatures <- featuresk[,-duplicates]
    }
    else
    {
      Xfeatures <- featuresk
    }
    model <- lm(Yk~as.matrix(Xfeatures)+0)
    z <- c(z,model$residuals)
  }
  z <- scale(z,center=F)
  #z <- z/sqrt(mean(z^2))
  X <- NNRTI_DATA$X
  eps <- sample(z, size = nrow(X), replace=TRUE)
  #duplicates <- which(colSums((cor(X)-1)^2 <1e-8)>1)
  #if(length(duplicates) > 0)
  #{
  #  X <- X[,-duplicates]
  #}
  X <- X[,!duplicated(as.list(as.data.frame(X)))]
  npar <- ncol(X)
  truesupp <- sample(npar, size=20, replace=FALSE)
  betatrue <- numeric(npar)
  betatrue[truesupp] <- rnorm(20, 0, 1) * 3.5
  Y <- as.matrix(X)[,truesupp] %*% betatrue[truesupp] + eps
  
  #npar <- ncol(X)
  X <- normc(X)
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
    temp <- which(abs(model_equiv$beta[i,])>0)
    if(length(temp) == 0) {Z_equiv[i] <- 0}
    else
    {
      temp <- min(temp)
      Z_equiv[i] <- model_equiv$lambda[temp]
    }
  }
  Z_equiv <- matrix(Z_equiv, ncol=2)
  W_equiv <- apply(Z_equiv, 1, max) * sign(Z_equiv[,1]- Z_equiv[,2])
  Thres <- thres_compute(W_equiv, type="knockoff", q = q)
  index_equiv <- which(W_equiv>= Thres)
  Thres <- thres_compute(W_equiv, type="knockoff+", q=q)
  index_equiv_plus <- which(W_equiv>= Thres)
  
  model_ls <- lm(Y~as.matrix(X)+0)
  Z_ls <- summary(model_ls)$coef[, "t value"]
  Thres <- thres_compute(Z_ls, type="BHq")
  index.BHq <- which(abs(Z_ls) >= Thres)
  fdr_trial <- sum(!(index_equiv %in% truesupp))/max(length(index_equiv),1)
  power_knockoff <- sum(index_equiv %in% truesupp)/20
  fdr_trial_plus<- sum(!(index_equiv_plus %in% truesupp))/max(length(index_equiv_plus),1)
  power_knockoff_plus <- sum(index_equiv_plus %in% truesupp)/20
  fdr_trial_BHq <- sum(!(index.BHq %in% truesupp))/max(length(index.BHq),1)
  power_BHq <- sum(index.BHq %in% truesupp)/20
  
  return(c(fdr_trial, fdr_trial_plus, fdr_trial_BHq, 
           power_knockoff, power_knockoff_plus, power_BHq))
}

result <- matrix(nrow =1000, ncol=6)
result <- as.data.frame(result)
names(result) <- c("fdr_trial", "fdr_trial_plus", "fdr_trial_BHq",
                               "power_knockoff","power_knockoff_plus", "power_BHq")
for(i in 1:1000)
{
  result[i,] <- simulation_nonGaussian()
  gc()
  if(i%%10==0){cat("Trial ",i," finished. ",as.character(Sys.time()), "\n", sep="")}
}



## Figure 7
counts <- t(rbind(cbind(selection_PI$discoveries_ko, selection_PI$falsediscoveries_ko), cbind(selection_PI$discoveries_BHq, selection_PI$falsediscoveries_BHq)))
counts[1,] <- counts[1,] - counts[2,]
counts <- as.data.frame(counts)
row.names(counts) <- c("In TSM", "Not in TSM")
names(counts) <- c(rep("Knockoff", 7), rep("BHq", 7))
pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "PI", ".pdf", sep=""), width=10,height=10)
par(mfrow=c(3,3))
par(mar = c(2,3,2,2)+0.1)
for(i in 1:7)
{
  if(i == 1)
  {
    barplot(as.matrix(counts[,c(i,i+7)]), ylim= c(0,40), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(PI_DATA$Y)[i]), border = "black", legend=c("In TSM", "Not in TSM"), 
            args.legend = list(x = 2.5, y = 33, cex=1.5))
  }
  else
  {
    barplot(as.matrix(counts[,c(i,i+7)]), ylim= c(0,40), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5,cex.names = 1.5,
            main = paste("resistance to",names(PI_DATA$Y)[i]), border = "black")
  }
  box()
  abline(h = length(tsm_PI), lty=2)
  text(x=1.6, y=length(tsm_PI)+2, paste("n=", result_PI[[i]]$nrows, ", p=", result_PI[[i]]$ncolumn), cex=2)
}
dev.off()

counts <- t(rbind(cbind(selection_PI_paper$discoveries_ko, selection_PI_paper$falsediscoveries_ko), cbind(selection_PI_paper$discoveries_BHq, selection_PI_paper$falsediscoveries_BHq)))
counts[1,] <- counts[1,] - counts[2,]
counts <- as.data.frame(counts)
row.names(counts) <- c("In TSM", "Not in TSM")
names(counts) <- c(rep("Knockoff", 7), rep("BHq", 7))
pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "PI_paper", ".pdf", sep=""), width=10,height=10)
par(mfrow=c(3,3))
par(mar = c(2,3,2,2)+0.1)
for(i in 1:7)
{
  if(i == 1)
  {
    barplot(as.matrix(counts[,c(i,i+7)]), ylim= c(0,40), col=c("darkblue","red"),  cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(PI_DATA$Y)[i]), border = "black", legend=c("In TSM", "Not in TSM"), 
            args.legend = list(x = 2.5, y = 33, cex=1.5))
  }
  else
  {
    barplot(as.matrix(counts[,c(i,i+7)]), ylim= c(0,40), col=c("darkblue","red"),  cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(PI_DATA$Y)[i]), border = "black")
  }
  box()
  abline(h = length(tsm_PI), lty=2)
  text(x=1.5, y=length(tsm_PI)+2, paste("n=", result_PI_paper[[i]]$nrows, ", p=", result_PI_paper[[i]]$ncolumn), cex=2)
}
dev.off()

counts <- t(rbind(cbind(selection_NRTI$discoveries_ko, selection_NRTI$falsediscoveries_ko), cbind(selection_NRTI$discoveries_BHq, selection_NRTI$falsediscoveries_BHq)))
counts[1,] <- counts[1,] - counts[2,]
counts <- as.data.frame(counts)
row.names(counts) <- c("In TSM", "Not in TSM")
names(counts) <- c(rep("Knockoff", 6), rep("BHq", 6))
pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "NRTI", ".pdf", sep=""), width=10,height=6.7)
par(mfrow=c(2,3))
par(mar = c(2,3,2,2)+0.1)
for(i in 1:6)
{
  if(i == 1)
  {
    barplot(as.matrix(counts[,c(i,i+6)]), ylim= c(0,31), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NRTI_DATA$Y)[i]), border = "black", legend=c("In TSM", "Not in TSM"), 
            args.legend = list(x = 2.5, y = 24, cex=1.5))
  }
  else
  {
    barplot(as.matrix(counts[,c(i,i+6)]), ylim= c(0,31), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NRTI_DATA$Y)[i]), border = "black")
  }
  box()
  abline(h = length(tsm_NRTI), lty=2)
  text(x=1, y=length(tsm_NRTI)+2, paste("n=", result_NRTI[[i]]$nrows, ", p=", result_NRTI[[i]]$ncolumn), cex=2)
}
dev.off()

counts <- t(rbind(cbind(selection_NRTI_paper$discoveries_ko, selection_NRTI_paper$falsediscoveries_ko), cbind(selection_NRTI_paper$discoveries_BHq, selection_NRTI_paper$falsediscoveries_BHq)))
counts[1,] <- counts[1,] - counts[2,]
counts <- as.data.frame(counts)
row.names(counts) <- c("In TSM", "Not in TSM")
names(counts) <- c(rep("Knockoff", 6), rep("BHq", 6))
pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "NRTI_paper", ".pdf", sep=""), width=10,height=6.7)
par(mfrow=c(2,3))
par(mar = c(2,3,2,2)+0.1)
for(i in 1:6)
{
  if(i == 1)
  {
    barplot(as.matrix(counts[,c(i,i+6)]), ylim= c(0,31), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NRTI_DATA$Y)[i]), border = "black", legend=c("In TSM", "Not in TSM"), 
            args.legend = list(x = 2.5, y = 24, cex=1.5))
  }
  else
  {
    barplot(as.matrix(counts[,c(i,i+6)]), ylim= c(0,31), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NRTI_DATA$Y)[i]), border = "black")
  }
  box()
  abline(h = length(tsm_NRTI), lty=2)
  text(x=1, y=length(tsm_NRTI)+2, paste("n=", result_NRTI_paper[[i]]$nrows, ", p=", result_NRTI_paper[[i]]$ncolumn), cex=2)
}
dev.off()

counts <- t(rbind(cbind(selection_NNRTI$discoveries_ko, selection_NNRTI$falsediscoveries_ko), cbind(selection_NNRTI$discoveries_BHq, selection_NNRTI$falsediscoveries_BHq)))
counts[1,] <- counts[1,] - counts[2,]
counts <- as.data.frame(counts)
row.names(counts) <- c("In TSM", "Not in TSM")
names(counts) <- c(rep("Knockoff", 3), rep("BHq", 3))
pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "NNRTI", ".pdf", sep=""), width=10,height=3.3)
par(mfrow=c(1,3))
par(mar = c(2,3,2,2)+0.1)
for(i in 1:3)
{
  if(i == 1)
  {
    barplot(as.matrix(counts[,c(i,i+3)]), ylim= c(0,36), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NNRTI_DATA$Y)[i]), border = "black", legend=c("In TSM", "Not in TSM"), 
            args.legend = list(x="topleft", cex=1.5))
  }
  else
  {
    barplot(as.matrix(counts[,c(i,i+3)]), ylim= c(0,36), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NNRTI_DATA$Y)[i]), border = "black")
  }
  box()
  abline(h = length(tsm_NNRTI), lty=2)
  text(x=1, y=length(tsm_NNRTI)+11, paste("n=", result_NNRTI[[i]]$nrows, ", p=", result_NNRTI[[i]]$ncolumn), cex=2)
}
dev.off()

counts <- t(rbind(cbind(selection_NNRTI_paper$discoveries_ko, selection_NNRTI_paper$falsediscoveries_ko), cbind(selection_NNRTI_paper$discoveries_BHq, selection_NNRTI_paper$falsediscoveries_BHq)))
counts[1,] <- counts[1,] - counts[2,]
counts <- as.data.frame(counts)
row.names(counts) <- c("In TSM", "Not in TSM")
names(counts) <- c(rep("Knockoff", 3), rep("BHq", 3))
pdf(file = paste("C:/Users/Peiran Liu/Documents/UW courses/Prelim/Tex Format/Img/", 
                 "NNRTI_paper", ".pdf", sep=""), width=10,height=3.3)
par(mfrow=c(1,3))
par(mar = c(2,3,2,2)+0.1)
for(i in 1:3)
{
  if(i == 1)
  {
    barplot(as.matrix(counts[,c(i,i+3)]), ylim= c(0,36), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NNRTI_DATA$Y)[i]), border = "black", legend=c("In TSM", "Not in TSM"), 
            args.legend = list(x="topleft", cex=1.5))
  }
  else
  {
    barplot(as.matrix(counts[,c(i,i+3)]), ylim= c(0,36), col=c("darkblue","red"), cex.main=1.5,cex.lab=1.5, cex.names = 1.5,
            main = paste("resistance to",names(NNRTI_DATA$Y)[i]), border = "black")
  }
  box()
  abline(h = length(tsm_NNRTI), lty=2)
  text(x=1, y=length(tsm_NNRTI)+11, paste("n=", result_NNRTI_paper[[i]]$nrows, ", p=", result_NNRTI_paper[[i]]$ncolumn), cex=2)
}
dev.off()


