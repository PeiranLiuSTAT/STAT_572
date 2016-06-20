##################
# DRMcv function #
##################

# author: Haley Hedlin
# email: hedlin@stanford.edu
# last updated: September 24, 2014

# function to run cross-validated ordinary least squares
# using data in PR_PhenoSense_DataSet.csv AFTER it has been cleaned 
# to remove duplicate patients, mixtures of major DRMs, etc.

# requires R package 'caret'
# lars option requires R package 'glmnet'


###

# To run this function:
# 1. Open R (Install the caret R package if not already installed)
# 2. Open the DRMcv.R file and run the code contained within.
# 3. Type DRMcv() at the R prompt
# 4. Two files will be output into your R working directory
# (Type getwd() at the prompt to get the path)


### arguments to DRMcv

# dataset specifies which dataset the function is to be run on
# must be either "PR", "NRTI", or "NNRTI"

# muts.in is a character vector of the mutations to be included as
# independent variables in the OLS
# each entry in the vector is 3 characters long
  # a letter representing the mutation in the first character
  # followed by two numbers representing the position

# drug is a string indicating the drug to be used for the dependent variable
# must be equal to one of the following: "FPV", "ATV", "IDV", "LPV", "NFV", 
# "RTV", "SQV", "TPV", "DRV"

# min.muts is the minimum number of sequences that a mutation must appear in.
# If a mutation appears in too few sequences, it is removed from the model.

# nfold is the number of folds in the cross-validatioin (CV)
# nrep is the number of times to repeat the CV

# confusion controls whether the confusion matrix should be output

# lars controls whether LARS estimates are output


### output from DRMcv.cleanin

# DRMcv outputs up to four files 
# the files will write to R's working directory
# type getwd() in R to see the working directory path

# the first file contains a matrix containing the OLS estimates
# the estimated coefficients for the input mutations are in the first column
# and the SEs of the estimates are in the second column

# the second file contains the mean square errors (MSEs) estimated from the CV folds

# if confusion=TRUE, the third file will contain a confusion matrix 
# with the actual class in the rows and the predicted class in the columns

# if lars=TRUE, the last file will contain the coefficients from LARS 


DRMcv <- function(dataset="PI", drug="LPV", min.muts=10, nfold=5, nrep=10,
                  muts.in=c("47A", "84A", "50V", "76V", "82A", "82F", "84V",
                            "84C", "82S", "82T", "82M", "32I", "47V", "54M", "54L",
                            "54V", "90M", "54A", "54S", "54T", "46I", "46L", "48V",
                            "48M", "24I", "82C", "33F", "10F", "73S", "73T", "73C",
                            "73A", "11I", "11L", "89V", "20T", "53L", "88S", "50L",
                            "24F", "30N", "43T", "46V", "58E", "83D", "88T", "85V"), 
                       confusion=FALSE,lars=FALSE){
  
  require(caret)
  
  # check that arguments are entered correctly
  if(!is.character(muts.in)){
    stop('The argument "muts.in" must be a character vector.')
  }
  if(!is.character(drug) | length(drug)!=1){
    stop('The argument "drug" must be a character vector of length 1.')
  }
  
  if(!is.character(dataset) | !dataset%in%c("PI","NRTI","NNRTI")){
    stop('The argument "dataset" must be a character vector equal to "PI", 
         "NRTI", or "NNRTI".')
  }
  
  if(drug=="3TC") drug <- "X3TC"
  
  muts.in.conv <- convert.muts(muts.in)
  check.muts(muts.in.conv)
  
  # get the amino acids and positions for the mutations to be included in the model
  mut <- ifelse(nchar(muts.in)==3,toupper(substr(muts.in,3,3)),
                toupper(substr(muts.in,4,4)))
  ps <- suppressWarnings(ifelse(nchar(muts.in)==3,as.numeric(substr(muts.in,1,2)),
                                as.numeric(substr(muts.in,1,3))))  
  
  ## automatically read in the data using url
  if(dataset=="PI"){ 
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt",
                      header=TRUE, sep="\t",
                      stringsAsFactors=FALSE)
    dat[dat=="."] <- NA
    posu <- dat[,10:108]
  }
  if(dataset=="NRTI"){
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)    
    dat[dat=="."] <- NA
    posu <- dat[,8:247]  
  }
  if(dataset=="NNRTI"){ 
    dat <- read.table("http://hivdb.stanford.edu/download/GenoPhenoDatasets/NNRTI_DataSet.txt",
                      header=TRUE, sep="\t", comment.char="@",
                      stringsAsFactors=FALSE)       
    dat[dat=="."] <- NA
    posu <- dat[,6:245]  
  }
  
  # construct design matrix for OLS
  X <- buildX(posu, mut, ps)
  
  # construct dependent variable
  drugcol <- which(names(dat)==drug)    
  Y <- as.numeric(dat[,drugcol])  # absolute measure
  Ylog10 <- log10(Y)
  df.log <- data.frame(Y=Ylog10, X=X)
  
  # remove all rows with missing values
  rem.rows <- unique(which(is.na(df.log),arr.ind=TRUE)[,1])
  df.log.cc <- df.log[-rem.rows,]  # complete case
  
  # remove mutations that are rare
  rare.muts <- which(colSums(df.log.cc[,-1])<min.muts)
  if(length(rare.muts)>0){
    message(paste0(muts.in[rare.muts],
                   " excluded from the model because it appears in fewer than ",
                   min.muts," sequences.\n"))
    df.log.cc <- df.log.cc[,-(rare.muts+1)]  
  }
  
  ### fit the models 
  
  fit <- lm(Y~., data=df.log.cc)
  
  CVout <- train(df.log.cc[,-1], df.log.cc$Y, method='lm', 
                 trControl=trainControl(number=nfold,repeats=nrep))
  
  if(confusion==TRUE){
    ### code for the confusion matrix
    cutoffmat <- matrix(NA, nrow=18, ncol=2)
    rownames(cutoffmat) <- c("FPV","ATV","IDV","LPV","NFV","SQV","TPV","DRV",
                             "X3TC","ABC","AZT","D4T","DDI","TDF",
                             "EFV","NVP","ETR","RPV")
    colnames(cutoffmat) <- c("lower","upper")
    cutoffmat[1,] <- c(3,15) # FPV
    cutoffmat[2,] <- c(3,15) # ATV
    cutoffmat[3,] <- c(3,15) # IDV
    cutoffmat[4,] <- c(9,55) # LPV
    cutoffmat[5,] <- c(3,6) # NFV
    cutoffmat[6,] <- c(3,15) # SQV
    cutoffmat[7,] <- c(2,8) # TPV
    cutoffmat[8,] <- c(10,90) # DRV
    cutoffmat[9,] <- c(5,25) # X3TC
    cutoffmat[10,] <- c(2,6) # ABC
    cutoffmat[11,] <- c(3,15) # AZT
    cutoffmat[12,] <- c(1.5,3) # D4T
    cutoffmat[13,] <- c(1.5,3) # DDI
    cutoffmat[14,] <- c(1.5,3) # TDF
    cutoffmat[15,] <- c(3,10) # EFV
    cutoffmat[16,] <- c(3,10) # NVP
    cutoffmat[17,] <- c(3,10) # ETR
    cutoffmat[18,] <- c(3,10) # RPV
    
    cutoff <- cutoffmat[which(rownames(cutoffmat)==drug),]
    
    # predicted and actual categories
    predicted <- cut(10^predict(fit),c(0,cutoff,Inf),labels=FALSE)
    actual <- cut(10^df.log.cc$Y,c(0,cutoff,Inf),labels=FALSE)
    
    conftab <- table(predicted,actual)
    rownames(conftab) <- colnames(conftab) <- c("susceptible",
                                                "intermediate-level resistant",
                                                "high-level resistant")
  }
  
  
  if(lars==TRUE){
    require(glmnet)
    
    larsfit <- cv.glmnet(as.matrix(df.log.cc[,-1]), df.log.cc$Y, nfolds=5)
    larscoef <- coef(larsfit,s="lambda.min")
  }
  
  # output model coefficients and SEs, and the MSE  
  write.table(summary(fit)$coefficients[,1:2],file="OLScoefs.txt") 
  write.table(CVout$resample$RMSE^2,file="CVmse.txt")
  if(confusion==TRUE) write.table(conftab,file="Confusion.txt")
  if(lars==TRUE) write.table(cbind(rownames(larscoef),matrix(larscoef)),file="LARScoef.txt") 
}


#### helper functions

# function to check that the mutations have been entered correctly
# if the letters are entered as lower case, they are converted to upper case

# convert insertions or deletions to # or ~
convert.muts <- function(muts.in){
  muts.in5 <- which(nchar(muts.in) == 5)
  muts.in6 <- which(nchar(muts.in) == 6)
  for(mi in muts.in5){
    postmp <- substr(muts.in[mi],1,2)
    if(substr(muts.in[mi],3,5) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],3,5) == "del") muts.in[mi] <- paste0(postmp,"~")   
  }
  for(mi in muts.in6){
    postmp <- substr(muts.in[mi],1,3)
    if(substr(muts.in[mi],4,6) == "ins") muts.in[mi] <- paste0(postmp,"#")
    if(substr(muts.in[mi],4,6) == "del") muts.in[mi] <- paste0(postmp,"~")   
  }
  return(muts.in)
}

check.muts <- function(muts.in){
  # all entries should be nchar=3
  if(any(nchar(muts.in) < 3) | any(nchar(muts.in) > 4))
    stop('All entries in argument "muts.in" should be between 3 and 4 characters long.')
  
  muts.in3 <- muts.in[nchar(muts.in) == 3]
  muts.in4 <- muts.in[nchar(muts.in) == 4]
  
  # all should have numbers first two or three characters and a letter for the last 
  if(!all(toupper(substr(muts.in3,3,3))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last 
         character.')
  if(any(is.na(as.numeric(substr(muts.in3,1,2)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
  
  if(!all(toupper(substr(muts.in4,4,4))%in%c(LETTERS,"#","~")))
    stop('All entries in argument "muts.in" must have a letter, #, or ~ in the last 
         character.')
  if(any(is.na(as.numeric(substr(muts.in4,1,3)))))
    stop('All entries in argument "muts.in" must begin in two or three digits.')
}



# function to create the design matrix X with the input mutations/positions
buildX <- function(dat, mut, ps){
  X <- matrix(NA, nrow=nrow(dat), ncol=length(mut))

  # loop through all positions
  for(p in unique(ps)){
    p1 <- substr(dat[,p],1,1)  # first mutation at this position
    p2 <- substr(dat[,p],2,2)
    for(ind in which(ps==p)){
      X[,ind] <-  as.numeric(p1==as.character(mut[ind]) | 
                               p2==as.character(mut[ind]))  
    }
  }  
  colnames(X) <- paste0(ps,mut)  
  return(X)
}


