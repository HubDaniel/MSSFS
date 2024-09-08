# used for both parametric bootstrap CI and simulations
# (1) for CI, set true values as in real data; same # of subjects, same missing value place for each subject
# (1) for CI, need to load in the real data, i.e., dataY first for missing values

library(data.table)
library(dplyr)
library(tidyr)

`%notin%` <- Negate(`%in%`)
setwd("~/SSM/use_covariates")


library(doSNOW)
library(foreach)
library(parallel)
library(xtable)

myF <- function(s){
  source("ACODE_rcpp_likelihood_GgVaZ1bW_fixRef0_expWeight.R")
  source("ACODE_rcpp_MKFIS_GgVaZ1bW_fixRef0_expWeight.R")
  source("ACODE_utils_GgVaZ1bW_fixRef0_expWeight.R")
  # univariate case with simplest settings
  # F=I ; fixed W0 W1; no beta
  # estimate: V, G0, G1, g0, g1, alpha0, alpha1, zeta0=0, zeta1
  # generate Y to be small value, if large, eta in lh function is large which leads to exp(-x)=0
  print(paste0("work: ",s))
  set.seed(s)
  
  ptm <- proc.time()
  
  HYPERS = hyperSetting()
  
  use_BFGS = HYPERS$use_BFGS
  
  generatedTheta <- vector("list",HYPERS$numS)
  generatedX <- vector("list",HYPERS$numS)
  generatedinit <- vector("list",HYPERS$numS)
  generatedpi <- vector("list",HYPERS$numS)
  generatedI <- vector("list",HYPERS$numS)
  generatedmissingt <- vector("list",HYPERS$numS)
  generatedY_full <- vector("list",HYPERS$numS)
  generatedY <- vector("list",HYPERS$numS)
  generatedI_full <- vector("list",HYPERS$numS)
  
  
  for(kk in 1:HYPERS$numS){
    num_miss = 0#sample(round(1:HYPERS$timeL/3),1)
    generatedData <- gdataus(num_miss)
    generatedY[[kk]] <- generatedData$generatedY
    # generatedY[[kk]][is.na(dataY[[kk]])] <- NA
    generatedX[[kk]] <- generatedData$generatedX
    generatedinit[[kk]] <- initialus(
      par = generatedData$true_par,
      Y = generatedData$generatedY,
      useTruePara = F)
  }
  
  
  generatedpar <- generatedData$true_par
  
  #####################
  ##### algorithm #####
  #####################
  
  ptm <- proc.time()
  
  estlagL1 <- HYPERS$estlagL1
  
  use_BFGS <- HYPERS$use_BFGS
  
  generatedinit <- vector("list",HYPERS$numS)
  
  for(kk in 1:HYPERS$numS){
    generatedinit[[kk]] <- initialus(
      par = NA,
      Y = generatedY[[kk]],
      useTruePara = F,
      HYPERS)
  }
  
  
  
  
  
  
  
  ##### EM algorithm ##### 
  
  ### store each EM iteration results
  Likelihoods <- c()
  thetaDiff <- c() # MSE of estimated theta and real theta
  
  
  #########################
  ### initial EM - no Z ###
  #########################
  ### MLE estimate parameters
  Generatedinit <- Reduce(`+`, generatedinit)/HYPERS$numS
  store = c()
  
  
  if(!use_BFGS){
    ests <- try(optim(par = Generatedinit[1:(length(Generatedinit)-1)],
                      fn=lhcpp,
                      hypers=HYPERS,
                      Y=generatedY,X=generatedX,Z=NULL,
                      init=1,control = list(fnscale = -1,maxit=10000)))
  }else{
    ests <- try( optim(par = Generatedinit[1:(length(Generatedinit)-1)],
                       fn=lhcpp,
                       Y=generatedY,X=generatedX,Z=NULL,
                       init=1,
                       hypers=HYPERS,
                       control = list(fnscale = -1,maxit=10000),
                       method="L-BFGS-B",
                       lower = c(rep(-10,6),rep(-Inf,2+2*HYPERS$xL)),
                       upper = c(rep(10,6),rep(Inf,2+2*HYPERS$xL)) ))
    
  }
  
  ### extract MLE estimates
  MLEs <- MLEus(ests,init=T)
  
  ### MKF and MFIS with estimated parameters
  smoothedThetar <- vector("list",HYPERS$numS)# the most current smooth estimators for all subjects
  filteredThetar <- vector("list",HYPERS$numS)# the most current filtered estimators for all subjects
  posteriorProbf <- vector("list",HYPERS$numS) # the most current filtered posterior probability for all subjects
  posteriorProb <- vector("list",HYPERS$numS) # the most current smoothed posterior probability for all subjects
  tranProb01 <- vector("list",HYPERS$numS)# the most current filtered estimators for all subjects
  tranProb11 <- vector("list",HYPERS$numS)# the most current filtered estimators for all subjects
  
  for(kk in 1:HYPERS$numS){
    temp <- MKFIScpp(MLES = MLEs,
                     Y = generatedY[[kk]],
                     X = generatedX[[kk]],
                     Z = matrix(NA),
                     init = 1,
                     hypers=HYPERS)
    smoothedThetar[[kk]] <- temp$ThetaSmooth
    posteriorProb[[kk]] <- as.numeric(temp$PrIt__T__1__all)
    posteriorProbf[[kk]] <-temp$PrI__t__t[2,]
    filteredThetar[[kk]] <- temp$ThetaForward
  }
  
  ### store results
  Likelihoods <- c(Likelihoods,ests$value)
  thetaDiff <- c(thetaDiff,mean((smoothedThetar[[1]]-generatedTheta[[1]][,-1])^2))
  oldMLEs <- MLEs
  oldMLEs$est_zeta1 <- 0  # add the initial value of zeta1 as the estimate of this first iteration
  
  ##################
  ### rest of EM ###
  ##################
  
  MLEsA = round(as.matrix(unlist(oldMLEs),ncol=1),4)
  numIter=1
  max_rel_change = 1e-3
  relChange = 1
  while( (numIter <= HYPERS$numEM) & (relChange > max_rel_change)){
    
    if(numIter==1){
      
      if(!use_BFGS){
        ests <- try(optim(par = c(ests$par,0),
                          fn=lhcpp,
                          Y=generatedY,X=generatedX,Z=smoothedThetar,
                          init=0,hypers=HYPERS,
                          control = list(fnscale = -1,maxit=10000)))
      }else{
        ests <- try(optim(par = c(ests$par,0),
                          fn=lhcpp,
                          Y=generatedY,X=generatedX,Z=smoothedThetar,
                          init=0,hypers=HYPERS,
                          control = list(fnscale = -1,maxit=10000),
                          method="L-BFGS-B",
                          lower = c(rep(-10,6),rep(-Inf,3+2*HYPERS$xL)),
                          upper = c(rep(10,6),rep(Inf,3+2*HYPERS$xL)) ))
      }
    }else{
      
      if(!use_BFGS){
        ests <- try(optim(par = ests$par,
                          fn=lhcpp,
                          Y=generatedY,X=generatedX,Z=smoothedThetar,
                          init=0,hypers=HYPERS,
                          control = list(fnscale = -1,maxit=10000)))
      }else{
        ests <- try(optim(par = ests$par,
                          fn=lhcpp,
                          Y=generatedY,X=generatedX,Z=smoothedThetar,
                          init=0,hypers=HYPERS,
                          control = list(fnscale = -1,maxit=10000),
                          method="L-BFGS-B",
                          lower = c(rep(-10,6),rep(-Inf,3+2*HYPERS$xL)),
                          upper = c(rep(10,6),rep(Inf,3+2*HYPERS$xL)) ))
      }
    }
    
    
    if(class(ests)=="try-error"){
      print(paste0("work: ",s," estimation error"))
      return(list(error="error",para=oldMLEs,Z=smoothedThetar,generatedY=generatedY,X=generatedX,generatedI=generatedI,generatedTheta=generatedTheta))
    }
    
    MLEs <- MLEus(ests,init=F)
    
    for(kk in 1:HYPERS$numS){
      temp <- MKFIScpp(MLES = MLEs,
                       Y = generatedY[[kk]],
                       X = generatedX[[kk]],
                       Z = smoothedThetar[[kk]],
                       init = 0,
                       hypers=HYPERS)
      
      smoothedThetar[[kk]] <- temp$ThetaSmooth
      posteriorProb[[kk]] <- as.numeric(temp$PrIt__T__1__all)
      posteriorProbf[[kk]] <-temp$PrI__t__t[2,]
      filteredThetar[[kk]] <- temp$ThetaForward
    }
    
    ### store results
    Likelihoods <- c(Likelihoods,ests$value)
    
    newEsts = c(MLEs$est_V,MLEs$est_W0,MLEs$est_W1,
                MLEs$est_delta,MLEs$est_G0,MLEs$est_G1,
                MLEs$est_alpha0,MLEs$est_alpha1,MLEs$est_zeta1,
                MLEs$est_beta0,MLEs$est_beta1)
    oldEsts = c(oldMLEs$est_V,oldMLEs$est_W0,oldMLEs$est_W1,
                oldMLEs$est_delta,oldMLEs$est_G0,oldMLEs$est_G1,
                oldMLEs$est_alpha0,oldMLEs$est_alpha1,oldMLEs$est_zeta1,
                oldMLEs$est_beta0,oldMLEs$est_beta1)
    
    relChange <- sum((newEsts-oldEsts)^2)/(sum(oldEsts^2) + 1e-6)
    
    print(paste0("work: ",s," numIter/max numEM: ",numIter," / ",HYPERS$numEM," relChange: ",relChange))
    
    
    oldMLEs <- MLEs
    
    MLEsA <- cbind(MLEsA,round(unlist(MLEs)[rownames(MLEsA)],4))
    numIter = numIter+1
  }
  
  
  
  MLES = unlist(MLEs)[-which(names(MLEs)%in%c("est_g0","est_g1","est_ref1"))]
  EstimatesR <- matrix(nrow=2,ncol=length(MLES),
                       dimnames = list( c("true","estimated"),names(MLES) ))
  EstimatesR[2,] <- MLES
  EstimatesR[1,] <- unlist(generatedpar)
  EstimatesR <- cbind(EstimatesR,c(3,lags))
  colnames(EstimatesR)[ncol(EstimatesR)] = "est_L"
  ttt <- proc.time() - ptm
  
  EstimatesR <- cbind(EstimatesR,c(relChange,NA),c(ttt[3],NA),
                      c(Likelihoods[length(Likelihoods)],NA),c(thetaDiff[length(thetaDiff)],NA))
  
  colnames(EstimatesR)[(ncol(EstimatesR)-3):ncol(EstimatesR)] <- c("RelChange","time","log-lh","MSE theta")
  
  
  Returns <- list(EstimatesR=EstimatesR)
  Returns$posteriorProb <- posteriorProb
  Returns$posteriorProbf <- posteriorProbf
  Returns$smoothedThetar <- smoothedThetar
  Returns$filteredThetar <- filteredThetar
  Returns$optimR <- ests
  Returns$generatedTheta <- generatedTheta
  Returns$generatedY <- generatedY
  Returns$generatedI <- generatedI
  Returns$generatedY_full <- generatedY_full
  Returns$generatedI_full <- generatedI_full
  Returns$MLEsA = MLEsA
  
  saveRDS(Returns,paste0("d_10_sub_300_N_withL/",s,".rds"))
  
  # for(i in 1:HYPERS$numS){
  #   par(mfrow=c(4,1))
  #   plot(generatedY[[i]][1,],type="b",ylab="O2",main=paste0("subject ",i))
  #   abline(v=which(generatedI[[i]]==1),col="red",lty=3)
  #   plot(generatedpi[[i]],type="b",ylab="transition infection probability")
  #   abline(v=which(generatedI[[i]]==1),col="red",lty=3)
  # 
  #   plot(as.numeric(smoothedThetar[[i]]),type="l",col="red",ylab="smooth est Theta")
  #   lines(as.numeric(generatedTheta[[i]][,-1]),type="l",col="black",lty=3,ylab="real Theta")
  #   lines(as.numeric(filteredThetar[[i]]),type="l",col="blue",lty=3,ylab="filter est Theta")
  #   abline(v=which(generatedI[[i]]==1),col="red",lty=3)
  # 
  #   plot(as.numeric(posteriorProb[[i]]),type="l",col="red",ylab="est posterior prob")
  #   lines(as.numeric(posteriorProbf[[i]]),type="l",col="blue",lty=3,ylab="est posterior prob")
  #   lines(as.numeric(generatedpi[[i]]),type="l",col="black",lty=3,ylab="real prob")
  #   abline(v=which(generatedI[[i]]==1),col="red",lty=3)
  # }
  
  
  
  
  return(Returns)
  
}

cl<-makeCluster(14, outfile="log.txt") 
registerDoSNOW(cl)

kks <- foreach(i = 1:100,.noexport = c("lhcpp","MKFIScpp")) %dopar% {
  myF(i)
} 

stopCluster(cl)



kkR <- kks


kkR <- list()
for(i in 1:100){
  kkR[[i]] <- readRDS(paste0("d_10_sub_300_N_withL/",i,".rds"))
}


vNames <- colnames(kkR[[1]]$EstimatesR)[1:(ncol(kkR[[1]]$EstimatesR)-4)]
nNames <- ncol(kkR[[1]]$EstimatesR)-4
# vNames <- colnames(kkR[[1]]$EstimatesR)[1:(ncol(kkR[[1]]$EstimatesR))]
# nNames <- ncol(kkR[[1]]$EstimatesR)

rmatrix <- matrix(ncol=nNames,nrow=length(kkR))
Times <- c()
loglhs <- c()
MSEThetas <- c()
RelChange <- c()
niter = c()


for(i in 1:length(kkR)){
  rmatrix[i,] <- kkR[[i]]$EstimatesR[2,1:nNames]
  niter = c(niter,ncol(kkR[[i]]$MLEsA))
  RelChange <- c(RelChange,kkR[[i]]$EstimatesR[1,nNames+1])
  Times <- c(Times,kkR[[i]]$EstimatesR[1,nNames+2])
  loglhs <- c(loglhs,kkR[[i]]$EstimatesR[1,nNames+3])
  MSEThetas <- c(MSEThetas,kkR[[i]]$EstimatesR[1,nNames+4])
}
colnames(rmatrix) = vNames



par(mfrow=c(2,round(nNames/2)))
par(mar = c(4, 3, 1, 1))
for(i in 1:ncol(rmatrix)){
  boxplot(rmatrix[,i],xlab=vNames[i],
          ylim=range(kkR[[1]]$EstimatesR[1,i],rmatrix[,i]))
  abline(h=kkR[[1]]$EstimatesR[1,i],col="red")
}





truePAR = c(0.1, 0.03, 0.3, 0, 10, 0.5, 0.5, -3, 4, 0.15, -0.2, -0.8, 0.5, -0.3)
MSE_all = Bias_all = Var_all = c()
for(i in 1:ncol(rmatrix)){
  tempcol = rmatrix[,i]
  temptrue = truePAR[i]
  tempMSE = mean((tempcol-temptrue)^2)
  tempbias = mean(tempcol-temptrue)
  tempvar = var(tempcol) * (99) / 100
  MSE_all = c(MSE_all,tempMSE)
  Bias_all = c(Bias_all,tempbias)
  Var_all = c(Var_all,tempvar)
}
names(MSE_all) = names(Bias_all) = names(Var_all) = vNames
RESULTS = rbind(MSE_all,Bias_all,Var_all)
saveRDS(RESULTS,"d_10_sub_500_N/MSE.rds")










readRDS("d_5_sub_100_N/MSE.rds")
readRDS("d_5_sub_300_N/MSE.rds")
readRDS("d_5_sub_500_N/MSE.rds")







