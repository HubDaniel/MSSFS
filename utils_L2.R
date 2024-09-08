# hyperparameter settings
hyperSetting <- function(lag=3,nS=86,nT=101,Ws="exp.ave",xL = 2,max_rel_change = 1e-3){
  timeL <- nT                  # num of observations
  thetaL <- 1                   # num of biomarkers
  numCov <- (1+thetaL)*thetaL/2 # num of elements to be estimated in the covariance matrices
  numS <- nS                     # num of subjects
  xL <- xL                      # num of covariates, same covariates in both p01 and p11 
  lagL0 <- 0                   # num of lags we used in generate data
  lagL1 <- lag                   # num of lags we used in generate data
  estlagL0 <- 0                # num of lags we used in estimation
  estlagL1 <- lag                # num of lags we used in estimation
  numEM <- 10                    # max num of iterations for EM algorithm
  ref1_greater <- 1
  all_sub_changed <- F         # if all generated subjects have a change
  # this will cause bias in est_alpha0 since we manually removed patients without
  # a change, and boost the prob. of infection (p_01)
  # a better way is to increase alpha0 to a large number so all patients will have a 
  # change
  ref0_in_switch <- 0  
  # 0 if we do not use abs(Z-ref0) in the switch equation
  
  max_rel_change = max_rel_change
  
  use_BFGS = T
  have_miss = T # whether have missing values
  
  Ws = Ws
  if(Ws=="sim.ave"){ # simple average
    ws = rep(1/estlagL1,estlagL1)
  }else if(Ws=="exp.ave"){ # exponential weighted average involving latest "lag" terms
    raw_ws = exp(0.5*(1:estlagL1)) # only a heuristic way
    ws = raw_ws / sum(raw_ws)
  }else if(Ws=="no.ave"){ # not average, only single term at t-lag
    # if(estlagL1 == 1){
    #   raw_ws = 1
    # }else if(estlagL1 == 2){
    #   raw_ws = c(1,1)
    # }else{
    #   raw_ws = c(1,rep(0,estlagL1-2),1)
    # }
    
    # raw_ws = dbeta(seq(0.01,0.99,length.out = estlagL1),0.5,0.5)
    
    raw_ws = dgamma(seq(0.01,10,length.out = estlagL1),shape = 2, scale = 2)
    ws = raw_ws / sum(raw_ws)
    
    
  }else{ # otherwise, provide numerical weights and pass in
    ws = Ws
  }
  
  
  
  HSR <- list(timeL=timeL,thetaL=thetaL,numCov=numCov,numS=numS,
              lagL0=lagL0,lagL1=lagL1,estlagL0=estlagL0,estlagL1=estlagL1,
              numEM=numEM,use_BFGS=use_BFGS,have_miss=have_miss,ref1_greater=ref1_greater,
              all_sub_changed=all_sub_changed,ref0_in_switch=ref0_in_switch,
              xL = xL,ws=matrix(ws,nrow=1),max_rel_change=max_rel_change)
  return(HSR)
}

# data generation and initialize parameters for optim function with constraints (only for initialization)
gdataus <- function(HYPERS=HYPERS,num_miss=0){
  
  
  
  ##### setting parameters in the model
  
  ### hyperparameters
  timeL <- HYPERS$timeL # num of observations
  thetaL <- HYPERS$thetaL  # num of biomarkers : temperature, O2 saturation, blood pressure: temperature negative O2, pressure negative O2, temperature positive pressure
  xL <- HYPERS$xL
  ws = as.numeric(HYPERS$ws)
  ### observation equation - F is identity matrix, v is sigma_e2*I and does not change across time
  # obs. error variance
  V <- 0.3840
  W0 <- 0
  W1 <- 0.6213
  
  ### switch equation
  
  
  ## external variables x and their coefficients alpha and beta
  # gender = sample(c(0,1),1)
  # age = rnorm(1,56,5)
  
  gender = sample(x=c(0,1),size=1,prob = c(0.395, 1-0.395))
  age = rnorm(1,0,1)
  
  
  # vintage = rep(rnorm(1,5,0.3),timeL)
  
  xall <- matrix(c(gender,age),nrow=1)
  
  
  
  alpha0 <- -3.7075           # intercept for p_01 in the switch equation
  alpha1 <- 1.0634              # intercept for p_11 in the switch equation
  zeta1 <- 2.6026 
  
  beta0 <- c(0.1646,-0.1900)
  beta1 <- c(-0.8674,0.8423) 
  ## history of the system and their coefficients
  lagL0 <- HYPERS$lagL0       # time lags for p_01 when generate data
  lagL1 <- HYPERS$lagL1       # time lags for p_11 when generate data
  
  
  ### reference levels and changing rate
  ref0 <- 0
  ref1 <- 1.0980
  r0 <- 0.9392 - 1
  r1 <- 0.6652 - 1
  
  # construct gamma0 and G0 using ref and r
  g0 <- - r0 * ref0
  g1 <- - r1 * ref1
  G0 <- r0+1
  G1 <- r1+1
  
  
  ### generate data
  ## initial theta values
  Y <- c()            # generated observations
  Theta <- c()        # generated latent variables
  Ii <- c()          # generated status: assume the first time no COVID, use ref0
  pii <- c()          # generated transition probabilities Pr(I_t|I_{t-1})
  piim0 <- c()         # generated marginal prob.
  piim1 <- c()         # generated marginal prob.
  sigma_t2 <- 0     # initial variance for state0
  
  theta0 <- rnorm(1,mean = ref0,sd=sqrt(sigma_t2))  # start from a value that is close to reference level
  Theta <- cbind(Theta,theta0)
  thetaPrev <- theta0          # keep track of the previous theta value to generated next theta value
  IiPrev <- 0                  # keep track of the previous Ii value to decide which transition probability should be used
  piimPrev0 <- 1               # keep track of the previous marginal prob., proportion of I=1 when t -> inf
  piimPrev1 <- 0               # keep track of the previous marginal prob., proportion of I=1 when t -> inf
  
  for (i in 1:timeL){
    
    pitemp01 <- 1/( 1 + exp( -alpha0 - sum(beta0 * xall) ) )
    
    if(i==1){
      
      pitemp11 <- 1/( 1 + exp( -alpha1 - sum(beta1 * xall) ) )
      
    }else if( i<(lagL1+1) ){
      
      ws_temp = ws[(length(ws)-i+2) : length(ws)]
      
      ws_temp = ws_temp/sum(ws_temp)
      
      pitemp11 = 1/( 1 + exp( - alpha1
                              - zeta1 * sum(Theta[,2:ncol(Theta)] * ws_temp )
                              - sum(beta1 * xall) ) )
      
    }else{
      
      pitemp11 = 1/( 1 + exp( -alpha1 
                              - zeta1 * sum(Theta[,(ncol(Theta)-lagL1+1):ncol(Theta)] * ws )
                              - sum(beta1 * xall) ))
      
    }
    
    
    
    if(IiPrev==0){ # decide which transition prob. to use
      pitemp <- pitemp01
    }else{
      pitemp <- pitemp11
    }
    
    
    
    Iitemp <- 1*(runif(1)<pitemp) # decide whether now is infected or not according to the transition prob.
    
    
    
    if(Iitemp==1){
      thetatemp <- g1 + G1 %*% thetaPrev + rnorm(n=1, 0, sqrt(W1))
    }else{
      thetatemp <- g0 + G0 %*% thetaPrev + rnorm(n=1, 0, sqrt(W0))
    }
    
    ytemp <- thetatemp + rnorm(1, 0, sqrt(V))   # assume the observational error are independent - V is diagonal
    
    
    # marginal prob.
    piim0temp <- piimPrev0 * (1 - pitemp01) + piimPrev1 * (1 - pitemp11)
    piim1temp <- piimPrev0 * pitemp01 + piimPrev1 * pitemp11
    
    
    # store results
    piim0 <- c(piim0,piim0temp)
    piim1 <- c(piim1,piim1temp)
    pii <- c(pii,pitemp) # this is transition prob. not marginal prob.
    Y <- cbind(Y,ytemp)
    Theta <- cbind(Theta,thetatemp)
    Ii <- c(Ii,Iitemp)            # append current status
    
    # update prev values
    thetaPrev <- thetatemp
    IiPrev <- Iitemp
    piimPrev0 <- piim0temp
    piimPrev1 <- piim1temp
  }
  
  
  Y_obs = Y
  Ii_obs = Ii
  
  if(HYPERS$have_miss){
    missing_t = sort(sample(1:HYPERS$timeL, num_miss))
    Y_obs[,missing_t] = NA
    Ii_obs[missing_t] = NA
  }
  
  
  # par(mfrow=c(3,2),mar=c(5,5,3,3))
  # plot(Y[1,],type="b",ylab="Y")
  # abline(v=which(Ii==1),col="red",lty=3)
  # print(xall)
  # plot(pii,type="b",ylab="transition infection probability")
  # abline(v=which(Ii==1),col="red",lty=3)
  # plot(Theta[1,2:ncol(Theta)],type="b",ylab="Theta")
  # abline(v=which(Ii==1),col="red",lty=2)
  # #cbind(pii,Ii)
  # plot(piim1,type="b",ylab="marginal infection probability")
  # abline(v=which(Ii==1),col="red",lty=2)
  # plot(Y_obs[1,],type="b",ylab="Y",xlab="dropped missing values")
  # abline(v=which(Ii==1),col="red",lty=3)
  # gender[1]
  # Bmi[1]
  
  
  
  gdataR <- list()
  gdataR$generatedTheta <- Theta
  gdataR$generatedY_full <- Y
  gdataR$generatedY <- Y_obs
  gdataR$generated_missingt <- missing_t
  gdataR$generatedX <- xall
  gdataR$generatedProb <- pii
  gdataR$generatedI_full <- Ii
  gdataR$generatedI <- Ii_obs
  
  gdataR$true_par <- list(V_true=V,W0_true=W0,W1_true=W1,
                          ref0_true=ref0,delta_true=ref1-ref0,
                          G0_true=G0, G1_true=G1,
                          alpha0_true=alpha0,alpha1_true=alpha1,
                          beta0_true = beta0,beta1_true = beta1,
                          zeta1_true=zeta1)
  
  
  
  return(gdataR)
}

# initialization used for optim function -- change to ref0,delta,G0,G1
initialus <- function(par,Y,useTruePara=F,HYPERS=hyperSetting()){ # par is the true parameter set, Y is observation set
  
  xL <- HYPERS$xL
  
  # use true parameters to initialize optim function
  if(useTruePara){
    Vinit <- log(par$V_true)   # use log because V is always positive and is exp() in lh function
    W0init <- log(par$W0_true)
    W1init <- log(par$W1_true)
    
    ref0init <- par$ref0_true
    deltainit <- log(par$delta_true)
    G0init <- log(par$G0_true/(1-par$G0_true)) #use log because g0 and g1 is always positive and is exp() in lh function
    G1init <- log(par$G1_true/(1-par$G1_true))
    
    alpha0init <- par$alpha0_true
    alpha1init <- par$alpha1_true
    zeta1init <- par$zeta1_true
    beta0init <- par$beta0_true
    beta1init <- par$beta1_true
    
    
    
  }else{
    Vinit <- log(1) # use log because V is always positive and is exp() in lh function
    W0init <- log(1)
    W1init <- log(1)
    
    deltainit <- log(max(Y,na.rm=T)-min(Y,na.rm=T)) # diff use max-min
    G0init <- log(0.5/(1-0.5))      # assume change rate is -0.5 
    G1init <- log(0.5/(1-0.5))
    
    alpha0init <- 0
    alpha1init <- 0
    zeta1init <- 0
    
    beta0init <- rep(0,xL)
    beta1init <- rep(0,xL)
    
    
  }
  
  
  initR <- c(Vinit,W0init,W1init,deltainit,G0init,G1init,
             alpha0init,alpha1init,beta0init,beta1init,zeta1init)
  
  return(initR)
}

# extract MLE estimates of all parameters after optim function with constraints -- change to ref0,delta,G0,G1
MLEus <- function(est,init=T,HYPERS=hyperSetting()){
  timeL <- HYPERS$timeL                 
  thetaL <- HYPERS$thetaL                   
  numCov <- HYPERS$numCov  
  numS <- HYPERS$numS   
  xL <- HYPERS$xL
  
  ## MLE's
  
  if(est$par[1]>250){
    est_V <- exp(250)
  }else if(est$par[1] < (-250)){
    est_V <- exp(-250)
  }else{
    est_V <- exp(est$par[1])
  }
  if(est$par[2]>250){
    est_W0 <- exp(250)
  }else if(est$par[2] < (-250)){
    est_W0 <- exp(-250)
  }else{
    est_W0 <- exp(est$par[2])
  }
  if(est$par[3]>250){
    est_W1 <- exp(250)
  }else if(est$par[3] < (-250)){
    est_W1 <- exp(-250)
  }else{
    est_W1 <- exp(est$par[3])
  }
  
  
  est_ref0 <- 0       # ref1 is positive for O2
  est_delta <- ifelse(est$par[4]<=250,exp(est$par[4]),exp(250))
  
  if(hyperSetting()$ref1_greater){
    est_ref1 <- est_ref0 + est_delta       # ref0 is greater than ref1 for O2
  }else{
    est_ref1 <- est_ref0 - est_delta       # ref0 is greater than ref1 for O2
  }
  
  est_par_G0 <- ifelse(est$par[5]>15,1/(exp(-15)+1),1/(exp(-est$par[5])+1))
  est_par_G1 <- ifelse(est$par[6]>15,1/(exp(-15)+1),1/(exp(-est$par[6])+1))
  est_g0 <- -est_ref0 * (est_par_G0-1)           # estimate of gamma0
  est_g1 <- -est_ref1 * (est_par_G1-1)           # estimate of gamma1
  
  est_alpha0 <- est$par[7]
  est_alpha1 <- est$par[8]
  
  est_beta0 <- est$par[9:(9+xL-1)]
  est_beta1 <- est$par[(9+xL):(9+2*xL-1)]
  
  
  if(!init){
    est_zeta1 <- est$par[9+2*xL]
  }
  
  ## assemble all estimated matrices use above parameters
  est_G0 <- est_par_G0
  est_G1 <- est_par_G1
  
  ## G0 and G1 cannot be 1 since r0 and r1 cannot be 0
  if(!init){
    MLER <- list(est_V=est_V,est_W0=est_W0,est_W1=est_W1,est_ref0=est_ref0,est_delta=est_delta,
                 est_G0=est_G0,est_G1=est_G1,est_alpha0=est_alpha0,est_alpha1=est_alpha1,
                 est_g0=est_g0,est_g1=est_g1,est_ref1=est_ref1,
                 est_beta0=est_beta0,est_beta1=est_beta1,
                 est_zeta1=est_zeta1)
  }else{
    MLER <- list(est_V=est_V,est_W0=est_W0,est_W1=est_W1,est_ref0=est_ref0,est_delta=est_delta,
                 est_G0=est_G0,est_G1=est_G1,est_alpha0=est_alpha0,est_alpha1=est_alpha1,
                 est_g0=est_g0,est_g1=est_g1,est_ref1=est_ref1,
                 est_beta0=est_beta0,est_beta1=est_beta1)
  }
  
  
  
  
  return(MLER)
}
