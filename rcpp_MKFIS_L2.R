# This is Rcpp version of the MKF and FIS process
# inputs:
# MLES: list of MLE estimates
# Y: list of observation matrix 1*timeL (univariate case, o.w. thetaL*timeL), 
#    each list represents a subject, missing values indicated by NA
# Z: smoothed estimates
# init: indicator, whether it is the initialization of EM
# hypers: list of hyperparameters

# outputs:
# likelihood value: Lv
library(Rcpp)
library(inline)


rcpp_inc = '
//include<iostream>
using namespace Rcpp;
using namespace arma;


List MKFIS(List& MLES, mat& Y, mat& X, mat& Z, double& init, List& hypers);


List MKFIS(List& MLES, mat& Y, mat& X, mat& Z, double& init, List& hypers){

  double xL,timeL,thetaL,numCov,numS,lagL0,lagL1,ref0_in_switch;
  double Vtemp,W0temp,W1temp,ref0temp,ref1temp,G0temp,G1temp,g0temp,g1temp;
  double est_alpha0,est_alpha1,est_zeta1;
  double PrI__0__0__0,PrI__0__0__1,PrIt_1t_10,PrIt_1t_11;
  double logLj,py;
  double middletemp,x1temp,x2temp;
  double PrItran__t__t_1__00,PrItran__t__t_1__01,PrItran__t__t_1__10,PrItran__t__t_1__11;
  double PrI__t_1t__t_1__00,PrI__t_1t__t_1__01,PrI__t_1t__t_1__10,PrI__t_1t__t_1__11;
  double pytIt_1It__00,pytIt_1It__01,pytIt_1It__10,pytIt_1It__11;
  double pIt_1It__00,pIt_1It__01,pIt_1It__10,pIt_1It__11;
  double PrI__t__t__0,PrI__t__t__1;
  double pIt1__T__0,pIt1__T__1;
  double pItIt1__00,pItIt1__01,pItIt1__10,pItIt1__11;
  double pIt__0,pIt__1,prI__t1__t__00,prI__t1__t__01,prI__t1__t__10,prI__t1__t__11;

  vec PrIt__T__0__all,PrIt__T__1__all;
  
  mat est_beta0,est_beta1;
  mat theta__0__0__0,theta__0__0__1,theta__t_1__t_1__0,theta__t_1__t_1__1;
  mat P__0__0__0,P__0__0__1,P__t_1__t_1__0,P__t_1__t_1__1;
  mat P__t__t_1__00,P__t__t_1__01,P__t__t_1__10,P__t__t_1__11;
  mat eta__t__t_1__00,eta__t__t_1__01,eta__t__t_1__10,eta__t__t_1__11;
  mat H__t__t_1__00,H__t__t_1__01,H__t__t_1__10,H__t__t_1__11;
  mat theta__t__t__00,theta__t__t__01,theta__t__t__10,theta__t__t__11;
  mat theta__t__t_1__00,theta__t__t_1__01,theta__t__t_1__10,theta__t__t_1__11;
  mat y__t__t_1__00,y__t__t_1__01,y__t__t_1__10,y__t__t_1__11;
  mat P__t__t__00,P__t__t__01,P__t__t__10,P__t__t__11;
  mat theta__t__t__0,theta__t__t__1,P__t__t__0,P__t__t__1;
  mat theta__tBt,P__tBt,theta__t1Bt,Sigma__term1,Sigma__term2;
  mat Sigma__t__t1,theta__t__t,P__t__t,PrI__t__t_1__0,PrI__t__t_1__1;
  mat theta__t1BT,P__t1BT;
  mat P__term1,P__term21,P__term22,P__t1Bt,theta__tBT,P__tBT;
  mat ws = hypers["ws"];
  
  xL = hypers["xL"];
  timeL = hypers["timeL"];                  
  thetaL = hypers["thetaL"];                    
  numCov = hypers["numCov"]; 
  numS = hypers["numS"]; 
  lagL0 = hypers["estlagL0"]; 
  lagL1 = hypers["estlagL1"];
  ref0_in_switch = hypers["ref0_in_switch"];
  
  
  
  // E(theta_t | psi_{t-1}, I_{t-1}, I_t) -- step 2
  // (o,p)=(0,0);(0,1);(1,0);(1,1) in 1st, 2nd, 3rd and 4th list
  List Theta__t__t_1__all(4);
  for(int ii = 0; ii < 4; ii++){
    List tempCreateList(timeL);
    for(int jj = 0; jj < timeL-1; jj++){
      mat tempCreateElement(thetaL,1);
      tempCreateList[jj] = tempCreateElement;
    }
    Theta__t__t_1__all[ii] = tempCreateList;
  }
  
  // Var(theta_t | psi_{t-1}, I_{t-1}, I_t) -- step 2
  // nested list for P_{t|t-1}^{(o,p)}, (o,p)=(0,0);(0,1);(1,0);(1,1) in 1st, 2nd, 3rd and 4th list
  List P__t__t_1__all(4);
  for(int ii = 0; ii < 4; ii++){
    List tempCreateList(timeL);
    for(int jj = 0; jj < timeL-1; jj++){
      mat tempCreateElement(thetaL,thetaL);
      tempCreateList[jj] = tempCreateElement;
    }
    P__t__t_1__all[ii] = tempCreateList;
  }
  
  
  // E(theta_t | psi_t, I_t) -- step 6
  // nested list for theta_{t|t}^{p}, p=0 in first list and p=1 in second list
  List Theta__t__t__all(2);
  for(int ii = 0; ii < 2; ii++){
    List tempCreateList(timeL);
    for(int jj = 0; jj < timeL-1; jj++){
      mat tempCreateElement(thetaL,1);
      tempCreateList[jj] = tempCreateElement;
    }
    Theta__t__t__all[ii] = tempCreateList;
  }  
  
  
  
  // Var(theta_t | psi_t, I_t) -- step 6
  // nested list for P_{t|t}^{p}, p=0 in first list and p=1 in second list
  List P__t__t__all(2);
  for(int ii = 0; ii < 2; ii++){
    List tempCreateList(timeL);
    for(int jj = 0; jj < timeL-1; jj++){
      mat tempCreateElement(thetaL,thetaL);
      tempCreateList[jj] = tempCreateElement;
    }
    P__t__t__all[ii] = tempCreateList;
  }  


  
  // Pr( I_{t} | psi_{t-1} )
  mat PrI__t__t_1(2,timeL);
  // Pr( I_t | psi_{t} )
  mat PrI__t__t(2,timeL);
  // Pr(I_t | I_{t-1}, psi_{t-1})
  mat PrItran__t__t_1__all(4,timeL);
  // (o,p)=(0,0);(0,1);(1,0);(1,1) in 1st, 2nd, 3rd and 4th row
  // Pr( I_{t-1}, I_t | psi_{t-1})
  mat PrI__t_1t__t_1__all(4,timeL);
  // (o,p)=(0,0);(0,1);(1,0);(1,1) in 1st, 2nd, 3rd and 4th row
  // E(theta_t|psi_t)
  mat ThetaForward(thetaL,timeL);
  // Var(theta_t|psi_t)
  cube PForward(thetaL,thetaL,timeL); // first 2 dimensions is the covariance matrix, last dimension is the time points
  
  //std::cout << "a" << std::endl;

  

  
  Vtemp = MLES["est_V"];
  W0temp = MLES["est_W0"];
  W1temp = MLES["est_W1"];
  G0temp = MLES["est_G0"];
  G1temp = MLES["est_G1"];
  g0temp = MLES["est_g0"];
  g1temp = MLES["est_g1"];
  ref0temp = MLES["est_ref0"];
  ref1temp = MLES["est_ref1"];
  
  est_alpha0 = MLES["est_alpha0"];
  est_alpha1 = MLES["est_alpha1"];
  
  if(xL > 0){
    est_beta0.set_size(xL, 1); 
    est_beta1.set_size(xL, 1);
    vec tempbeta0 = MLES["est_beta0"];
    vec tempbeta1 = MLES["est_beta1"];
    est_beta0.submat(0,0,xL-1,0) = tempbeta0;
    est_beta1.submat(0,0,xL-1,0) = tempbeta1;
  }
  
  
  //std::cout << "b" << std::endl;
  
  
  // convert parameters into 1 * 1 matrices
  mat est_V(1, 1); 
  mat est_W0(1, 1); 
  mat est_W1(1, 1); 
  mat est_ref0(1, 1); 
  mat est_ref1(1, 1); 
  mat est_G0(1, 1);
  mat est_G1(1, 1);
  mat est_g0(1, 1);
  mat est_g1(1, 1);
  
  est_V.submat( 0, 0, 0, 0 ) = Vtemp;
  est_W0.submat( 0, 0, 0, 0 ) = W0temp;
  est_W1.submat( 0, 0, 0, 0 ) = W1temp;
  est_ref0.submat( 0, 0, 0, 0 ) = ref0temp;
  est_ref1.submat( 0, 0, 0, 0 ) = ref1temp;
  est_G0.submat( 0, 0, 0, 0 ) = G0temp;
  est_G1.submat( 0, 0, 0, 0 ) = G1temp;
  est_g0.submat( 0, 0, 0, 0 ) = g0temp;
  est_g1.submat( 0, 0, 0, 0 ) = g1temp;
  
  
  
  //std::cout << "c" << std::endl;
  
  // algorithm start ---------------------------------------------------------------------------------------------
  // initial for theta: mean is est_ref0/est_ref1 and use diffuse covariance matrix
  mat diagL(thetaL, thetaL, fill::eye);
  
  theta__0__0__0 = est_ref0; // first two zero indicate t|t and the last zero indicate p
  P__0__0__0 = 0 * diagL; // first two zero indicate t|t and the last zero indicate p

  theta__0__0__1 = est_ref1; // first two zero indicate t|t and the last zero indicate p
  P__0__0__1 = 0 * diagL; // first two zero indicate t|t and the last zero indicate p
  
  
  // initial for transition probability PrI.t.s.r := Pr(I_{t}=r|psi_{s})
  PrI__0__0__0 = 1-(1e-15);                      // assume all patients are healthy at the beginning 
  PrI__0__0__1 = 1-PrI__0__0__0;


  // recursion of MKF
  theta__t_1__t_1__0 = theta__0__0__0;  // previous step estimate used in current step
  theta__t_1__t_1__1 = theta__0__0__1;
  P__t_1__t_1__0 = P__0__0__0;
  P__t_1__t_1__1 = P__0__0__1;
  PrIt_1t_10 = PrI__0__0__0;
  PrIt_1t_11 = PrI__0__0__1;
  
  mat yi = Y;                      // observed data of current subject
  
  
  
  //std::cout << "d" << std::endl;
  

  for(int i = 0; i < timeL; i++){
  
    //std::cout << i << std::endl;
  
    // one step forward prediction (step 2 of MKF)
    theta__t__t_1__00 = est_g0 + est_G0 * theta__t_1__t_1__0;
    theta__t__t_1__01 = est_g1 + est_G1 * theta__t_1__t_1__0;
    theta__t__t_1__10 = est_g0 + est_G0 * theta__t_1__t_1__1;
    theta__t__t_1__11 = est_g1 + est_G1 * theta__t_1__t_1__1;
    P__t__t_1__00 = est_G0 * P__t_1__t_1__0 * est_G0.t() + est_W0;
    P__t__t_1__01 = est_G1 * P__t_1__t_1__0 * est_G1.t() + est_W1;
    P__t__t_1__10 = est_G0 * P__t_1__t_1__1 * est_G0.t() + est_W0;
    P__t__t_1__11 = est_G1 * P__t_1__t_1__1 * est_G1.t() + est_W1;
    y__t__t_1__00 = theta__t__t_1__00;
    y__t__t_1__01 = theta__t__t_1__01;
    y__t__t_1__10 = theta__t__t_1__10;
    y__t__t_1__11 = theta__t__t_1__11;
    
    vec tempyi = yi.submat( 0, i, 0, i );
   
    
    if(!NumericVector::is_na( tempyi[0] )){  // whether the i-th element of yi vector is NA
      // if not missing : regular MKF
      // one step observation prediction error and prediction error variance (step 3 of MKF)
      eta__t__t_1__00 = yi.submat( 0, i, 0, i ) - y__t__t_1__00;
      eta__t__t_1__01 = yi.submat( 0, i, 0, i ) - y__t__t_1__01;
      eta__t__t_1__10 = yi.submat( 0, i, 0, i ) - y__t__t_1__10;
      eta__t__t_1__11 = yi.submat( 0, i, 0, i ) - y__t__t_1__11;
      H__t__t_1__00 = P__t__t_1__00 + est_V;
      H__t__t_1__01 = P__t__t_1__01 + est_V;
      H__t__t_1__10 = P__t__t_1__10 + est_V;
      H__t__t_1__11 = P__t__t_1__11 + est_V;
      
      // posterior of states and its variances (step 4 of MKF)
      theta__t__t__00 = theta__t__t_1__00 + P__t__t_1__00 * inv_sympd(P__t__t_1__00 + est_V) * eta__t__t_1__00;
      theta__t__t__01 = theta__t__t_1__01 + P__t__t_1__01 * inv_sympd(P__t__t_1__01 + est_V) * eta__t__t_1__01;
      theta__t__t__10 = theta__t__t_1__10 + P__t__t_1__10 * inv_sympd(P__t__t_1__10 + est_V) * eta__t__t_1__10;
      theta__t__t__11 = theta__t__t_1__11 + P__t__t_1__11 * inv_sympd(P__t__t_1__11 + est_V) * eta__t__t_1__11;
      P__t__t__00 = P__t__t_1__00 - P__t__t_1__00 * inv_sympd(P__t__t_1__00 + est_V) * P__t__t_1__00;
      P__t__t__01 = P__t__t_1__01 - P__t__t_1__01 * inv_sympd(P__t__t_1__01 + est_V) * P__t__t_1__01;
      P__t__t__10 = P__t__t_1__10 - P__t__t_1__10 * inv_sympd(P__t__t_1__10 + est_V) * P__t__t_1__10;
      P__t__t__11 = P__t__t_1__11 - P__t__t_1__11 * inv_sympd(P__t__t_1__11 + est_V) * P__t__t_1__11;
    }else{
      // step 3 of MKF: skip 
      // step 4 of MKF
      theta__t__t__00 = theta__t__t_1__00;
      theta__t__t__01 = theta__t__t_1__01;
      theta__t__t__10 = theta__t__t_1__10;
      theta__t__t__11 = theta__t__t_1__11;
      P__t__t__00 = P__t__t_1__00;
      P__t__t__01 = P__t__t_1__01;
      P__t__t__10 = P__t__t_1__10;
      P__t__t__11 = P__t__t_1__11;
    }
      
 

    mat xi = X;

    if(init==0){
    
      //std::cout << "g1" << std::endl;
    
    
      mat Zi = Z;
      est_zeta1 = MLES["est_zeta1"];
      
      
      
      
      
      x1temp = -est_alpha0 - as_scalar(xi * est_beta0);
      if(x1temp > 20 | x1temp < -20){
        if(x1temp > 20){
          x1temp = 20;
        }else{
          x1temp = -20;
        }
      }
      PrItran__t__t_1__01 = 1/( 1 + exp(x1temp) );
      
      if(i == 0){

        x2temp = -est_alpha1 - as_scalar(xi * est_beta1);

        if(x2temp > 20 | x2temp < -20){
          if(x2temp > 20){
            x2temp = 20;
          }else{
            x2temp = -20;
          }
        }
        

      }else if(i < lagL1){


        mat ws_temp = ws.submat(0, ws.n_cols-i, 0, ws.n_cols-1);
        ws_temp = ws_temp/as_scalar(sum(ws_temp,1));
        x2temp = -est_alpha1-est_zeta1 * as_scalar( Zi.submat( 0, 0, 0, i-1 ) * ws_temp.t() ) - as_scalar(xi * est_beta1) ;

        if(x2temp > 20 | x2temp < -20){
          if(x2temp > 20){
            x2temp = 20;
          }else{
            x2temp = -20;
          }
        }

      }else if (i >= lagL1){
      
        x2temp = -est_alpha1 - est_zeta1 *  as_scalar( Zi.submat( 0, i-lagL1, 0, i-1 ) * ws.t() ) - as_scalar(xi * est_beta1);

        if(x2temp > 20 | x2temp < -20){
          if(x2temp > 20){
            x2temp = 20;
          }else{
            x2temp = -20;
          }
        }
        
      }
      
      PrItran__t__t_1__11 = 1/( 1 + exp(x2temp) );      
      
  
    
    }else{
    

      x1temp = -est_alpha0 - as_scalar(xi * est_beta0);
      if((x1temp > 20) | (x1temp < -20)){
        if(x1temp > 20){
          x1temp = 20;
        }else{
          x1temp = -20;
        }
      }
      x2temp = -est_alpha1 - as_scalar(xi * est_beta1);
      if((x2temp > 20) | (x2temp < -20)){
        if(x2temp > 20){
          x2temp = 20;
        }else{
          x2temp = -20;
        }
      }
      PrItran__t__t_1__01 = 1/( 1 + exp(x1temp) );
      PrItran__t__t_1__11 = 1/( 1 + exp(x2temp) );
    }

    
    //std::cout << "h" << std::endl;
    

    PrItran__t__t_1__00 = 1 - PrItran__t__t_1__01;
    PrItran__t__t_1__10 = 1 - PrItran__t__t_1__11;
    
    
    // step 5(b) of MKF:   PrI.t_1t.t_1.ab := Pr( I_{t-1}=a,I_{t}=b | psi_{t-1} )
    PrI__t_1t__t_1__00 = PrItran__t__t_1__00 * PrIt_1t_10;
    PrI__t_1t__t_1__01 = PrItran__t__t_1__01 * PrIt_1t_10;
    PrI__t_1t__t_1__10 = PrItran__t__t_1__10 * PrIt_1t_11;
    PrI__t_1t__t_1__11 = PrItran__t__t_1__11 * PrIt_1t_11;
    
    
    // marginal out I_{t-1} which is used for MFIS:   PrI.t.t_1.a := Pr( I_{t}=a | psi_{t-1} )
    PrI__t__t_1__0 = PrI__t_1t__t_1__00 + PrI__t_1t__t_1__10;
    PrI__t__t_1__1 = PrI__t_1t__t_1__01 + PrI__t_1t__t_1__11;
    
    
    
    //std::cout << "i" << std::endl;

    
    // step 5(c)~5(e) of MKF:
    if(!NumericVector::is_na( tempyi[0] )){
      // not missing
      // step 5(c) of MKF
      middletemp = as_scalar(-0.5*eta__t__t_1__00.t()*inv_sympd(H__t__t_1__00)*eta__t__t_1__00);
      if((middletemp > 500) | (middletemp < -500)){
        if(middletemp > 500){
          middletemp = 500;
        }else{
          middletemp = -500;
        }
      }
      pytIt_1It__00 = (1/sqrt(det(H__t__t_1__00))) * exp(middletemp) * PrI__t_1t__t_1__00;
      
      middletemp = as_scalar(-0.5*eta__t__t_1__01.t()*inv_sympd(H__t__t_1__01)*eta__t__t_1__01);
      if((middletemp > 500) | (middletemp < -500)){
        if(middletemp > 500){
          middletemp = 500;
        }else{
          middletemp = -500;
        }
      }
      pytIt_1It__01 = (1/sqrt(det(H__t__t_1__01))) * exp(middletemp) * PrI__t_1t__t_1__01;
      
      middletemp = as_scalar(-0.5*eta__t__t_1__10.t()*inv_sympd(H__t__t_1__10)*eta__t__t_1__10);
      if((middletemp > 500) | (middletemp < -500)){
        if(middletemp > 500){
          middletemp = 500;
        }else{
          middletemp = -500;
        }
      }
      pytIt_1It__10 = (1/sqrt(det(H__t__t_1__10))) * exp(middletemp) * PrI__t_1t__t_1__10;
      
      middletemp = as_scalar(-0.5*eta__t__t_1__11.t()*inv_sympd(H__t__t_1__11)*eta__t__t_1__11);
      if((middletemp > 500) | (middletemp < -500)){
        if(middletemp > 500){
          middletemp = 500;
        }else{
          middletemp = -500;
        }
      }
      pytIt_1It__11 = (1/sqrt(det(H__t__t_1__11))) * exp(middletemp) * PrI__t_1t__t_1__11;
      
      // step 5(d) of MKF
      py = pytIt_1It__00 + pytIt_1It__01 + pytIt_1It__10 + pytIt_1It__11;
      pIt_1It__00 = pytIt_1It__00 / py;
      pIt_1It__01 = pytIt_1It__01 / py;
      pIt_1It__10 = pytIt_1It__10 / py;
      pIt_1It__11 = pytIt_1It__11 / py;
      
      // step 5(e) of MKF - collapse
      PrI__t__t__0 = pIt_1It__00 + pIt_1It__10;
      PrI__t__t__1 = pIt_1It__01 + pIt_1It__11;
      
      // Likelihood of subject i
      logLj = logLj + log(py);
      
    }else{
      // when missing - skip 5(c)
      // 5(d) is the same as 5(b)          
      pIt_1It__00 = PrI__t_1t__t_1__00;
      pIt_1It__01 = PrI__t_1t__t_1__01;
      pIt_1It__10 = PrI__t_1t__t_1__10;
      pIt_1It__11 = PrI__t_1t__t_1__11;
      // 5(e) can be calculated from 5(b) directly
      PrI__t__t__0 = PrI__t_1t__t_1__00 + PrI__t_1t__t_1__10;
      PrI__t__t__1 = PrI__t_1t__t_1__01 + PrI__t_1t__t_1__11;
      
      // do not update likelihood since no observation
    }


    //std::cout << "j" << std::endl;
    
    
    // step 6 of MKF - collapse
    // no need to distinguish missing and non-missing since 
    // we have updated them accordingly above
    
    theta__t__t__0 = (pIt_1It__00 * theta__t__t__00 + pIt_1It__10 * theta__t__t__10) / PrI__t__t__0;
    theta__t__t__1 = (pIt_1It__01 * theta__t__t__01 + pIt_1It__11 * theta__t__t__11) / PrI__t__t__1;
    P__t__t__0 = (pIt_1It__00 * (P__t__t__00 + (theta__t__t__00 - theta__t__t__0) * (theta__t__t__00 - theta__t__t__0).t() ) +
                  pIt_1It__10 * (P__t__t__10 + (theta__t__t__10 - theta__t__t__0) * (theta__t__t__10 - theta__t__t__0).t() )) / PrI__t__t__0;
    P__t__t__1 = (pIt_1It__01 * (P__t__t__01 + (theta__t__t__01 - theta__t__t__1) * (theta__t__t__01 - theta__t__t__1).t() ) +
                  pIt_1It__11 * (P__t__t__11 + (theta__t__t__11 - theta__t__t__1) * (theta__t__t__11 - theta__t__t__1).t() )) / PrI__t__t__1;
    
    
    
    
    //std::cout << "k" << std::endl;
    
    
    
    
    // store filtered estimate of theta E(theta_t|psi_t)
    // step 7 of MKF - collapse again - not needed when do estimation, only needed after estimation
    theta__t__t = PrI__t__t__0 * theta__t__t__0 + PrI__t__t__1 * theta__t__t__1;
    
    //std::cout << "k1" << std::endl;
    
    P__t__t = PrI__t__t__0 * (P__t__t__0 + (theta__t__t__0-theta__t__t) * (theta__t__t__0-theta__t__t).t()) +
      PrI__t__t__1 * (P__t__t__1 + (theta__t__t__1-theta__t__t) * (theta__t__t__1-theta__t__t));
    
    //std::cout << "k2" << std::endl;

    ThetaForward.submat(0,i,ThetaForward.n_rows-1,i) = theta__t__t;
    
    //std::cout << "k3" << std::endl;
    
    PForward.subcube( 0, 0, i, PForward.n_rows-1, PForward.n_cols-1, i ) = P__t__t;
    
    
    //std::cout << "l" << std::endl;
    
    
    // store MKF results for MFIS - use templist for the following nested list operations
    
    //Theta__t__t_1__all[0][i] = theta__t__t_1__00;
    //Theta__t__t_1__all[1][i] = theta__t__t_1__01;
    //Theta__t__t_1__all[2][i] = theta__t__t_1__10;
    //Theta__t__t_1__all[3][i] = theta__t__t_1__11;
    
    List templist;
    
    //std::cout << "l1" << std::endl;
    
    templist = Theta__t__t_1__all[0];
    
    //std::cout << "l2" << std::endl;
    
    templist[i] = theta__t__t_1__00;
    
    //std::cout << "l3" << std::endl;
    
    Theta__t__t_1__all[0] = templist;
    
    //std::cout << "l4" << std::endl;
    
    templist = Theta__t__t_1__all[1];
    templist[i] = theta__t__t_1__01;
    Theta__t__t_1__all[1] = templist;
    
    templist = Theta__t__t_1__all[2];
    templist[i] = theta__t__t_1__10;
    Theta__t__t_1__all[2] = templist;
    
    templist = Theta__t__t_1__all[3];
    templist[i] = theta__t__t_1__11;
    Theta__t__t_1__all[3] = templist;    
    
    //std::cout << "m" << std::endl;
    
    //P__t__t_1__all[0][i] = P__t__t_1__00;
    //P__t__t_1__all[1][i] = P__t__t_1__01;
    //P__t__t_1__all[2][i] = P__t__t_1__10;
    //P__t__t_1__all[3][i] = P__t__t_1__11;
    
    templist = P__t__t_1__all[0];
    templist[i] = P__t__t_1__00;
    P__t__t_1__all[0] = templist;
    
    templist = P__t__t_1__all[1];
    templist[i] = P__t__t_1__01;
    P__t__t_1__all[1] = templist;
    
    templist = P__t__t_1__all[2];
    templist[i] = P__t__t_1__10;
    P__t__t_1__all[2] = templist;
    
    templist = P__t__t_1__all[3];
    templist[i] = P__t__t_1__11;
    P__t__t_1__all[3] = templist;  
    
    //std::cout << "n" << std::endl;
    
    
    //Theta__t__t__all[0][i] = theta__t__t__0;
    //Theta__t__t__all[1][i] = theta__t__t__1;
    
    templist = Theta__t__t__all[0];
    templist[i] = theta__t__t__0;
    Theta__t__t__all[0] = templist; 
    
    templist = Theta__t__t__all[1];
    templist[i] = theta__t__t__1;
    Theta__t__t__all[1] = templist;  
    
    
    //P__t__t__all[0][i] = P__t__t__0;
    //P__t__t__all[1][i] = P__t__t__1;
    
    templist = P__t__t__all[0];
    templist[i] = P__t__t__0;
    P__t__t__all[0] = templist; 
    
    templist = P__t__t__all[1];
    templist[i] = P__t__t__1;
    P__t__t__all[1] = templist;
    
    
    
    //std::cout << "o" << std::endl;
    
    
    
    
    PrI__t__t_1.submat(0,i,0,i) = PrI__t__t_1__0;
    PrI__t__t_1.submat(1,i,1,i) = PrI__t__t_1__1;
    
    PrI__t__t.submat(0,i,0,i) = PrI__t__t__0;
    PrI__t__t.submat(1,i,1,i) = PrI__t__t__1;
    
    PrItran__t__t_1__all.submat(0,i,0,i) = PrItran__t__t_1__00;
    PrItran__t__t_1__all.submat(1,i,1,i) = PrItran__t__t_1__01;
    PrItran__t__t_1__all.submat(2,i,2,i) = PrItran__t__t_1__10;
    PrItran__t__t_1__all.submat(3,i,3,i) = PrItran__t__t_1__11;
    
    
    // step 5(b) of MKF:   PrI__t_1t__t_1__ab := Pr( I_{t-1}=a,I_{t}=b | psi_{t-1} )
    PrI__t_1t__t_1__all.submat(0,i,0,i) = PrI__t_1t__t_1__00;
    PrI__t_1t__t_1__all.submat(1,i,1,i) = PrI__t_1t__t_1__01;
    PrI__t_1t__t_1__all.submat(2,i,2,i) = PrI__t_1t__t_1__10;
    PrI__t_1t__t_1__all.submat(3,i,3,i) = PrI__t_1t__t_1__11;

    // update previous for next recursion
    theta__t_1__t_1__0 = theta__t__t__0;
    theta__t_1__t_1__1 = theta__t__t__1;
    P__t_1__t_1__0 = P__t__t__0;
    P__t_1__t_1__1 = P__t__t__1;
    PrIt_1t_10 = PrI__t__t__0;
    PrIt_1t_11 = PrI__t__t__1;
      
  }
  
  
  
  
  
  
  
  
  
  
  ////// MFIS ------------------------------------------------------------------------------------------------------------------------------
  // stored MKF results for MFIS
  //// Theta__t__t_1__all stores from t=1 to T:   theta_{1|0}^{.,.}        ~ theta_{T|T-1}^{.,.}
  // P.t.t_1.all     stores from t=1 to T:   P_{1|0}^{.,.}            ~ P_{T|T-1}^{.,.}
  //// Theta__t__t__all   stores from t=1 to T:   theta_{1|1}^{.}          ~ theta_{T|T}^{.}
  //// P.t.t.all       stores from t=1 to T:   P_{1|1}^{.}              ~ P_{T|T}^{.}
  // PrI__t__t_1       stores from t=1 to T:   Pr(I_1|psi_{0})          ~ Pr(I_T|psi_{T-1})
  // PrI__t__t         stores from t=1 to T:   Pr(I_1|psi_{1})          ~ Pr(I_T|psi_{T})
  // PrItran__t__t_1__all stor from t=1 to T:   Pr(I_1 | I_{0}, psi_{0}) ~ Pr(I_T | I_{T-1}, psi_{T-1})
  //// PrI__t_1t__t_1__all store from t=1 to T:   Pr(I_0, I_1 | psi_{0})   ~ Pr(I_{T-1}, I_{T} | psi_{T-1})
  //// ThetaForward    stores from t=1 to T:   theta{1|1}               ~ theta{T|T}
  //// PForward        stores from t=1 to T:   P{1|1}                   ~ P{T|T}
  
  //// start with t=T-1
  pIt1__T__0 = as_scalar(PrI__t__t.submat(0,timeL-1,0,timeL-1));
  pIt1__T__1 = as_scalar(PrI__t__t.submat(1,timeL-1,1,timeL-1));
  
  theta__t1BT = ThetaForward.submat(0,timeL-1,ThetaForward.n_rows-1,timeL-1); // theta_{T|T}, t1 is t+1, B is |, t1Bt1 means last t|t
  P__t1BT = PForward.subcube( 0, 0, timeL-1, PForward.n_rows-1, PForward.n_cols-1, timeL-1 );   // P_{T|T}
  
  // storing unconditional estimate of theta and probabilities of interest
  mat ThetaSmooth(thetaL,timeL);
  
  PrIt__T__0__all = zeros( timeL );
  PrIt__T__1__all = zeros( timeL );
  
  PrIt__T__0__all(0) = as_scalar(PrI__t__t.submat(0,timeL-1,0,timeL-1)); // Pr(I_t|psi_T) for t=T,...,1
  PrIt__T__1__all(0) = as_scalar(PrI__t__t.submat(1,timeL-1,1,timeL-1));
  mat prI__t1__t__all(4,timeL); // Pr(I_t|I_{t-1},psi_T) for t=T,...,2, we never have t=1 for these transition prob.
  
  
  mat ThetaPred(1,timeL);
  
  for(int i = timeL-2; i >= 0; i--){
    
    theta__tBt = ThetaForward.submat(0,i,ThetaForward.n_rows-1,i);
    
    P__tBt = PForward.subcube( 0, 0, i, PForward.n_rows-1, PForward.n_cols-1, i );
    
    //theta__t1Bt = Theta__t__t_1__all[0][i] * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
    //  Theta__t__t_1__all[1][i] * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
    //  Theta__t__t_1__all[2][i] * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
    //  Theta__t__t_1__all[3][i] * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1);
      
    
    List templist0, templist1, templist2, templist3;
    mat templist0i, templist1i, templist2i, templist3i;
    templist0 = Theta__t__t_1__all[0];
    templist1 = Theta__t__t_1__all[1];
    templist2 = Theta__t__t_1__all[2];
    templist3 = Theta__t__t_1__all[3];
    templist0i = Rcpp::as<double>(templist0[i+1]);
    templist1i = Rcpp::as<double>(templist1[i+1]);
    templist2i = Rcpp::as<double>(templist2[i+1]);
    templist3i = Rcpp::as<double>(templist3[i+1]);    
    
    theta__t1Bt = templist0i * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
      templist1i * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
      templist2i * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
      templist3i * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1);
    
    ThetaPred.submat(0,i+1,0,i+1) = theta__t1Bt;
    
    //Sigma__term1 = 
    //  (Theta__t__t__all[0][i] * est_g0.t() + 
    //     Theta__t__t__all[0][i] * Theta__t__t__all[0][i].t() * est_G0.t() + 
    //     P__t__t__all[0][i] * est_G0.t()) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
    //  (Theta__t__t__all[0][i] * est_g1.t() + 
    //     Theta__t__t__all[0][i] * Theta__t__t__all[0][i].t() * est_G1.t() + 
    //     P__t__t__all[0][i] * est_G1.t()) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) +
    //  (Theta__t__t__all[1][i] * est_g0.t() + 
    //     Theta__t__t__all[1][i] * Theta__t__t__all[1][i].t() * est_G0.t() + 
    //     P__t__t__all[1][i] * est_G0.t()) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) +
    //  (Theta__t__t__all[1][i] * est_g1.t() + 
    //     Theta__t__t__all[1][i] * Theta__t__t__all[1][i].t() * est_G1.t() + 
    //     P__t__t__all[1][i] * est_G1.t()) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1);
         
    templist0 = Theta__t__t__all[0];
    templist1 = Theta__t__t__all[1];
    templist2 = P__t__t__all[0];
    templist3 = P__t__t__all[1];
    templist0i = Rcpp::as<double>(templist0[i]);
    templist1i = Rcpp::as<double>(templist1[i]);
    templist2i = Rcpp::as<double>(templist2[i]);
    templist3i = Rcpp::as<double>(templist3[i]); 
         
    Sigma__term1 = 
      (templist0i * est_g0.t() + 
         templist0i * templist0i.t() * est_G0.t() + 
         templist2i * est_G0.t()) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
      (templist0i * est_g1.t() + 
         templist0i * templist0i.t() * est_G1.t() + 
         templist2i * est_G1.t()) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) +
      (templist1i * est_g0.t() + 
         templist1i * templist1i.t() * est_G0.t() + 
         templist3i * est_G0.t()) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) +
      (templist1i * est_g1.t() + 
         templist1i * templist1i.t() * est_G1.t() + 
         templist3i * est_G1.t()) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1);         
         
 
    Sigma__term2 = theta__tBt * theta__t1Bt.t();
    Sigma__t__t1 = Sigma__term1 - Sigma__term2;
    
    
    
    //P__term1 = 
    //  (est_G0 * P__t__t__all[0][i] * est_G0.t() + est_W0) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
    //  (est_G1 * P__t__t__all[0][i] * est_G1.t() + est_W1) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
    //  (est_G0 * P__t__t__all[1][i] * est_G0.t() + est_W0) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
    //  (est_G1 * P__t__t__all[1][i] * est_G1.t() + est_W1) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1); 
    
    P__term1 = 
      (est_G0 * templist2i * est_G0.t() + est_W0) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
      (est_G1 * templist2i * est_G1.t() + est_W1) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
      (est_G0 * templist3i * est_G0.t() + est_W0) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
      (est_G1 * templist3i * est_G1.t() + est_W1) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1);     
    
 
    //P__term21 = 
    //  (est_g0 * est_g0.t() + est_g0 * Theta__t__t__all[0][i].t() * est_G0.t() + 
    //  (est_g0 * Theta__t__t__all[0][i].t() * est_G0.t()).t() + 
    //   est_G0 * Theta__t__t__all[0][i] * (est_G0 * Theta__t__t__all[0][i]).t()) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
    //  (est_g1 * est_g1.t() + est_g1 * Theta__t__t__all[0][i].t() * est_G1.t() + 
    //  (est_g1 * Theta__t__t__all[0][i].t() * est_G1.t()).t() + 
    //   est_G1 * Theta__t__t__all[0][i] * (est_G1 * Theta__t__t__all[0][i]).t()) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
    //  (est_g0 * est_g0.t() + est_g0 * Theta__t__t__all[1][i].t() * est_G0.t() + 
    //  (est_g0 * Theta__t__t__all[1][i].t() * est_G0.t()).t() + 
    //   est_G0 * Theta__t__t__all[1][i] * (est_G0 * Theta__t__t__all[1][i]).t()) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
    //  (est_g1 * est_g1.t() + est_g1 * Theta__t__t__all[1][i].t() * est_G1.t() + 
    //  (est_g1 * Theta__t__t__all[1][i].t() * est_G1.t()).t() + 
    //   est_G1 * Theta__t__t__all[1][i] * (est_G1 * Theta__t__t__all[1][i]).t()) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1); 
    
    
    P__term21 = 
      (est_g0 * est_g0.t() + est_g0 * templist0i.t() * est_G0.t() + 
      (est_g0 * templist0i.t() * est_G0.t()).t() + 
       est_G0 * templist0i * (est_G0 * templist0i).t()) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) + 
      (est_g1 * est_g1.t() + est_g1 * templist0i.t() * est_G1.t() + 
      (est_g1 * templist0i.t() * est_G1.t()).t() + 
       est_G1 * templist0i * (est_G1 * templist0i).t()) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
      (est_g0 * est_g0.t() + est_g0 * templist1i.t() * est_G0.t() + 
      (est_g0 * templist1i.t() * est_G0.t()).t() + 
       est_G0 * templist1i * (est_G0 * templist1i).t()) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
      (est_g1 * est_g1.t() + est_g1 * templist1i.t() * est_G1.t() + 
      (est_g1 * templist1i.t() * est_G1.t()).t() + 
       est_G1 * templist1i * (est_G1 * templist1i).t()) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1); 
        
    //P__term22 = 
    //  (est_g0 + est_G0 * Theta__t__t__all[0][i]) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) +
    //  (est_g1 + est_G1 * Theta__t__t__all[0][i]) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
    //  (est_g0 + est_G0 * Theta__t__t__all[1][i]) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
    //  (est_g1 + est_G1 * Theta__t__t__all[1][i]) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1); 
      
    P__term22 = 
      (est_g0 + est_G0 * templist0i) * PrI__t_1t__t_1__all.submat(0,i+1,0,i+1) +
      (est_g1 + est_G1 * templist0i) * PrI__t_1t__t_1__all.submat(1,i+1,1,i+1) + 
      (est_g0 + est_G0 * templist1i) * PrI__t_1t__t_1__all.submat(2,i+1,2,i+1) + 
      (est_g1 + est_G1 * templist1i) * PrI__t_1t__t_1__all.submat(3,i+1,3,i+1);      
      
      
    
    P__t1Bt = P__term1 + P__term21 - P__term22 * P__term22.t();
    
    
    
    theta__tBT = theta__tBt + Sigma__t__t1 * inv_sympd(P__t1Bt) * (theta__t1BT - theta__t1Bt);
    P__tBT = P__tBt - Sigma__t__t1 * inv_sympd(P__t1Bt) * Sigma__t__t1.t() + 
             Sigma__t__t1 * inv_sympd(P__t1Bt) * P__t1BT * inv_sympd(P__t1Bt) * Sigma__t__t1.t();
    
    
    
    
    // step 5(a) of MFIS
    pItIt1__00 = as_scalar(pIt1__T__0 * PrI__t__t.submat(0,i,0,i) * PrItran__t__t_1__all.submat(0,i+1,0,i+1) / PrI__t__t_1.submat(0,i+1,0,i+1));
    pItIt1__01 = as_scalar(pIt1__T__1 * PrI__t__t.submat(0,i,0,i) * PrItran__t__t_1__all.submat(1,i+1,1,i+1) / PrI__t__t_1.submat(1,i+1,1,i+1));
    pItIt1__10 = as_scalar(pIt1__T__0 * PrI__t__t.submat(1,i,1,i) * PrItran__t__t_1__all.submat(2,i+1,2,i+1) / PrI__t__t_1.submat(0,i+1,0,i+1));
    pItIt1__11 = as_scalar(pIt1__T__1 * PrI__t__t.submat(1,i,1,i) * PrItran__t__t_1__all.submat(3,i+1,3,i+1) / PrI__t__t_1.submat(1,i+1,1,i+1));
    
    
    // step 5(b) of MFIS
    pIt__0 = pItIt1__00 + pItIt1__01;
    pIt__1 = pItIt1__10 + pItIt1__11;
    
    // this is posterior transition probability, return and compare with posterior marginal probability
    prI__t1__t__00 = pItIt1__00 / pIt__0;
    prI__t1__t__01 = pItIt1__01 / pIt__0;
    prI__t1__t__10 = pItIt1__10 / pIt__1;
    prI__t1__t__11 = pItIt1__11 / pIt__1;
    
    
    
    
    // store unconditional smoothing estimate for next EM iteration
    ThetaSmooth.submat(0,i,ThetaSmooth.n_rows-1,i) = theta__tBT;
    
    PrIt__T__0__all(timeL-i-1) = pIt__0;         // Pr(I_t|psi_T) for t=T,...,1
    PrIt__T__1__all(timeL-i-1) = pIt__1;
    prI__t1__t__all.submat(0,i+1,0,i+1) = prI__t1__t__00;    // Pr(I_t|I_{t-1},psi_T) for t=T,...,1
    prI__t1__t__all.submat(1,i+1,1,i+1) = prI__t1__t__01;
    prI__t1__t__all.submat(2,i+1,2,i+1) = prI__t1__t__10;
    prI__t1__t__all.submat(3,i+1,3,i+1) = prI__t1__t__11;
    
    // update for next recursion
    theta__t1BT = theta__tBT;
    P__t1BT = P__tBT;
    pIt1__T__0 = pIt__0;
    pIt1__T__1 = pIt__1;
    
  }
  
  // the first smoothing estimate is the collapsed last MKF estimate, same formula as in MFIS step 8
  ThetaSmooth.submat(0,timeL-1,ThetaSmooth.n_rows-1,timeL-1) = ThetaForward.submat(0,timeL-1,ThetaSmooth.n_rows-1,timeL-1);
  
  
  List templist0, templist1, templist2, templist3;
  mat templist0i, templist1i, templist2i, templist3i;
  templist0 = Theta__t__t_1__all[0];
  templist1 = Theta__t__t_1__all[1];
  templist2 = Theta__t__t_1__all[2];
  templist3 = Theta__t__t_1__all[3];
  templist0i = Rcpp::as<double>(templist0[0]);
  templist1i = Rcpp::as<double>(templist1[0]);
  templist2i = Rcpp::as<double>(templist2[0]);
  templist3i = Rcpp::as<double>(templist3[0]);    
  
  theta__t1Bt = templist0i * PrI__t_1t__t_1__all.submat(0,0,0,0) + 
    templist1i * PrI__t_1t__t_1__all.submat(1,0,1,0) + 
    templist2i * PrI__t_1t__t_1__all.submat(2,0,2,0) + 
    templist3i * PrI__t_1t__t_1__all.submat(3,0,3,0);
  
  ThetaPred.submat(0,0,0,0) = theta__t1Bt;
  
  
  
  
  

  
  return List::create(Named("ThetaSmooth") = ThetaSmooth, Named("ThetaForward") = ThetaForward,
                      Named("PrIt__T__1__all") = reverse(PrIt__T__1__all), Named("prI__t1__t__all") = prI__t1__t__all,
                      Named("PrI__t__t") = PrI__t__t, Named("PrI__t__t_1") = PrI__t__t_1,
                      Named("ThetaPred") = ThetaPred);
                      
  
  
  
  
}  
'  



src = '
double init;
List MLES, hypers;
mat Y, X, Z;

init = as<double>(init_in);
MLES = as<List>(MLES_in);
hypers = as<List>(hypers_in);
Y = as<mat>(Y_in);
X = as<mat>(X_in);
Z = as<mat>(Z_in);

return wrap(MKFIS(MLES, Y, X, Z, init, hypers));
'

MKFIScpp = cxxfunction(signature( 
  init_in="numeric",MLES_in='numeric',hypers_in='numeric',
  Y_in='numeric',X_in='numeric',Z_in='numeric'), 
  includes = rcpp_inc,
  body = src, 
  plugin = "RcppArmadillo")






