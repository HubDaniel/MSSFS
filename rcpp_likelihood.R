# This is Rcpp version of the likelihood function
# inputs:
# par: parameter vector
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
#include<iostream>
using namespace Rcpp;
using namespace arma;


double lv(vec& par, List& Y, List& X, List& Z, double& init, List& hypers);


double lv(vec& par, List& Y, List& X, List& Z, double& init, List& hypers){

  double xL,timeL,thetaL,numCov,numS,lagL0,lagL1,minlag,maxlag,ref1_greater,ref0_in_switch;
  double Vtemp,W0temp,W1temp,ref0temp,deltatemp,G0temp,G1temp;
  double alpha0,alpha1,zeta1;
  double PrI__0__0__0,PrI__0__0__1,PrIt_1t_10,PrIt_1t_11;
  double logL,logLj,py;
  double middletemp,x1temp,x2temp;
  double PrItran__t__t_1__00,PrItran__t__t_1__01,PrItran__t__t_1__10,PrItran__t__t_1__11;
  double PrI__t_1t__t_1__00,PrI__t_1t__t_1__01,PrI__t_1t__t_1__10,PrI__t_1t__t_1__11;
  double pytIt_1It__00,pytIt_1It__01,pytIt_1It__10,pytIt_1It__11;
  double pIt_1It__00,pIt_1It__01,pIt_1It__10,pIt_1It__11;
  double PrI__t__t__0,PrI__t__t__1;


  mat ref1,g0,g1,beta0,beta1;
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
  
  xL = hypers["xL"];
  timeL = hypers["timeL"];                  
  thetaL = hypers["thetaL"];                    
  numCov = hypers["numCov"]; 
  numS = hypers["numS"]; 
  lagL0 = hypers["estlagL0"]; 
  lagL1 = hypers["estlagL1"];
  ref1_greater = hypers["ref1_greater"];
  ref0_in_switch = hypers["ref0_in_switch"];
  mat ws = hypers["ws"]; 
  
  // parameters
  Vtemp = exp(par(0));
  W0temp = exp(par(1));
  W1temp = exp(par(2));
  ref0temp = 0;  
  deltatemp = exp(par(3));
  G0temp = 1/(exp(-par(4))+1); // need to makes sure 0 <= G0 < 1 cannot be 1 since r cannot be 0
  G1temp = 1/(exp(-par(5))+1);
  
  // convert parameters into matrices
  mat V(1, 1); 
  mat W0(1, 1); 
  mat W1(1, 1); 
  mat ref0(1, 1); 
  mat delta(1, 1); 
  mat G0(1, 1);
  mat G1(1, 1);
  V.submat( 0, 0, 0, 0 ) = Vtemp;
  W0.submat( 0, 0, 0, 0 ) = W0temp;
  W1.submat( 0, 0, 0, 0 ) = W1temp;
  ref0.submat( 0, 0, 0, 0 ) = ref0temp;
  delta.submat( 0, 0, 0, 0 ) = deltatemp;
  G0.submat( 0, 0, 0, 0 ) = G0temp;
  G1.submat( 0, 0, 0, 0 ) = G1temp;
  
  
  
  //if(hyperSetting()$ref1_greater){
  if(ref1_greater==1){
    ref1 = ref0 + delta;       // ref0 is lower than ref0 for temperature
  }else{
    ref1 = ref0 - delta; 
  }
  
  g0 = -ref0 * (G0-1);           // gamma0
  g1 = -ref1 * (G1-1);           // gamma1
  
  
  alpha0 = par(6);                               // alpha0 in switch equation
  alpha1 = par(7); 
  
  if(xL > 0){
    beta0.set_size(xL, 1); 
    beta1.set_size(xL, 1);
    beta0.submat(0,0,xL-1,0) = par.subvec(8,8+xL-1);
    beta1.submat(0,0,xL-1,0) = par.subvec(8+xL,8+2*xL-1);
  }
  
  if(init==0){
    zeta1 = par(8+2*xL);
  }
  
  //std::cout << "a" << std::endl;
  
  
  // algorithm start ---------------------------------------------------------------------------------------------
  // initial for theta: mean is ref0/ref1 and use diffuse covariance matrix
  mat diagL(thetaL, thetaL, fill::eye);
  
  theta__0__0__0 = ref0; // first two zero indicate t|t and the last zero indicate p
  P__0__0__0 = 0 * diagL; // first two zero indicate t|t and the last zero indicate p

  theta__0__0__1 = ref1; // first two zero indicate t|t and the last zero indicate p
  P__0__0__1 = 0 * diagL; // first two zero indicate t|t and the last zero indicate p
  
  
  // initial for transition probability PrI.t.s.r := Pr(I_{t}=r|psi_{s})
  PrI__0__0__0 = 1-(1e-15);                      // assume all patients are healthy at the beginning 
  
  //std::cout << PrI__0__0__0 << std::endl;
  
  PrI__0__0__1 = 1-PrI__0__0__0;
  
  logL = 0;                       // total log-likelihood for all subjects
  
  //std::cout << "b" << std::endl;
  
  
  
  
  for(int j = 0; j < numS; j++){               // add here for multiple subjects
  
    // recursion of MKF
    theta__t_1__t_1__0 = theta__0__0__0;  // previous step estimate used in current step
    theta__t_1__t_1__1 = theta__0__0__1;
    P__t_1__t_1__0 = P__0__0__0;
    P__t_1__t_1__1 = P__0__0__1;
    PrIt_1t_10 = PrI__0__0__0;
    PrIt_1t_11 = PrI__0__0__1;
    
    logLj = 0;                    // log-likelihood for current subject
    

    mat yi = Y[j];                      // observed data of current subject

    mat xi = X[j];                      // covariates of current subject


    for(int i = 0; i < timeL; i++){
      // one step forward prediction (step 2 of MKF)
      theta__t__t_1__00 = g0 + G0 * theta__t_1__t_1__0;
      theta__t__t_1__01 = g1 + G1 * theta__t_1__t_1__0;
      theta__t__t_1__10 = g0 + G0 * theta__t_1__t_1__1;
      theta__t__t_1__11 = g1 + G1 * theta__t_1__t_1__1;
      P__t__t_1__00 = G0 * P__t_1__t_1__0 * G0.t() + W0;
      P__t__t_1__01 = G1 * P__t_1__t_1__0 * G1.t() + W1;
      P__t__t_1__10 = G0 * P__t_1__t_1__1 * G0.t() + W0;
      P__t__t_1__11 = G1 * P__t_1__t_1__1 * G1.t() + W1;
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
        H__t__t_1__00 = P__t__t_1__00 + V;
        H__t__t_1__01 = P__t__t_1__01 + V;
        H__t__t_1__10 = P__t__t_1__10 + V;
        H__t__t_1__11 = P__t__t_1__11 + V;
        
        
        // posterior of states and its variances (step 4 of MKF)
        theta__t__t__00 = theta__t__t_1__00 + P__t__t_1__00 * inv_sympd(P__t__t_1__00 + V) * eta__t__t_1__00;
        theta__t__t__01 = theta__t__t_1__01 + P__t__t_1__01 * inv_sympd(P__t__t_1__01 + V) * eta__t__t_1__01;
        theta__t__t__10 = theta__t__t_1__10 + P__t__t_1__10 * inv_sympd(P__t__t_1__10 + V) * eta__t__t_1__10;
        theta__t__t__11 = theta__t__t_1__11 + P__t__t_1__11 * inv_sympd(P__t__t_1__11 + V) * eta__t__t_1__11;
        P__t__t__00 = P__t__t_1__00 - P__t__t_1__00 * inv_sympd(P__t__t_1__00 + V) * P__t__t_1__00;
        P__t__t__01 = P__t__t_1__01 - P__t__t_1__01 * inv_sympd(P__t__t_1__01 + V) * P__t__t_1__01;
        P__t__t__10 = P__t__t_1__10 - P__t__t_1__10 * inv_sympd(P__t__t_1__10 + V) * P__t__t_1__10;
        P__t__t__11 = P__t__t_1__11 - P__t__t_1__11 * inv_sympd(P__t__t_1__11 + V) * P__t__t_1__11;
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
        
      // Hamilton Filter (step 5 of MKF) 
      

      mat tempxi = xi;
      
      if(init==0){
        mat Zi = Z[j];
      
        x1temp = -alpha0 - as_scalar(tempxi * beta0);
        if(x1temp > 20 | x1temp < -20){
          if(x1temp > 20){
            x1temp = 20;
          }else{
            x1temp = -20;
          }
        }
        PrItran__t__t_1__01 = 1/( 1 + exp(x1temp) );
        
        
        if(i == 0){

          x2temp = -alpha1 - as_scalar(tempxi * beta1);

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
          x2temp = -alpha1 - zeta1 * as_scalar( Zi.submat( 0, 0, 0, i-1 ) * ws_temp.t() ) - as_scalar(tempxi * beta1) ;

          if(x2temp > 20 | x2temp < -20){
            if(x2temp > 20){
              x2temp = 20;
            }else{
              x2temp = -20;
            }
          }

        }else if (i >= lagL1){
        
          x2temp = -alpha1 - zeta1 *  as_scalar( Zi.submat( 0, i-lagL1, 0, i-1 ) * ws.t() ) - as_scalar(tempxi * beta1) ;

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
        x1temp = -alpha0 - as_scalar(tempxi * beta0);
        if((x1temp > 20) | (x1temp < -20)){
          if(x1temp > 20){
            x1temp = 20;
          }else{
            x1temp = -20;
          }
        }
        x2temp = -alpha1 - as_scalar(tempxi * beta1);
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
        
  
      PrItran__t__t_1__00 = 1 - PrItran__t__t_1__01;
      PrItran__t__t_1__10 = 1 - PrItran__t__t_1__11;
      
      
      // step 5(b) of MKF:   PrI.t_1t.t_1.ab := Pr( I_{t-1}=a,I_{t}=b | psi_{t-1} )
      PrI__t_1t__t_1__00 = PrItran__t__t_1__00 * PrIt_1t_10;
      PrI__t_1t__t_1__01 = PrItran__t__t_1__01 * PrIt_1t_10;
      PrI__t_1t__t_1__10 = PrItran__t__t_1__10 * PrIt_1t_11;
      PrI__t_1t__t_1__11 = PrItran__t__t_1__11 * PrIt_1t_11;
      
      
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
  
      
      // step 6 of MKF - collapse
      // no need to distinguish missing and non-missing since 
      // we have updated them accordingly above
      
      theta__t__t__0 = (pIt_1It__00 * theta__t__t__00 + pIt_1It__10 * theta__t__t__10) / PrI__t__t__0;
      theta__t__t__1 = (pIt_1It__01 * theta__t__t__01 + pIt_1It__11 * theta__t__t__11) / PrI__t__t__1;
      P__t__t__0 = (pIt_1It__00 * (P__t__t__00 + (theta__t__t__00 - theta__t__t__0) * (theta__t__t__00 - theta__t__t__0).t() ) +
                    pIt_1It__10 * (P__t__t__10 + (theta__t__t__10 - theta__t__t__0) * (theta__t__t__10 - theta__t__t__0).t() )) / PrI__t__t__0;
      P__t__t__1 = (pIt_1It__01 * (P__t__t__01 + (theta__t__t__01 - theta__t__t__1) * (theta__t__t__01 - theta__t__t__1).t() ) +
                    pIt_1It__11 * (P__t__t__11 + (theta__t__t__11 - theta__t__t__1) * (theta__t__t__11 - theta__t__t__1).t() )) / PrI__t__t__1;
      
      // update previous for next recursion
      theta__t_1__t_1__0 = theta__t__t__0;
      theta__t_1__t_1__1 = theta__t__t__1;
      P__t_1__t_1__0 = P__t__t__0;
      P__t_1__t_1__1 = P__t__t__1;
      PrIt_1t_10 = PrI__t__t__0;
      PrIt_1t_11 = PrI__t__t__1;
        
    }
    
    logL = logL + logLj;
  }
  
  return logL;
}  
'  



src = '
double init;
vec par;
List Y, X, Z, hypers;

init = as<double>(init_in);
par = as<vec>(par_in);
Y = as<List>(Y_in);
X = as<List>(X_in);
Z = as<List>(Z_in);
hypers = as<List>(hypers_in);

return wrap(lv(par, Y, X, Z, init, hypers));
'


lhcpp = cxxfunction(signature( init_in="numeric",par_in='numeric',Y_in='numeric',
                               X_in='numeric',Z_in='numeric',hypers_in='numeric'), 
                    includes = rcpp_inc, 
                    body = src, 
                    plugin = "RcppArmadillo")






