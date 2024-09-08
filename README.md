# MSSFS
The code used for the model Multiprocess State Space Model with Feedback and Switching (MSSFS)

sim.R  
The main function for simulations

rcpp_likelihood.R  
The Rcpp function used to compute the Q functino in the EM algorithm. Optimized using optim() function in R to get estimate of parameters.

rcpp_MKFIS.R  
The Rcpp function used to do Multiprocess Kalman Filter and Multiprocess Fixed Interval Smoothing. 

utils.R  
Other functions needed for simulaltion setup, collecting results, etc.

rcpp_likelihood_L2.R  
Added L2 penalty to rcpp_likelihood.R

rcpp_MKFIS_L2.R  
The same as rcpp_MKFIS.R. For version control use only.

utils_L2.R  
The same as utils.R. For version control use only.
