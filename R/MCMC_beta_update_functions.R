


## This file contains MCMC functions for updating the beta parameters.



                                                        ## update beta in uncorrelated case

uncorrelated_beta_update <- function(K, s.e., tau, de, Z, X)
{
    s2 = (s.e.)^2
    TAU = rep(tau,K)
    indx = which(Z==1)
    if(length(indx) > 0) TAU[indx] = tau/de

    sigma2.inv <- (1/s2) + (1/TAU^2)
    sigma <- sqrt(1/sigma2.inv)
    mean <- (sigma^2/s2)*X
    beta = rnorm(K,mean,sigma) ;
    return(beta) ;
}



                                                         ## Update the beta parameters in correlated case                                              

correlated_beta_update <- function(K, tau, de, Z, X, Sig1.inv)
{

    TAU = rep(tau,K)
    indx = which(Z==1)
    if(length(indx) > 0) TAU[indx] = tau/de
    
    Sig2.inv = diag(1/(TAU^2))

    combo.inv <- solve(Sig1.inv+Sig2.inv) ;

    mean <- combo.inv %*% Sig1.inv %*% X ;

    Sigma <- combo.inv ;

    beta <- mvrnorm(1,mean,Sigma)

}



