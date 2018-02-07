
 ##=================================**********************************************===============================##
 ## This file contains some functions that are called by both of the uncorrelated and correlated versions of CPBayes
 ## For example, to initiate the parameter values in the MCMC. Update Z, q, de in the MCMC. Note that, the full conditional  
 ## posterior distributions of these parameters in the MCMC for both uncorrelated and correlated versions are the same.
 ##=================================**********************************************===============================## 



################################ Function implementing the FDR procedure #############################

BH_selection <- function(pv, traits, level){

  K <- length(pv)                                  ## Number of phenotypes
  pval <- matrix(0, K, 2)
  pval[ ,1] <- seq(1,K)
  pval[ ,2] <- pv
  pval <- pval[order(pval[,2]), ]

  cutoff <- seq(1,K) * (1/K) * level
  comparison <- pval[ ,2] - cutoff

  selected_traits <- NULL

  if(any(comparison <= 0) == TRUE){
    BH_posi <- max(which(comparison <= 0))
    select_index <- sort(pval[1:BH_posi,1])
    selected_traits <- traits[select_index]
  }

  return(selected_traits)
}
##################################***************************************###################################



                         ##*************** initialization for the CpgBayes MCMC ****************##


initiate_MCMC <- function( K, X, s.e. )
{
                                                                   
                                                                   ## minimum and maximum choices of the intial value of q 
     min.q = 0.05
     max.q = 0.95
     
     q.mode = min.q
                                                               
                                                               ## level of FDR correction in the BY or BH procedure
     level = 0.01 
     pv = pchisq( (X/s.e.)^2, df=1, lower.tail=F )
     traits <- paste0("Trait", 1:K)
     
     select <- BH_selection(pv, traits, level)
     K1.FDR <- length(select)
     
     nonnull.set = 0                                           
     Z = rep(0,K)                                              ## Introduce the allocations 
     
     if( K1.FDR > 0 )
     {   
         nonnull.set <- match(select, traits)
         Z[nonnull.set] <- rep(1, K1.FDR)
         q.mode = K1.FDR/K
     }
     
     if( q.mode == 1 ) q.mode <- max.q
     
     beta = X                                                  ## beta = X, under continuous spike, beta is always non-zero
     
     data <- list( beta = beta, Z = Z, q = q.mode, K1_FDR = K1.FDR )

}
                                                        

                         ##*************** Update allocations (Z) for both uncorrelated and correlated cases ****************##
                         ##*************** for the model in which q is included in the MCMC to be updated    ****************##

Z_update <- function(K, q, tau, de, beta)
{
    Z = rep(0,K) ;

    log.ratio = log(q) - log(1-q) + log(de) - ( ((beta^2)/(2*(tau^2))) * ((de^2)-1) ) ;
    pr0 = 1/(1+exp(log.ratio)) ;
    u = runif(K,0,1) ;
    non_zero = which(u>pr0) ;
    if(length(non_zero)>0) Z[non_zero] = 1

    data <- list(prob = pr0, Z = Z)
    return(data)
}


                         ##*************** Update q for both of uncorrelated and correlated versions of CPBayes ****************##

q_update <- function(K, Z, shape1, shape2)
{
    k1 = sum(Z)
    k2 = K-k1
    sh1 = shape1+k1
    sh2 = shape2+k2
    q = rbeta(1,sh1,sh2)
    return(q)
}


                         ##*************** Update allocations (Z) for both uncorrelated and correlated cases ****************##
                         ##*************** for the model in which q is integrated out from the model         ****************##

Z_integrated_update <- function(K, log_ratio_p1, tau_const, de, beta)
{
    Z = rep(0,K) ;

    log.ratio = log_ratio_p1 + log(de) + ( ((beta^2)/tau_const) * (1-(de^2)) )
    pr0 = 1/(1+exp(log.ratio))
    u = runif(K,0,1)
    non_zero = which(u > pr0)
    if(length(non_zero) > 0) Z[non_zero] = 1

    data <- list(prob = pr0, Z = Z)
    return(data)
}


                         ##*************** Function to update the 'de' parameter ****************##
                         ##*************** for both uncorrelated and correlated versions for CPBayes ***************##


                         ## Updating the 'de' parameter when length(Z==1) > 0


de_update1 <- function(min_de, max_de, shape1_de, beta, Z, tau)
{
    K1 = sum(Z)
    non_zero = which(Z>0)
    beta2 = beta^2
    const = (1/(2*(tau^2))) * sum(beta2[non_zero])

    y1 = (min_de^2)*2*const
    y2 = (max_de^2)*2*const

    u = runif(1,0,1)
    df = shape1_de+K1

    prob1 = pchisq(y1, df=df, ncp=0, lower.tail=TRUE)
    prob2 = pchisq(y2, df=df, ncp=0, lower.tail=TRUE)

    prob = prob2-prob1
    p = prob1 + (u*prob)
    y = qchisq(p, df=df, ncp=0, lower.tail=TRUE)

    d = sqrt(y/(2*const))

    return(d)
}


                         ## Updating the 'de' parameter when length(Z==1) = 0

de_update0 <- function(min_de, max_de, shape1_de)
{
    u = runif(1,0,1)

    prob1 = pbeta(min_de, shape1_de, 1, lower.tail = TRUE)
    prob2 = pbeta(max_de, shape1_de, 1, lower.tail = TRUE)

    prob = prob2-prob1
    p = prob1 + (u*prob)
    d = qbeta(p, shape1_de, 1, lower.tail = TRUE)

    return(d)
}






