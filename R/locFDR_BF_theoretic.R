
##### DEPENDENCY: source("common_MCMC_functions.R"); library("mvtnorm"); library("purrr");

source("common_MCMC_functions.R"); library("mvtnorm"); library("purrr");

###########  Adding additional function to compute locFDR in straightforward way for correlated case.  #############

########### choose the shape parameters of the Beta distribution for 'q' prior.
choose_shape_parameters = function(K, X, s.e.){

  ## an informed initialization 
  level = 0.01
  pv = pchisq( (X/s.e.)^2, df=1, lower.tail=FALSE)
  traits <- paste0("Trait", 1:K)

  select <- BH_selection(pv, traits, level)
  nA <- length(select)

  ## specify the shape parameters of 'q' prior
  shape2 <- 1                                               ## shape2 parameter for Beta prior of q
  LB <- 0.1; UB <- 0.5;
  qm <- nA/K
  if(qm < LB) qm <- LB
  if(qm > UB) qm <- UB
  shape1 <- (qm/(1-qm)) * shape2
  shapes <- list(shape1 = shape1, shape2 = shape2)
  shapes
}





###################################################################################################
################################# UNCORRELATED SUMMARY STATISTICS ##################################
####################################################################################################

## writing a quotient for computing locFDR in the uncorrelated case

logexponent = function(x, s, v){
  ## x = betahat, s = standard error, v = prior sd
  s2 = s^2
  v2 = v^2
  sigma2 = 1/((1/s2) + (1/v2))
  sigma = sqrt(sigma2)
  mu = (sigma2/s2) * x
  part1 = -log(sqrt(2*pi)) + log(sigma) - log(s) - log(v)
  part2 = (mu^2/(2*sigma2)) - (x^2/(2*s2))
  log_exponent = part1 + part2
}


## code for computing the locFDR (Bayes factor) analytically for uncorrelated summary statistics


 analytic_locFDR_uncor = function(X, s.e., spikevar = 0.0001, slabvar = 0.8){

   K = length(X)

   shapes = choose_shape_parameters(K, X, s.e.)                    ## get the shape params
   shape1 = shapes$shape1; shape2 = shapes$shape2;

   p1 = shape1/(shape1+shape2) 
   p0 = 1-p1                                                  ## P(z=0)
   logp0 = log(p0)
   logp1 = log(p1)

   tau <- sqrt(spikevar)
   de <- sqrt(spikevar/slabvar)                              ## 1/de = ratio of slab sd and spike sd

   ## compute the log numerator and denominator
   sum0 = K*logp0
   sum1 = 0

   for(j in 1:K){

     b = logexponent(x = X[j], s = s.e.[j], v = tau)
     sum0 = sum0 + b

     c0 = b
     c1 = logexponent(x = X[j], s = s.e.[j], v = (tau/de))
     a0 = exp(logp0+c0)
     a1 = exp(logp1+c1)
     sum1 = sum1+log(a0+a1)
   }

   log_nume = sum0
   log_deno = sum1
   locFDR = exp(log_nume - log_deno)

   prior_prob_null = exp(K*logp0);
   logBF = log(1-locFDR) - log(locFDR) + log(prior_prob_null) - log(1-prior_prob_null)
   log10_BF = log10(exp(logBF))

   if(log10_BF == Inf) log10_BF <- 300

   pleio_measure = list(locFDR = locFDR, log10_BF = log10_BF)
   pleio_measure

 }







################################## Most recent code for computing locFDR theoretically ##########################




locFDR_uncor_analytic = function(X, s.e., spikevar=0.0001, slabvar=0.8){

  K = length(X)    ## number of traits
  shapes = choose_shape_parameters(K, X, s.e.)                    ## get the shape params
  shape1 = shapes$shape1; shape2 = shapes$shape2;

  p = shape1/(shape1+shape2) ;
  logp = log(p); logq = log(1-p);

  ## log null density under all traits null.
  trait_index = 1:K; zeromean = numeric(K); Variance = (s.e.)^2
  var_null = Variance + rep(spikevar, K); sd_null = sqrt(var_null);

  null_log_density = dnorm(X, mean = zeromean, sd = sd_null, log = TRUE)
  null_log_density = sum(null_log_density)
  null_log_density = null_log_density + (K*logq)      ## combine with null prior probability

  ## compute the complete data likelihood 

  log_compL = as.list(1:2)
  var_nonnull = Variance + rep(slabvar, K); sd_nonnull = sqrt(var_nonnull);

  log_compL[[1]] = logp + dnorm(X, mean = zeromean, sd = sd_nonnull, log = TRUE)   ## under non-null
  log_compL[[2]] = logq + dnorm(X, mean = zeromean, sd = sd_null, log = TRUE)      ## under null
  
  total = exp(log_compL[[1]]) + exp(log_compL[[2]]) 
  log_density = log(total)
  logL = sum(log_density)

  ## compute local FDR.
  locFDR = exp(null_log_density - logL)

  ## compute log10BF
  prior_prob_null = exp(K*logq);
  logBF = log(1-locFDR) - log(locFDR) + log(prior_prob_null) - log(1-prior_prob_null)
  log10_BF = log10(exp(logBF))

  if(log10_BF == Inf) log10_BF <- 300

  pleio_measure = list(locFDR = locFDR, log10_BF = log10_BF)
  pleio_measure
 
}














###################################################################################################
################################## CORRELATED SUMMARY STATISTICS ##################################
####################################################################################################


########### Compute the log-likelihood under a causal configuration of traits' association status. 

non_null_density_computation = function(c, K, spikevar, slabvar, zeromean, X, S){
  TAU = rep(spikevar, K); TAU[c] = slabvar;
  Sigma = (diag(TAU)) + S
  log_density = dmvnorm(X, mean=zeromean, sigma=Sigma, log=TRUE)
  log_density
}

########### Main function to compute the local FDR and the Bayes factor

locFDR_cor_analytic_full = function(X, s.e., corln, spikevar=0.0001, slabvar=0.8){


  K = length(X)    ## number of traits
  shapes = choose_shape_parameters(K, X, s.e.)                    ## get the shape params
  shape1 = shapes$shape1; shape2 = shapes$shape2;

  p = shape1/(shape1+shape2) ;
  logp = log(p); logq = log(1-p);

  ## If the covariance matrix is not positive definite, the diagonal elements are incremented to make it PD 
  S <- diag(s.e.) %*% corln %*% diag(s.e.)
  epsilon = 10^(-5)
  increment = rep(epsilon,K)
  while(det(S) <= 0) diag(S) = diag(S)+increment

  ### start computing log-likelihood under different causal configurations

  trait_index = 1:K; zeromean = numeric(K);

  ## log null density under all traits null.
  Sigma = ( spikevar * diag(K) ) + S              ## the covariance matrix
  null_log_density = dmvnorm(X, mean=zeromean, sigma=Sigma, log=TRUE)
  null_log_density = null_log_density + (K*logq)  ## combine with the prior prob of Z=0

  total_density = exp(null_log_density)           ## total data likelihood

  ## computing the log likelihood under the non-null configuration of the causal status of the traits
  for(k in 1:K){
    cn = combn(trait_index, k, simplify = FALSE)
    log_density = map_dbl(cn, non_null_density_computation, K=K, spikevar=spikevar, slabvar=slabvar, zeromean=zeromean, X=X, S=S)
    log_prior_prob = (k*logp) + ((K-k)*logq)    
    log_density = log_density + log_prior_prob    ## combine with the prior prob of the causal configuration.
    total_density = total_density + sum(exp(log_density))
  }

  locFDR = exp(null_log_density - log(total_density))  

  prior_prob_null = exp(K*logq);
  logBF = log(1-locFDR) - log(locFDR) + log(prior_prob_null) - log(1-prior_prob_null)
  log10_BF = log10(exp(logBF))

  if(log10_BF == Inf) log10_BF <- 300
 
  pleio_measure = list(locFDR = locFDR, log10_BF = log10_BF)
  pleio_measure

}







