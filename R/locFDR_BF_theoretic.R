
## DEPENDENCY: source("common_MCMC_functions.R"); library("mvtnorm"); library("purrr");

## source("common_MCMC_functions.R"); library("mvtnorm"); library("purrr");

## Adding additional function to compute locFDR in straightforward way for correlated case.

##============= choose the shape parameters of the Beta distribution for 'q' prior. ==================
choose_shape_parameters = function(K, X, SE){

  ## an informed initialization
  level = 0.01
  pv = pchisq( (X/SE)^2, df=1, lower.tail=FALSE)
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


################################## Most recent code for computing locFDR theoretically ##########################
##-------------------------------------------------------------------------------------------------------------##
#' Analytic calculation of the local FDR & Bayes factor for uncorrelated summary statistics.
#'
#' Run the \code{\link{analytic_locFDR_BF_uncor}} function to analytically compute the local FDR & Bayes factor (BF)
#' that quantifies the evidence of aggregate-level pleiotropic association for uncorrelated summary statistics.
#' Here a fixed value of slab variance is considred instead of a range of it in \code{\link{cpbayes_uncor}}.
#' @param BetaHat A numeric vector of length K where K is the number of phenotypes. It
#'  contains the beta-hat values across studies/traits. No default.
#' @param SE A numeric vector with the same dimension as BetaHat providing the standard errors
#'  corresponding to BetaHat. Every element of SE must be positive. No default.
#' @param SpikeVar Variance of spike (normal distribution with small variance) representing the null effect distribution.
#' Default is 10^(-4).
#' @param SlabVar Variance of slab normal distribution representing the non-null effect distribution.
#' Default is 0.8.
#' @return The output produced by the function is a list which consists of the local FDR and log10(Bayes factor).
#'    \item{locFDR}{It provides the analytically computed local false discovery rate (posterior probability of null association) under CPBayes model
#'     (a Bayesian analog of the p-value) which is a measure of the evidence of the
#'     aggregate-level pleiotropic association. Bayes factor is adjusted for prior odds, but
#'      locFDR is solely a function of the posterior odds.}
#'    \item{log10_BF}{It provides the analytically computed log10(Bayes factor) produced by CPBayes that measures the
#'     evidence of the overall pleiotropic association.}
#'
#' @references Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018) An efficient Bayesian meta analysis approach for studying cross-phenotype genetic associations. PLoS Genet 14(2): e1007139.
#'
#' @seealso \code{\link{cpbayes_uncor}}, \code{\link{analytic_locFDR_BF_cor}}, \code{\link{cpbayes_cor}}, \code{\link{estimate_corln}}, \code{\link{post_summaries}}, \code{\link{forest_cpbayes}}
#'
#' @examples
#' data(ExampleDataUncor)
#' BetaHat <- ExampleDataUncor$BetaHat
#' BetaHat
#' SE <- ExampleDataUncor$SE
#' SE
#' result <- analytic_locFDR_BF_uncor(BetaHat, SE)
#' str(result)
#'
#' @export
analytic_locFDR_BF_uncor = function(BetaHat, SE, SpikeVar=0.0001, SlabVar=0.8){

  # Check whether any of the primary arguments is missing
  if(missing(BetaHat) || missing (SE))
    stop("BetaHat or SE vector is missing!", call. = FALSE)
  # Argument 1 :: BetaHat
  BetaHat <- checkPrimaryVar(BetaHat, "BetaHat")
  # Argument 2 :: SE
  SE <- checkPrimaryVar(SE, "SE")
  # Check whether all entries are strictly positive
  if(!all(SE > 0))
    stop("One or more elements in the SE vector are not positive!", call. = FALSE)
  # Argument 1 and 2 ::
  if(length(BetaHat) != length(SE))
    stop("BetaHat and SE vectors must have the same number of elements!", call. = FALSE)

  # Argument 3 and 4 ::
  chkVar = chkSlabVarSpikeVar(SpikeVar, SlabVar)

  spikevar = chkVar[["SpikeVar"]]
  slabvar = chkVar[["SlabVar"]]
  X = BetaHat
  K = length(X)    ## number of traits
  shapes = choose_shape_parameters(K, X, SE)                    ## get the shape params
  shape1 = shapes$shape1; shape2 = shapes$shape2;

  p = shape1/(shape1+shape2) ;
  logp = log(p); logq = log(1-p);

  ## log null density under all traits null.
  trait_index = 1:K; zeromean = numeric(K); Variance = (SE)^2
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
########### Main function to compute the local FDR and the Bayes factor for correlated summary statistics ##########
#' Analytic calculation of the local FDR & Bayes factor for correlated summary statistics.
#'
#' Run the \code{\link{analytic_locFDR_BF_cor}} function to analytically compute the local FDR & Bayes factor (BF)
#' that quantifies the evidence of aggregate-level pleiotropic association for correlated summary statistics.
#' Here a fixed value of slab variance is considred instead of a range of it in \code{\link{cpbayes_cor}}.
#' @param BetaHat A numeric vector of length K where K is the number of phenotypes. It
#'  contains the beta-hat values across studies/traits. No default.
#' @param SE A numeric vector with the same dimension as BetaHat providing the standard errors
#'  corresponding to BetaHat. Every element of SE must be positive. No default.
#' @param Corln A numeric square matrix of order K by K providing the correlation matrix of BetaHat.
#'  The number of rows & columns of Corln must be the same as the length of BetaHat. No default
#'   is specified. See \code{\link{estimate_corln}}.
#' @param SpikeVar Variance of spike (normal distribution with small variance) representing the null effect distribution.
#' Default is 10^(-4).
#' @param SlabVar Variance of slab normal distribution representing the non-null effect distribution.
#' Default is 0.8.
#' @return The output produced by the function is a list which consists of the local FDR and log10(Bayes factor).
#'    \item{locFDR}{It provides the analytically computed local false discovery rate (posterior probability of null association) under CPBayes model
#'     (a Bayesian analog of the p-value) which is a measure of the evidence of the
#'     aggregate-level pleiotropic association. Bayes factor is adjusted for prior odds, but
#'      locFDR is solely a function of the posterior odds.}
#'    \item{log10_BF}{It provides the analytically computed log10(Bayes factor) produced by CPBayes that measures the
#'     evidence of the overall pleiotropic association.}
#'
#' @references Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018) An efficient Bayesian meta analysis approach for studying cross-phenotype genetic associations. PLoS Genet 14(2): e1007139.
#'
#' @seealso \code{\link{cpbayes_cor}}, \code{\link{estimate_corln}}, \code{\link{analytic_locFDR_BF_uncor}}, \code{\link{cpbayes_uncor}}, \code{\link{post_summaries}}, \code{\link{forest_cpbayes}}
#'
#' @examples
#' data(ExampleDataCor)
#' BetaHat <- ExampleDataCor$BetaHat
#' BetaHat
#' SE <- ExampleDataCor$SE
#' SE
#' cor <- ExampleDataCor$cor
#' cor
#' result <- cpbayes_cor(BetaHat, SE, cor)
#' str(result)
#'
#' @export
analytic_locFDR_BF_cor = function(BetaHat, SE, Corln, SpikeVar=0.0001, SlabVar=0.8){
  # Check whether any of the primary arguments is missing
  if(missing(BetaHat) || missing (SE))
    stop("BetaHat or SE vector is missing!", call. = FALSE)
  if(missing(Corln))
    stop("Correlation matrix is missing!", call. = FALSE)
  # Argument 1 :: BetaHat
  BetaHat <- checkPrimaryVar(BetaHat, "BetaHat")
  # Argument 2 :: SE
  SE <- checkPrimaryVar(SE, "SE")
  # Check whether all entries are strictly positive
  if(!all(SE > 0))
    stop("One or more elements in the SE vector are not positive!", call. = FALSE)
  # Argument 1 and 2 ::
  if(length(BetaHat) != length(SE))
    stop("BetaHat and SE vectors must have the same number of elements!", call. = FALSE)

  # Argument 3 :: Correlation
  corln <- checkCorln(Corln, BetaHat)

  # Argument 4 and 5 ::
  chkVar = chkSlabVarSpikeVar(SpikeVar, SlabVar)
  spikevar = chkVar[["SpikeVar"]]
  slabvar = chkVar[["SlabVar"]]

  X = BetaHat
  K = length(X)    ## number of traits
  shapes = choose_shape_parameters(K, X, SE)                    ## get the shape params
  shape1 = shapes$shape1; shape2 = shapes$shape2;

  p = shape1/(shape1+shape2) ;
  logp = log(p); logq = log(1-p);

  ## If the covariance matrix is not positive definite, the diagonal elements are incremented to make it PD
  S <- diag(SE) %*% corln %*% diag(SE)
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



