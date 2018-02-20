## Function to call CPBayes_cor (CPBayes function for correlated summary statistics)
#' Run correlated version of CPBayes.
#' 
#' Run correlated version of CPBayes when the main genetic effect (beta/log(odds ratio)) estimates across
#'  studies/traits are correlated.
#' @param BetaHat A numeric vector of length K where K is the number of phenotypes. It
#'  contains the beta-hat values across studies/traits. No default is specified.
#' @param SE A numeric vector with the same dimension as BetaHat providing the standard errors
#'  corresponding to BetaHat. Every element of SE must be positive. No default is specified.
#' @param Corln A numeric square matrix of order K by K providing the correlation matrix of BetaHat.
#'  The number of rows & columns of Corln must be the same as the length of BetaHat. No default
#'   is specified. See \code{\link{estimate_corln}}.
#' @param Phenotypes A character vector of the same length as BetaHat providing the name
#'  of the phenotypes. Default is specified as trait1, trait2, . . . , traitK. Note that
#'   BetaHat, SE, Corln, and Phenotypes must be in the same order.
#' @param Variant A character vector of length one providing the name of the genetic variant.
#'  Default is `Variant'.
#' @param UpdateSlabVar A logical vector of length one. If TRUE, the variance of the slab distribution
#'  that presents the prior distribution of non-null effects is
#'  updated at each MCMC iteration in a range (MinSlabVar -- MaxSlabVar) (see next). If FALSE, 
#'   it is fixed at (MinSlabVar + MaxSlabVar)/2. Default is TRUE.
#' @param MinSlabVar A numeric value greater than 0.01 providing the minimum value of the
#'  variance of the slab distribution. Default is 0.6.
#' @param MaxSlabVar A numeric value smaller than 10.0 providing the maximum value of the
#'  variance of the slab distribution. Default is 1.0. **Note that,
#'  a smaller value of the slab variance will increase the sensitivity of CPBayes while selecting the optimal
#'   subset of associated traits but at the expense of lower specificity. Hence the slab variance
#'   parameter in CPBayes is inversely related to the level of false discovery rate (FDR) in a frequentist
#'   FDR controlling procedure. For a specific dataset, an user
#'    can experiment different choices of these three arguments: UpdateSlabVar, MinSlabVar, and MaxSlabVar.
#' @param MCMCiter A positive integer greater than or equal to 10,000 providing the total number of
#'  iterations in the MCMC. Default is 20,000.
#' @param Burnin A positive integer greater than or equal to 2,000 providing the burn in period 
#' in the MCMC. Default is 5,000. Note that the MCMC sample size (MCMCiter - Burnin) must be at least 5,000.
#' @return The output produced by \code{\link{cpbayes_cor}} is a list which consists of various components.
#'    \item{variantName}{It is the name of the genetic variant provided by the user. If 
#'    not specified by the user, default name is `Variant'.} 
#'    \item{log10_BF}{It provides the log10(Bayes factor) produced by CPBayes that measures the
#'     evidence of the overall pleiotropic association.}
#'    \item{locFDR}{It provides the posterior probability of null association produced by
#'     CPBayes (a Bayesian analog of the p-value) which is another measure of the evidence
#'      of aggregate-level pleiotropic association. Bayes factor is adjusted for prior odds, but
#'      locFDR (PPNA) is solely a function of the posterior odds. locFDR can sometimes be small
#'      indicating an association, but log10_BF may not indicate an association. Hence, always check both log10_BF and locFDR.}
#'    \item{subset}{It provides the optimal subset of associated/non-null traits selected by
#'     CPBayes. It is NULL if no phenotype is selected.}
#'    \item{important_traits}{It provides the traits which yield a trait-specific posterior 
#'    probability of association (PPAj) > 20\%. Even if a phenotype is not 
#'    selected in the optimal subset of non-null traits, it can produce a non-negligible
#'     value of PPAj. Note that, `important_traits' is expected to include the traits already
#'      contained in `subset'. It provides both the name of the important traits and their
#'       corresponding value of PPAj. Always check 'important_traits' even if 'subset' contains
#'        a single trait. It helps to better explain an observed pleiotropic signal.}
#'    \item{auxi_data}{It contains supplementary data including the MCMC data which is used later 
#'     by \code{\link{post_summaries}} and \code{\link{forest_cpbayes}}:
#'        \enumerate{
#'            \item traitNames: Name of all the phenotypes.
#'            \item K: Total number of phenotypes.
#'            \item mcmc.samplesize: MCMC sample size.
#'            \item PPAj: Trait-specific posterior probability of association for all the traits. 
#'            \item Z.data: MCMC data on the latent association status of all the traits (Z).
#'            \item sim.beta: MCMC data on the unknown true genetic effect (beta) on each trait.
#'            \item betahat: The beta-hat vector provided by the user which will be used by \code{\link{forest_cpbayes}}.
#'            \item se: The standard error vector provided by the user which will be used by \code{\link{forest_cpbayes}}.
#'        }
#'    }
#'    \item{uncor_use}{'Yes' or 'No'. Whether the combined strategy of CPBayes (implemented for correlated
#'     summary statistics) used the uncorrelated version or not.}
#'    \item{runtime}{It provides the runtime (in seconds) taken by \code{\link{cpbayes_cor}}. It will help the user to
#'     plan the whole analysis.}
#'    
#' @references Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018) An efficient Bayesian meta analysis approach for studying cross-phenotype genetic associations. PLoS Genet 14(2): e1007139.
#'    
#' @seealso  \code{\link{estimate_corln}}, \code{\link{post_summaries}}, \code{\link{forest_cpbayes}}, \code{\link{cpbayes_uncor}}
#' 
#' @examples
#' data(ExampleDataCor)
#' BetaHat <- ExampleDataCor$BetaHat
#' BetaHat
#' SE <- ExampleDataCor$SE
#' SE
#' cor <- ExampleDataCor$cor
#' cor
#' traitNames <- paste("Disease", 1:10, sep = "")
#' SNP1 <- "rs1234"
#' result <- cpbayes_cor(BetaHat, SE, cor, Phenotypes = traitNames, Variant = SNP1)
#' str(result)
#' 
#' 
#' @export
cpbayes_cor <- function(BetaHat, SE, Corln, Phenotypes, Variant, UpdateSlabVar = TRUE, MinSlabVar = 0.6, MaxSlabVar = 1.0, MCMCiter = 20000, Burnin = 5000)
{
  
  UpdateDE <- UpdateSlabVar
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
    COR <- checkCorln(Corln, BetaHat)
  
  # Argument 4 :: Phenotype names
  if(!missing(Phenotypes))
    checkPhen(Phenotypes, BetaHat)
  else Phenotypes = paste("trait", 1:length(BetaHat), sep = "")
  
  # Argument 5 :: Variant name
  if(!missing(Variant))
  {
    Variant <- checkVarName(Variant)
    variantName <- unname(Variant)           # Assignment
  }
  else variantName <- "Variant"
  
    # Argument 6 :: Update model parameter DE
    # Check whether argument 6 is a vector of length 1
    if(!is.vector(UpdateDE) || (length(UpdateDE) != 1))
    {
      warning("UpdateSlabVar is not a vector of length 1 (default option used).", call. = FALSE)
      UpdateDE <- TRUE
    }
    if(!is.logical(UpdateDE))
    {
      warning("UpdateSlabVar not provided as logical (default option used).", call. = FALSE)
      UpdateDE <- TRUE
    }

    # Argument 7 & 8:: Minimum and maximum value of slab variance
    MinSlabVarDefault <- 0.6
    MaxSlabVarDefault <- 1.0
    MinSlabVarBound <- 0.01
    MaxSlabVarBound <- 10.0
    
    # Check whether argument 7 is a vector of length 1
    if(!is.vector(MinSlabVar) || (length(MinSlabVar) != 1))
    {
      warning("MinSlabVar is not a vector of length 1 (default option used).", call. = FALSE)
      MinSlabVar <- MinSlabVarDefault
    }
    # Check whether argument 7 is numeric
    if(!is.numeric(MinSlabVar)){
      warning("MinSlabVar is not numeric (default option used).", call. = FALSE)
      MinSlabVar <- MinSlabVarDefault
    }
    # Check whether argument 7 is more than its minimum bound
    if(MinSlabVar < MinSlabVarBound){
      warning("MinSlabVar should not be smaller than ", MinSlabVarBound, ", so it is assigned to ", MinSlabVarBound, ".", call. = FALSE)
      MinSlabVar <- MinSlabVarBound
    }
    
    # Argument 8 :: Maximum value of slab variance    
    
    # Check whether argument 8 is a vector of length 1
    if(!is.vector(MaxSlabVar) || (length(MaxSlabVar) != 1))
    {
      warning("MaxSlabVar is not a vector of length 1 (default option used).", call. = FALSE)
      MaxSlabVar <- MaxSlabVarDefault
    }
    # Check whether argument 8 is numeric
    if(!is.numeric(MaxSlabVar)){
      warning("MaxSlabVar is not numeric (default option used).", call. = FALSE)
      MaxSlabVar <- MaxSlabVarDefault
    }
    # Check whether argument 8 is less than a maximum bound
    if(MaxSlabVar > MaxSlabVarBound){
      warning("MaxSlabVar should not be bigger than ", MaxSlabVarBound, ", so it is assigned to ", MaxSlabVarBound, ".", call. = FALSE)
      MaxSlabVar <- MaxSlabVarBound
    }
    
    # Argument 7 and 8 :: Checking MinSlabVar < MaxSlabVar
    if(MinSlabVar >= MaxSlabVar){
      warning("MaxSlabVar is not bigger than MinSlabVar! (default option used).", call. = FALSE)
      MinSlabVar <- MinSlabVarDefault
      MaxSlabVar <- MaxSlabVarDefault
    }
    

    # Argument 9 :: Number of MCMC iteration
    mcmcDefault <- 20000
    # Check whether argument 9 is a vector of length 1
    if(!is.vector(MCMCiter) || (length(MCMCiter) != 1))
    {
      warning("MCMCiter is not a vector of length 1 (default option used).", call. = FALSE)
      MCMCiter <- mcmcDefault
    }
    # Check whether argument 9 is numeric and integer
    if(!is.numeric(MCMCiter) || MCMCiter%%1 != 0)
    {
      warning("MCMCiter not provided as integer (default option used).", call. = FALSE)
      MCMCiter <- mcmcDefault
    }
    # Check whether argument 9 is more than 10000
    if(MCMCiter < 10000)
    {
      warning("MCMCiter should be at least 10000 (default option used).", call. = FALSE)
      MCMCiter <- mcmcDefault
    }
    
    # Argument 10 :: Burnin
    BurninDefault <- 5000
    # Check whether argument 10 is a vector of length 1
    if(!is.vector(Burnin) || (length(Burnin) != 1))
    {
      warning("Burnin is not a vector of length 1 (default option used).", call. = FALSE)
      Burnin <- BurninDefault
    }
    # Check whether argument 10 is numeric and integer
    if(!is.numeric(Burnin) || Burnin%%1 != 0)
    {
      warning("Burnin not provided as integer (default option used).", call. = FALSE)
      Burnin <- BurninDefault
    }
    # Check whether argument 10 is more than 2000
    if(Burnin < 2000)
    {
      warning("Burnin should be at least 2000 (default option used).", call. = FALSE)
      Burnin <- BurninDefault
    }
    
    # Argument 9 and 10
    if((MCMCiter - Burnin) < 5000)
    {
      warning("MCMC sample size (MCMCiter - Burnin) provided less than 5000 (default options used)", call. = FALSE)
      MCMCiter <- mcmcDefault
      Burnin <- BurninDefault
    }
    
    # Call CPBayes function
    RESULT <- combined_CPBayes( variantName, Phenotypes, BetaHat, SE, COR, UpdateDE, MinSlabVar, MaxSlabVar, MCMCiter, Burnin )
    print_result(RESULT)
    invisible(RESULT)
}



