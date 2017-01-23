## Function to estimate correlation using case-control sample overlap matrix
#' Estimate correlation structure of beta-hat vector for multiple overlapping case-control studies
#' or a cohort study using sample-overlap matrix.
#'
#' Compute an approximate correlation matrix of the beta-hat vector for multiple overlapping
#'  case-control studies or a cohort study using the sample-overlap matrices which describe
#' the number of cases or controls shared between studies/traits, and the number of subjects
#' who are case for one study/trait but control for another study/trait. This approximation is more accurate
#' when none of the diseases/traits is associated with the environmental covariates present in the study.
#'   
#'***Important note on the estimation of correlation structure of correlated beta-hat vector:***
#' In general, environmental covariates are expected to be present in a study and associated
#' with the phenotypes of interest. Hence the above approximation of the correlation matrix
#' may not be completely accurate. So, in presence of environmental covariates, we recommend
#' an alternative strategy to estimate the correlation matrix using the genome-wide summary
#' statistics data across traits as follows. First, extract all the SNPs for each of which the
#' trait-specific univariate association p-value across all the traits are > 0.1. The
#' trait-specific univariate association p-values can be obtained based on the beta-hat
#' and standard error for each trait. Each of the SNPs selected in this way is either weakly
#' or not associated with any of the phenotypes (null SNP). Next, select a set of independent
#' null SNPs from the initial set of null SNPs by using a threshold of r^2 < 0.01 (r: the 
#' correlation between the genotypes at a pair of SNPs). Finally, compute the correlation
#' matrix of the effect estimates (beta-hat vector) as the sample correlation matrix of the
#' beta-hat vector across all the selected independent null SNPs. This strategy is more
#' general and applicable to a cohort study or multiple overlapping studies for binary or
#' quantitative traits with arbitrary distributions. Misspecification of the correlation
#' structure can affect the results produced by CPBayes to some extent. Hence, if 
#' genome-wide summary statistics data across traits is available, we recommend to use
#' this alternative strategy to estimate the correlation matrix of the beta-hat vector.
#' See our paper for more details at:  http://biorxiv.org/content/early/2017/01/18/101543. 
#'  
#' @param n11 An integer square matrix (number of rows must be the same as the number of
#' studies/traits) providing the
#' number of cases shared between all possible pairs of studies/traits. So (k,l)-th element of n11
#' is the number of subjects who are case for both k-th and l-th study/trait. Note that the diagonal elements of
#' n11 are the number of cases across studies/traits. In case, no case is shared between studies/traits,
#' the off-diagonal elements of n11 will be zero. No default is specified.
#' @param n00 An integer square matrix (number of rows must be the same as the
#' number of studies/traits) providing the
#' number of controls shared between all possible pairs of studies/traits. So (k,l)-th element of n00
#' is the number subjects who are control for both k-th and l-th study/trait. Note that the diagonal
#' elements of n00 are the number of controls across studies/traits. In case, no control is
#'  shared between studies/traits,
#' the off-diagonal elements will be zero. No default is specified.
#' @param n10 An integer square matrix (number of rows must be the same as the
#'  number of studies/traits) providing the
#' number of subjects who are case for one study/trait and control for another study/trait.
#'  Clearly, the diagonal elements 
#' will be zero. An off diagonal element, e.g., (k,l)-th element of n10 is the number of subjects who
#' are case for k-th study/trait and control for l-th study/trait. If there is no such overlap,
#' all the elements
#' of n10 will be zero. No default is specified.
#' @return This function returns an approximate correlation matrix of the beta-hat vector for 
#' multiple overlapping case-control studies or a cohort study. See the example below.
#'
#' @references Arunabha Majumdar, Tanushree Haldar, Sourabh Bhattacharya, John Witte.
#'  An efficient Bayesian meta-analysis 
#'  approach for studying cross-phenotype genetic associations (submitted). Available
#'  at: http://biorxiv.org/content/early/2017/01/18/101543.
#'  
#' @seealso \code{\link{cpbayes_cor}}
#' 
#' @examples
#' data(SampleOverlapMatrix)
#' n11 <- SampleOverlapMatrix$n11
#' n11
#' n00 <- SampleOverlapMatrix$n00
#' n00
#' n10 <- SampleOverlapMatrix$n10
#' n10
#' cor <- estimate_corln(n11, n00, n10)
#' cor
#' 
#' @export
estimate_corln <- function(n11, n00, n10)
{
  ## Checking of the three matrices n11, n00, n10
  ## checking of n11 - diagonal elements > 0, (symmetric) 
  ## checking of n00 - diagonal elements > 0  (symmetric)
  ## checking of n10 - diagonal elements = 0
  ## each matrix will be integer square matrix, of same dimension
  
  if(missing(n11) || missing(n00) || missing(n10))
    stop("n11, n00 or n10 matrix is missing!", call. = FALSE)
  else
    chkEstCorln(n11, n00, n10)
  
  
  
  
  

  
  
  
  
    
  
  n01 <- t(n10)                                          #n01 = transpose(n10)
  
  n1 <- diag(n11)                                        #number of cases for different traits
  n0 <- diag(n00)                                        #number of control for different traits
  n <- n1+n0                                             #total sample size of different studies
  sqrt_n <- sqrt(n)                                      #square root of sample size of different studies
  sqrt_samp_size_prod <- sqrt_n%*%t(sqrt_n)              #sqrt of sample size product matrix
  sqrt_n0 <- sqrt(n0)                                    #sqrt of number of controls across studies
  t_sqrt_n0 <- t(sqrt_n0)                                #transpose of sqrt of number of controls column vector
  sqrt_n1 <- sqrt(n1)                                    #sqrt of number of cases across studies
  t_sqrt_n1 <- t(sqrt_n1)                                #transpose of number of cases column vector
  
  sqrt_n0_n0_prod <- sqrt_n0 %*% t_sqrt_n0
  sqrt_n1_n1_prod <- sqrt_n1 %*% t_sqrt_n1
  sqrt_n1_n0_prod <- sqrt_n1 %*% t_sqrt_n0
  sqrt_n0_n1_prod <- sqrt_n0 %*% t_sqrt_n1
  
  a1 <- n11/sqrt_samp_size_prod                          #1st component of 1st part of correlation
  a2 <- sqrt_n0_n0_prod/sqrt_n1_n1_prod                  #2nd component of 1st part of correlation
  a <- a1*a2
  
  b1 <- n10/sqrt_samp_size_prod                          #1st component of 2nd part of correlation
  b2 <- sqrt_n0_n1_prod/sqrt_n1_n0_prod                  #2nd component of 2nd part of correlation
  b <- b1*b2
  
  c1 <- n01/sqrt_samp_size_prod                          #1st component of 3rd part of correlation
  c2 <- sqrt_n1_n0_prod/sqrt_n0_n1_prod                  #2nd component of 3rd part of correlation
  c <- c1*c2
  
  d1 <- n00/sqrt_samp_size_prod                          #1st component of 4th part of correlation
  d2 <- sqrt_n1_n1_prod/sqrt_n0_n0_prod                  #2nd component of 4th part of correlation
  d <- d1*d2
  
  corr <- a-b-c+d
  diag(corr) <- rep(1,nrow(n11))
  
  return(corr)
}
