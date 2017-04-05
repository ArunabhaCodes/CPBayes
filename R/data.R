# Documentation of the example data set to demonstrate how to run the uncorrelated version of CPBayes
#'An example data for uncorrelated summary statistics.
#'
#'ExampleDataUncor is a list which has two components: BetaHat, SE. The numeric vector 
#'ExampleDataUncor$BetaHat contains the main genetic effect (beta/log(odds ratio)) estimates
#'for a single nucleotide polymorphism (SNP) obtained from 10 separate case-control studies
#'for 10 different diseases. In each case-control study comprising a
#'distinct set of 7000 cases and 10000 controls, we fit a logistic regression
#'of the case-control status on the genotype coded as the minor allele count
#'for all the individuals in the sample. One can also include various covariates,
#'such as, age, gender, principal components (PCs) of ancestries in the logistic regression.
#'From each logistic regression for a
#'disease, we obtain the estimate of the main genetic association parameter
#'(beta/log(odds ratio)) along with the corresponding standard error.
#'Since the studies do not have any overlapping subject, the beta-hat across the traits are uncorrelated.
#'ExampleDataUncor$SE is the second numeric vector that contains the standard errors corresponding to
#'the uncorrelated beta-hat vector. 
#' @usage data(ExampleDataUncor)
#'  
#' @format A list of two numeric vectors each of length 10 (for 10 studies):
#' \describe{
#'     \item{BetaHat}{beta hat vector of length 10.}
#'     \item{SE}{standard error vector corresponding to beta-hat vector.}
#' }
#' 
#' @examples
#' data(ExampleDataUncor)
#' BetaHat <- ExampleDataUncor$BetaHat
#' BetaHat
#' SE <- ExampleDataUncor$SE
#' SE
#' \donttest{cpbayes_uncor(BetaHat, SE)}
"ExampleDataUncor"

#' An example data for correlated summary statistics.
#' 
#'ExampleDataCor is a list consisting of three components: BetaHat, SE, cor. ExampleDataCor$BetaHat is a
#'numeric vector that contains the main genetic effect (beta/log(odds ratio)) estimates
#'for a SNP across 10 overlapping case-control studies for 10 different diseases. Each of the 10 studies has
#'a distinct set of 7000 cases and a common set of 10000 controls shared across all the studies.
#'In each case-control study, we fit a logistic regression of the case-control status on the genotype
#'coded as the minor allele count for all the individuals in the sample. One can also include various
#'covariates, such as, age, gender, principal components (PCs) of ancestries in the logistic regression.
#'From each logistic regression for a disease, we obtain the estimate of the main genetic association
#'parameter (beta/log(odds ratio)) along with the corresponding standard error. Since the studies have overlapping
#'subjects, the beta-hat across traits are correlated. ExampleDataCor$SE contains the standard error vector
#'corresponding to the correlated beta-hat vector. ExampleDataCor$cor is a numeric square matrix providing
#'the correlation matrix of the correlated beta-hat vector. 
#' 
#' @usage data(ExampleDataCor)
#' 
#' @format A list consisting of two numeric vectors (each of length 10) and a numeric square matrix of
#'dimension 10 by 10: 
#'\describe{
#'  \item{BetaHat}{beta hat vector of length 10.}
#'  \item{SE}{standard error vector corresponding to the beta-hat vector.}
#'  \item{cor}{correlation matrix of the beta-hat vector.}
#'}
#' 
#' @examples 
#' data(ExampleDataCor)
#' BetaHat <- ExampleDataCor$BetaHat
#' BetaHat
#' SE <- ExampleDataCor$SE
#' SE
#' cor <- ExampleDataCor$cor
#' cor
#' \donttest{cpbayes_cor(BetaHat, SE, cor)}
"ExampleDataCor"

#' An example data of sample-overlap matrices.
#' 
#' An example data of sample-overlap matrices for five different diseases in the Kaiser 
#' GERA cohort (a real data).
#' SampleOverlapMatrix is a list that contains an example of the sample overlap matrices
#'  for five different diseases in the Kaiser GERA cohort. SampleOverlapMatrix$n11
#'   provides the number of cases shared between all possible pairs of diseases.
#'    SampleOverlapMatrix$n00 provides the number of controls shared between all possible
#'     pairs of diseases. SampleOverlapMatrix$n10 provides the number of subjects who
#'      are case for one disease and control for another disease.
#' 
#' @usage data(SampleOverlapMatrix)
#' 
#' @format A list consisting of three integer square matrices (each of dimension 5 by 5): 
#'\describe{
#'  \item{n11}{number of cases shared between all possible pairs of diseases.}
#'  \item{n00}{number of controls shared between all possible pairs of diseases.}
#'  \item{n10}{number of subjects who are case for one disease and control for another disease.}
#'}
#' @examples 
#' data(SampleOverlapMatrix)
#' n11 <- SampleOverlapMatrix$n11
#' n11
#' n00 <- SampleOverlapMatrix$n00
#' n00
#' n10 <- SampleOverlapMatrix$n10
#' n10
#'\donttest{estimate_corln(n11,n00,n10)}
"SampleOverlapMatrix"