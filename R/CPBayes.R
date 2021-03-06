#' CPBayes: An R-package implemeting a Bayesian meta analysis method for studying cross-phenotype
#' genetic associations.
#'
#' Simultaneous analysis of genetic associations with multiple phenotypes may reveal shared
#' genetic susceptibility across traits (pleiotropy). CPBayes is a Bayesian meta analysis
#' method for studying cross-phenotype genetic associations. It uses summary-level data
#' across multiple phenotypes to simultaneously measure the evidence of aggregate-level
#' pleiotropic association and estimate an optimal subset of traits associated with the
#' risk locus. CPBayes is based on a spike and slab prior.
#'
#' The package consists of following functions:
#'\code{\link{analytic_locFDR_BF_uncor}}, \code{\link{cpbayes_uncor}}; \code{\link{analytic_locFDR_BF_cor}}, \code{\link{cpbayes_cor}}; \code{\link{post_summaries}}, \code{\link{forest_cpbayes}}, \code{\link{estimate_corln}}.
#'
#' @section Functions:
#' \describe{
#' \item{\code{\link{analytic_locFDR_BF_uncor}}}{It analytically computes the local FDR (locFDR)
#' and Bayes factor (BF) quantifying the evidence of aggregate-level pleiotropic association
#' for uncorrelated summary statistics.}
#' \item{\code{\link{cpbayes_uncor}}}{It implements CPBayes (based on MCMC) for uncorrelated summary statistics to
#' figure out the optimal subset of non-null traits underlying a pleiotropic signal and other insights.
#' The summary statistics across traits/studies are uncorrelated when the studies
#' have no overlapping/genetically related subjects.}
#' \item{\code{\link{analytic_locFDR_BF_cor}}}{It analytically computes the local FDR (locFDR)
#'  and Bayes factor (BF) for correlated summary statistics.}
#' \item{\code{\link{cpbayes_cor}}}{It implements CPBayes (based on MCMC) for correlated summary statistics to figure out
#' the optimal subset of non-null traits underlying a pleiotropic signal and other insights. The summary statistics across
#' traits/studies are correlated when the studies have overlapping/genetically related subjects
#' or the phenotypes were measured in a cohort study.}
#' \item{\code{\link{post_summaries}}}{It summarizes the MCMC data produced by
#'  \code{\link{cpbayes_uncor}} or \code{\link{cpbayes_cor}}.
#'   It computes additional summaries to provide a better insight into a pleiotropic signal.
#'    It works in the same way for both \code{\link{cpbayes_uncor}} and \code{\link{cpbayes_cor}}.}
#' \item{\code{\link{forest_cpbayes}}}{It creates a forest plot presenting the pleiotropy result obtained by
#' \code{\link{cpbayes_uncor}} or \code{\link{cpbayes_cor}}. It works in the same way for
#'  both \code{\link{cpbayes_uncor}} and \code{\link{cpbayes_cor}}.}
#' \item{\code{\link{estimate_corln}}}{It computes an approximate correlation matrix of
#'  the beta-hat vector for multiple overlapping case-control studies using the
#'  sample-overlap count matrices.}
#' }
#' @references Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018) An efficient Bayesian meta analysis approach for studying cross-phenotype genetic associations. PLoS Genet 14(2): e1007139.
#'
#'
#' @docType package
#'
#' @name CPBayes
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom purrr map_dbl
#' @importFrom MASS mvrnorm
#' @importFrom utils combn
#' @importFrom stats runif rnorm rbeta quantile qchisq qbeta pchisq pbeta p.adjust dist aggregate sd dnorm qnorm
#' @importFrom forestplot forestplot fpColors
#' @importFrom grDevices dev.off pdf
NULL
