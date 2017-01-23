
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

Simultaneous analysis of genetic associations with multiple phenotypes may reveal shared genetic susceptibility across traits (pleiotropy). CPBayes is a Bayesian meta analysis approach for studying cross-phenotype genetic associations. It uses summary-level data across multiple phenotypes to simultaneously measure the evidence of aggregate-level pleiotropic association and estimate an optimal subset of traits associated with the risk locus. CPBayes is based on a spike and slab prior and is implemented by Markov chain Monte Carlo (MCMC) technique Gibbs sampling.

This R-package consists of four main functions:

1.  cpbayes\_uncor: This function implements CPBayes for uncorrelated summary statistics. The summary statistics across traits/studies are uncorrelated when the studies have no overlapping subject.
2.  cpbayes\_cor: This function implements CPBayes for correlated summary statistics. The summary statistics across traits/studies are correlated when the studies have overlapping subjects or the phenotypes were measured in a cohort study.
3.  post\_summaries: This function summarizes the MCMC data produced by the two main functions cpbayes\_uncor or cpbayes\_cor listed above. It computes additional summaries to provide a better insight into a pleiotropic signal. It works in the same way for both cpbayes\_uncor and cpbayes\_cor.
4.  estimate\_corln: This function computes an approximate correlation matrix of the beta-hat vector for multiple overlapping case-control studies or a cohort study using the sample-overlap matrices.

Installation
============

You can install CPBayes from CRAN.

``` r
install.packages("CPBayes")
library("CPBayes")
```

An example demonstrating how to run CPBayes for uncorrelated summary statistics.
================================================================================

Get the path to the data.

``` r
library("CPBayes")
# Load the beta hat vector
BetaHatfile <- system.file("extdata", "BetaHat.rda", package = "CPBayes")
load(BetaHatfile)
BetaHat
```

BetaHat contains an example data of the main genetic effect (beta/log(odds ratio)) estimates for a single nucleotide polymorphism (SNP) obtained from 10 separate case-control studies for 10 different diseases. Since the studies do not have any overlapping subject, beta-hat across the diseases are uncorrelated.

``` r
# Load the standard error vector
SEfile <- system.file("extdata", "SE.rda", package = "CPBayes")
load(SEfile)
SE
```

SE contains the standard errors corresponding to the above beta hat vector across 10 separate case-control studies.

Next we specify the name of the diseases/phenotypes and the genetic variant.

``` r
# Specify the name of the traits and the genetic variant.
traitNames <- paste("Disease", 1:10, sep = "")
SNP1 <- "rs1234"
traitNames
SNP1
```

Now we implement CPBayes for this example data. Since the studies are non-overlapping, the summary statistics across traits are uncorrelated. Hence we run the the cpbayes\_uncor function.

``` r
# Run the uncorrelated version of CPBayes.
result <- cpbayes_uncor(BetaHat, SE, Phenotypes = traitNames, Variant = SNP1)
```

After running cpbayes\_uncor, it prints the log10(Bayes factor) (denoted as log10\_BF) and the subset of associated/non-null traits (denoted as subset) produced by CPBayes. The Bayes factor evaluates the overall pleiotropic association and the subset of non-null traits are the most important phenotypes that underlie the pleiotropic signal. However, the printed outputs are only a part of 'result' which is a list that constitutes of various components. An overall summary of 'result' can be seen by using the str function (as shown below).

``` r
# Overall summary of the primary results produced by cpbayes_uncor.
str(result)
```

A detailed interpretation of all the outputs are described in the Value section of cpbayes\_uncor in the CPBayes manual.

The post\_summaries function provides important insights into an obseved pleiotropic signal, e.g., the direction of associations, posterior mean/median and 95% credible interval (Bayesian analog of the confidence interval) of the unknown true genetic effect (beta/odds ratio) on each trait, etc.

``` r
# Post summary of the MCMC data produced by cpbayes_uncor.
PleioSumm <- post_summaries(result, level = 0.05)  
str(PleioSumm)
```

So we have to pass the list 'result' returned by cpbayes\_uncor as the first argument and the 'level' as the second argument into the post\_summaries function. If 'level' is not specified, the default value is 0.05. For detailed description of different outputs provided by this function, please see the Value section of post\_summaries in the CPBayes manual.

An example demonstrating how to run CPBayes for correlated summary statistics.
==============================================================================

Next we demonstrate how to run CPBayes for correlated summary statistics. Get the path to the data.

``` r
# Load the beta-hat vector
datafile <- system.file("extdata", "cBetaHat.rda", package = "CPBayes")
load(datafile)
cBetaHat
```

Here 'c' in cBetaHat stands for correlated case. cBetaHat contains an example data of the main genetic association parameter (beta/log(odds ratio)) estimates for a SNP across 10 overlapping case-control studies for 10 different diseases. Each of the 10 studies has a distinct set of 7000 cases and a common set of 10000 controls shared across all the studies. Since the studies have overlapping subjects, beta-hat across the diseases are correlated.

``` r
# Load the standard error vector
datafile <- system.file("extdata", "cSE.rda", package = "CPBayes")
load(datafile)
cSE
```

cSE contains the standard errors corresponding to the above beta hat vector across 10 overlapping case-control studies.

``` r
# Load the correlation matrix of the beta-hat vector (cBetaHat)
datafile <- system.file("extdata", "cor.rda", package = "CPBayes")
load(datafile)
cor
```

The correlation matrix of the beta-hat vector (cBetaHat) is given by 'cor' which we estimated by employing the estimate\_corln function (demonstrated later in this tutorial) using the sample-overlap matrix (explained later in this tutorial). Next we run the correlated version of CPBayes for this example data.

``` r
# Run the correlated version of CPBayes.
result <- cpbayes_cor(cBetaHat, cSE, cor, Phenotypes = traitNames, Variant = SNP1)
```

After running cpbayes\_cor, it prints the log10(Bayes factor) (denoted as log10\_BF) and the subset of non-null/associated traits (denoted as subset) produced by CPBayes. However, the printed outputs are only a part of 'result' which is a list that constitutes of various components. An overall summary of 'result' can be seen by using the str function (as shown below).

``` r
# Overall summary of the primary results produced by cpbayes_cor.
str(result)
```

A detailed interpretation of all the outputs are described in the Value section of cpbayes\_cor in the CPBayes manual.

The post\_summaries function provides important insights into an observed pleiotropic signal, e.g., the direction of associations, posterior mean/median and 95% credible interval (Bayesian analog of the confidence interval) of the unknown true genetic effect (beta/odds ratio) on each trait, etc.

``` r
# Post summary of the MCMC data produced by cpbayes_cor.
PleioSumm <- post_summaries(result, level = 0.05)  
str(PleioSumm)
```

post\_summaries works exactly in the same way for both cpbayes\_cor and cpbayes\_uncor. Thus, we need to pass the list 'result' returned by cpbayes\_cor as the first argument and 'level' as the second argument into the post\_summaries function. For detailed description of different outputs provided by post\_summaries, please see the Value section of post\_summaries in the CPBayes manual.

An example how to run estimate\_corln
=====================================

The function estimate\_corln estimates the correlation matrix of the beta-hat vector for multiple overlapping case-control studies or a cohort study using the sample-overlap matrices which describe the number of cases or controls shared between studies/traits, and the number of subjects who are case for one study/trait but control for another study/trait.

``` r
# Example data of sample-overlap matrices
SampleOverlapMatrixFile <- system.file("extdata", "SampleOverlapMatrix.rda", package = "CPBayes")
load(SampleOverlapMatrixFile)
SampleOverlapMatrix
```

SampleOverlapMatrix is a list that contains an example of the sample overlap matrices for five different diseases in the Kaiser GERA cohort (a real data). The list constitutes of three matrices as follows. SampleOverlapMatrix$n11 provides the number of cases shared between all possible pairs of studies/traits. SampleOverlapMatrix$n00 provides the number of controls shared between all possible pairs of studies/traits. SampleOverlapMatrix$n10 provides the number of subjects who are case for one study/trait and control for another study/trait. For more detailed explanation, see the Arguments section of estimate\_corln in the CPBayes manual.

``` r
# Estimate the correlation matrix of correlated beta-hat vector
n11 <- SampleOverlapMatrix$n11
n00 <- SampleOverlapMatrix$n00
n10 <- SampleOverlapMatrix$n10
cor <- estimate_corln(n11, n00, n10)
cor
```

The function estimate\_corln computes an approximate correlation matrix of the correlated beta-hat vector obtained from multiple overlapping case-control studies or a cohort study. While demonstrating cpbayes\_cor, we used simulated data for 10 overlapping case-control studies with each study having a distinct set of 7000 cases and a common set of 10000 controls shared across all the studies. We used the estimate\_corln function to estimate the correlation matrix of the correlated beta-hat vector using the sample-overlap matrix.

Getting more details
====================

Please see the vignettes and the manual of CPBayes for more detailed description about the package. If you encounter any problem/issue related to the CPBayes package, please email to <statgen.arunabha@gmail.com> and <tanushree.haldar@gmail.com>. Also, see our paper for more details: Arunabha Majumdar, Tanushree Haldar, Sourabh Bhattacharya, John Witte. An efficient Bayesian meta-analysis approach for studying cross-phenotype genetic associations (submitted). Available at: <http://biorxiv.org/content/early/2017/01/18/101543>.
