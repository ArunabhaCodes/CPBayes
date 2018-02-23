
<!-- README.md is generated from README.Rmd. Please edit that file -->
Overview
========

Simultaneous analysis of genetic associations with multiple phenotypes may reveal shared genetic susceptibility across traits (pleiotropy). CPBayes is a Bayesian meta analysis method for studying cross-phenotype genetic associations. It uses summary-level data across multiple phenotypes to simultaneously measure the evidence of aggregate-level pleiotropic association and estimate an optimal subset of traits associated with the risk locus. CPBayes is based on a spike and slab prior and is implemented by Markov chain Monte Carlo (MCMC) technique Gibbs sampling.

This R-package consists of five main functions:

1.  cpbayes\_uncor: It implements CPBayes for uncorrelated summary statistics. The summary statistics across traits/studies are uncorrelated when the studies have no overlapping subject.
2.  cpbayes\_cor: It implements CPBayes for correlated summary statistics. The summary statistics across traits/studies are correlated when the studies have overlapping subjects or the phenotypes were measured in a cohort study.
3.  post\_summaries: It summarizes the MCMC data produced by cpbayes\_uncor or cpbayes\_cor. It computes additional summaries to provide a better insight into a pleiotropic signal. It works in the same way for both cpbayes\_uncor and cpbayes\_cor.
4.  forest\_cpbayes: It creates a forest plot presenting the pleiotropy result obtained by cpbayes\_uncor or cpbayes\_cor. It works in the same way for both cpbayes\_uncor and cpbayes\_cor.
5.  estimate\_corln: It computes an approximate correlation matrix of the beta-hat vector for multiple overlapping case-control studies using the sample-overlap matrices.

Installation
============

You can install CPBayes from CRAN.

``` r
install.packages("CPBayes")
library("CPBayes")
```

How to run CPBayes for uncorrelated summary statistics.
=======================================================

Get the path to the data.

``` r
library("CPBayes")
# Load the beta hat vector
BetaHatfile <- system.file("extdata", "BetaHat.rda", package = "CPBayes")
load(BetaHatfile)
BetaHat
```

BetaHat contains an example data of the main genetic effect (beta/log(odds ratio)) estimates for a single nucleotide polymorphism (SNP) obtained from 10 separate case-control studies for 10 different diseases. Since the studies do not have any overlapping subject, beta-hat across the diseases can be assumed uncorrelated.

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

Now we implement CPBayes for this example data. Since the studies are non-overlapping, we run the the cpbayes\_uncor function.

``` r
# Run the uncorrelated version of CPBayes.
result <- cpbayes_uncor(BetaHat, SE, Phenotypes = traitNames, Variant = SNP1)
```

After running cpbayes\_uncor, it prints the local false discovery rate (denoted as locFDR) and the subset of associated/non-null traits (denoted as subset) produced by CPBayes. The locFDR evaluates the overall pleiotropic association and the subset of non-null traits are the most important phenotypes that underlie the pleiotropic signal. However, the printed outputs are only a part of 'result' which is a list that constitutes of various components. An overall summary of 'result' can be seen by using the str() function (as shown below).

``` r
# Overall summary of the primary results produced by cpbayes_uncor.
str(result)
```

A detailed interpretation of all the outputs are described in the Value section of cpbayes\_uncor in the CPBayes manual.

The post\_summaries function provides important insights into an obseved pleiotropic signal, e.g., the direction of associations, trait-specific posterior probability of associations (PPAj), posterior mean/median and 95% credible interval (Bayesian analog of the confidence interval) of the unknown true genetic effect (beta/odds ratio) on each trait, etc.

``` r
# Post summary of the MCMC data produced by cpbayes_uncor.
PleioSumm <- post_summaries(result, level = 0.05)  
str(PleioSumm)
```

So we have to pass the list 'result' returned by cpbayes\_uncor as the first argument and the 'level' as the second argument into the post\_summaries function. If 'level' is not specified, the default value is 0.05. For detailed description of different outputs provided by this function, see the Value section of post\_summaries in the CPBayes manual.

Next we run the forest\_cpbayes function to create a forest plot that presents the pleiotropy result produced by cpbayes\_uncor.

``` r
# Forest plot for the pleiotropy result obtained by cpbayes_uncor.
forest_cpbayes(result, level = 0.05)
```

Similarly as for the post\_summaries function, we need to pass the same list \`result' returned by cpbayes\_uncor as the first argument into the function. Second argument is the level whose default value is 0.05. In the forest plot, (1-level)% confidence interval of the beta/log odds ratio parameter is plotted for each trait. For more details, see the section of forest\_cpbayes function in the CPBayes manual.

How to run CPBayes for correlated summary statistics.
=====================================================

Next we demonstrate how to run CPBayes for correlated summary statistics. Get the path to the data.

``` r
# Load the beta-hat vector
datafile <- system.file("extdata", "cBetaHat.rda", package = "CPBayes")
load(datafile)
cBetaHat
```

Here 'c' in cBetaHat stands for correlated case. cBetaHat contains an example data of the main genetic association parameter (beta/log odds ratio) estimates for a SNP across 10 overlapping case-control studies for 10 different diseases. Each of the 10 studies has a distinct set of 7000 cases and a common set of 10000 controls shared across all the studies. Since the studies have overlapping subjects, beta-hat across the diseases are correlated.

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

The correlation matrix of the beta-hat vector (cBetaHat) is given by 'cor' which we estimated by employing the estimate\_corln function (demonstrated later in this tutorial) using the sample-overlap matrices (explained later in this tutorial). Next we run the correlated version of CPBayes for this example data.

``` r
# Run the correlated version of CPBayes.
result <- cpbayes_cor(cBetaHat, cSE, cor, Phenotypes = traitNames, Variant = SNP1)
```

After running cpbayes\_cor, it prints the local false discovery rate (denoted as locFDR) and the subset of non-null/associated traits (denoted as subset) produced by CPBayes. However, the printed outputs are only a part of 'result' which is a list that constitutes of various components. An overall summary of 'result' can be seen by using the str() function (as shown below).

``` r
# Overall summary of the primary results produced by cpbayes_cor.
str(result)
```

A detailed interpretation of all the outputs are described in the Value section of cpbayes\_cor in the CPBayes manual.

The post\_summaries function provides important insights into an observed pleiotropic signal, e.g., the direction of associations, trait-specific posterior probability of associations (PPAj), posterior mean/median and 95% credible interval (Bayesian analog of the confidence interval) of the unknown true genetic effect (beta/odds ratio) on each trait, etc.

``` r
# Post summary of the MCMC data produced by cpbayes_cor.
PleioSumm <- post_summaries(result, level = 0.05)  
str(PleioSumm)
```

post\_summaries works exactly in the same way for both cpbayes\_cor and cpbayes\_uncor. For detailed description of different outputs provided by post\_summaries, see the Value section of post\_summaries in the CPBayes manual.

Next we run the forest\_cpbayes function to create a forest plot that presents the pleiotropy result produced by cpbayes\_cor.

``` r
# Forest plot for the pleiotropy result obtained by cpbayes_cor.
forest_cpbayes(result, level = 0.05)
```

Note that, forest\_cpbayes works exactly in the same way for both cpbayes\_cor and cpbayes\_uncor. For more details, see the section of forest\_cpbayes function in the CPBayes manual.

How to run estimate\_corln.
===========================

The function estimate\_corln estimates the correlation matrix of the beta-hat vector for multiple overlapping case-control studies using the sample-overlap matrices which describe the number of cases or controls shared between studies/traits, and the number of subjects who are case for one study/trait but control for another study/trait. For a cohort study, the phenotypic correlation matrix should be a reasonable substitute of this correlation matrix.

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

The function estimate\_corln computes an approximate correlation matrix of the correlated beta-hat vector obtained from multiple overlapping case-control studies using the sample-overlap matrices. Note that for a cohort study, the phenotypic correlation matrix should be a reasonable substitute of this correlation matrix. These approximations of the correlation structure are accurate when none of the diseases/traits is associated with the environmental covariates and genetic variant. While demonstrating cpbayes\_cor, we used simulated data for 10 overlapping case-control studies with each study having a distinct set of 7000 cases and a common set of 10000 controls shared across all the studies. We used the estimate\_corln function to estimate the correlation matrix of the correlated beta-hat vector using the sample-overlap matrices.

***Important note on the estimation of correlation structure of correlated beta-hat vector:*** In general, environmental covariates are expected to be present in a study and associated with the phenotypes of interest. Also, a small proportion of genome-wide genetic variants are expected to be associated. Hence the above approximations of the correlation matrix may not be accurate. So in general, we recommend an alternative strategy to estimate the correlation matrix using the genome-wide summary statistics data across traits as follows. First, extract all the SNPs for each of which the trait-specific univariate association p-value across all the traits are &gt; 0.1. The trait-specific univariate association p-values are obtained using the beta-hat and standard error for each trait. Each of the SNPs selected in this way is either weakly or not associated with any of the phenotypes (null SNP). Next, select a set of independent null SNPs from the initial set of null SNPs by using a threshold of r^2 &lt; 0.01 (r: the correlation between the genotypes at a pair of SNPs). In the absence of in-sample linkage disequilibrium (LD) information, one can use the reference panel LD information for this screening. Finally, compute the correlation matrix of the effect estimates (beta-hat vector) as the sample correlation matrix of the beta-hat vector across all the selected independent null SNPs. This strategy is more general and applicable to a cohort study or multiple overlapping studies for binary or quantitative traits with arbitrary distributions. It is also useful when the beta-hat vector for multiple non-overlapping studies become correlated due to genetically related individuals across studies. Misspecification of the correlation structure can affect the results produced by CPBayes to some extent. Hence, if genome-wide summary statistics data across traits is available, we recommend this alternative strategy to estimate the correlation matrix of the beta-hat vector.

Getting more details
====================

Please see the vignettes (tutorial) and the manual of CPBayes for more detailed description about the package. If you encounter any problem/issue related to the package, please email to <statgen.arunabha@gmail.com> and <tanushree.haldar@gmail.com>. Also, see our paper for more details: Majumdar A, Haldar T, Bhattacharya S, Witte JS (2018) An efficient Bayesian meta analysis approach for studying cross-phenotype genetic associations. PLoS Genet 14(2): e1007139.
