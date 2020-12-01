# CPBayes 1.1.0
Fifth version of CPBayes R-package

## New features

None.

## Bug fixes and minor improvements

The default choices of the burn in period and the MCMC sample size are modified. Computational runtime is reduced significantly.

# CPBayes 1.0.0
Fourth version of CPBayes R-package

## New features

Two new functions are introduced to analytically compute the local false discovery rate (locFDR) & Bayes factor (BF) that quantifies the evidence of aggregate-level pleiotropic association for uncorrelated and correlated summary statistics.

## Bug fixes and minor improvements

Instead of locFDR and optimal subset of non-null traits, the cpbayes_uncor() and cpbayes_cor() functions now print a list of important traits underlying an overall pleiotropic association.

# CPBayes 0.3.0
Third version of CPBayes R-package

## New features

An empirical Bayes approach to choose the hyper-parameters in the prior distribution of proportion of non-null traits is included in the Bayesian statistical framework underlying the meta-analysis method.

## Bug fixes and minor improvements

An additional argument is included in the forest_cpbayes() function that allows to plot a subset of selected traits having the trait-specific posterior probability of association above a specified threshold.  

Instead of the Bayes factor, the cpbayes_uncor() and cpbayes_cor() functions now print the local false discovery rate (locFDR) as the primary measure of overall pleiotropic association.

# CPBayes 0.2.0
Second version of CPBayes R-package

## New features

New forest_cpbayes() function to make a forest plot that provides a graphical presentation of the pleiotropy results obtained by CPBayes.

## Bug fixes and minor improvements

In the cpbayes_uncor() and cpbayes_cor() functions, we have renamed the argument `UpdateDE' as `UpdateSlabVar' (which is an abbreviation of Update Slab Variance).

New arguments option included in the cpbayes_uncor() and cpbayes_cor() functions. An user can now specify the minimum and maximum value of the slab variance parameter through the 'MinSlabVar' and 'MaxSlabVar' arguments.

In the previous version, the row names and column names of the correlation matrix of the beta hat vector in the cpbayes_cor() function had to be the same to pass the symmetricity checking. Now it is rectified and the row names and column names of the correlation matrix can be different.



# CPBayes 0.1.0

First version of CPBayes R-package.



