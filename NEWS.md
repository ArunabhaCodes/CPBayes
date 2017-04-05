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



