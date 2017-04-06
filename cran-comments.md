## Resubmission
This is a resubmission. In this version:

* I have reduced the test run time to 27-28 sec (less than 3 minutes) [checked with "R CMD check --as-cran"]

## Test environments
* local OS X El Capitan, R 3.2.3, R 3.3.3
* x86_64-pc-linux-gnu (on travis-ci), R 3.3.3
* win-builder (devel(3.4.0 alpha) and release (3.3.3))

## R CMD check results

We found 0 error, 0 warning and 0 note while checking on Linux and local mac.

The result of build_win() for both devel(3.4.0 alpha) and release (3.3.3) is 
0 errors | 0 warnings | 1 note

Note is regarding the possible misspelling in the 'DESCRIPTION' file. Words listed here are: 'CPBayes', 'phenotypes', 'pleiotropic'. We confirm that these words are correctly spelled.

## Reverse dependencies

There is no reverse dependencies, checked with devtools::revdep_check()

