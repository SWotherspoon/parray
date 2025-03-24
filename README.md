
# parray

<!-- badges: start -->
<!-- badges: end -->

The parray package provides facilities for manipulating arrays of conditional probabilities.

## Installation

The current version of parray can be installed from GitHub using the remotes package. 
```r
# install.packages("remotes")
remotes::install_github("SWotherspoon/parray")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(parray)
## Probability of cylinder and carburetor combinations conditional on numbers of gears
ptabs(~carb+cyl|gear,data=mtcars)
## Joint probability of numbers of cyclinders, carburetors and gears
ptabs(~carb+cyl+gear,data=mtcars)
## Joint probability as product of conditional and marginal probabilities
product(ptabs(~cyl+carb|gear,data=mtcars),ptabs(~gear,data=mtcars))
```

