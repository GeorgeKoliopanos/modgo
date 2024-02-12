
# modgo

<!-- badges: start -->
<!-- badges: end -->

modgo is an R package Mock Data Generation (modgo) that may be used for
simulating data from existing study data for continuous, ordinal categorical,
and dichotomous variables[reference].

## Installation

To install modgo from GitHub, run: 

``` r
library("devtools")
install_github("https://github.com/GeorgeKoliopanos/modgo")
```

## Usage

To see how to use modgo see ?modgo() in R . Below, we present a simple example
on how to run.

``` r
library(modgo)
data("Cleveland",package="modgo")
test_modgo <- modgo(data = Cleveland,
     bin_variables = c("CAD","HighFastBloodSugar","Sex","ExInducedAngina"),
     categ_variables =c("Chestpaintype"))
```

