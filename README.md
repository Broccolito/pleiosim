
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PleioSim <a href='https://github.com/Broccolito/PleioSim'><img src='man/figures/logo.png' align="right" height="140"/></a>

<!-- badges: start -->

<br> <!-- badges: end -->

## Installation

``` r
# Install devtools if it is not yet installed
if(!require("devtools")){
  install.packages("devtools")
}

# Install senlinplot
if(!require("pleiosim")){
  devtools::install_github("Broccolito/pleiosim")
  library("pleiosim")
}
```

## Usage

Load the PleioSim package, set seed for the simulation

``` r
library(pleiosim)
set.seed(492357816)
```
