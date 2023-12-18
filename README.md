
<!-- README.md is generated from README.Rmd. Please edit that file -->

# biomuts

<!-- badges: start -->
<!-- badges: end -->

## installation

You can install the development version of `biomuts` like so:

``` r
devtools::install_github("william-swl/biomuts")
```

## utils

- correct the gaps location to fit codon frames

``` r
DNA_gaps_corr("AT---G")
#> [1] "ATG---"
```

- translate DNA to AA (gaps allowed)

``` r
dna2aa("ATG---AAA")
#> [1] "M-K"
```
