
<!-- README.md is generated from README.Rmd. Please edit that file -->

<img src="https://i.imgur.com/kO2FL3r.png" alt="HodgesTools_HexLogo" width="200"/>

# **HodgesTools**

## **Overview**

HodgesTools is an R package built by hodges lab members for current and
future hodges lab members as well as anyone else who might find the
tools useful. This package contains functions that the lab uses everyday
to analyze various genomic datasets. Critically, this package only
contains general use functions; functions specific to a given technique
are reserved for a separate package. As the lab grows, we expect to
continue adding functions to the package to build on previous lab
members code.

## **Examples**

Below is an example for the read\_bed function.

#### *read\_bed*

The *read\_bed* function reads a tab-delimited BED formatted file into
R. It auto-detects BED3, BED6, BED3+, and BED6+ formats. Below is an
example for a BED6 file.

``` r
# load package
devtools::load_all() #this is development code, use library(HodgesTools)
#> i Loading HodgesTools

#load data from extdata directory
BED6 <- system.file(package = "HodgesTools", "extdata", "test_BED6.bed")

#read 6 column bedfile into R
bed6 <- read_bed(BED6)
#> User-supplied file read as 6-column BED.
head(bed6)
#> # A tibble: 4 x 6
#>   chr   start   end name  score strand
#>   <chr> <dbl> <dbl> <chr> <dbl> <chr> 
#> 1 chr1    123   456 NAME1   123 .     
#> 2 chr1    789  1908 NAME2   123 -     
#> 3 chr1   1789  2908 NAME3   123 +     
#> 4 chr1   3789  3908 NAME4   123 -
```

## **Installation**

This package can be installed via CRAN:

``` r
install.packages('HodgesTools')
```

## **Tools**

The following functions are provided in this release:

  - read\_bed
  - plot\_HOMERTFs
  - cpg\_analysis
  - createManhattandQQ
  - append\_section\_to\_ini

*Note: Usage and arguments for each function can be viewed in the
manual.*
