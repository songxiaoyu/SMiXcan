
# `SMiXcan: XXX`

A full description of the method can be found at: XXX

Reference: XXX.

## Installation

With R, users can install the QuadST package directly from GitHub with
devtools:

``` r
if (!require("devtools", quietly = TRUE)) install.packages("devtools")
library(devtools)
devtools::install_github("songxiaoyu/SMiXcan")
```

## Example of using QuadST analysis pipeline

``` r
library(SMiXcan)
sessionInfo() 
```

    ## R version 4.4.1 (2024-06-14)
    ## Platform: aarch64-apple-darwin20
    ## Running under: macOS Sonoma 14.4
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRblas.0.dylib 
    ## LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/Los_Angeles
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] SMiXcan_0.1.0  devtools_2.4.5 usethis_3.1.0 
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyr_1.3.1       generics_0.1.4    shape_1.4.6.1     lattice_0.22-7   
    ##  [5] digest_0.6.37     magrittr_2.0.3    grid_4.4.1        evaluate_1.0.4   
    ##  [9] ACAT_0.91         iterators_1.0.14  pkgload_1.4.0     fastmap_1.2.0    
    ## [13] Matrix_1.7-3      glmnet_4.1-10     foreach_1.5.2     doParallel_1.0.17
    ## [17] pkgbuild_1.4.8    sessioninfo_1.2.3 backports_1.5.0   survival_3.8-3   
    ## [21] urlchecker_1.0.1  promises_1.3.3    purrr_1.1.0       codetools_0.2-20 
    ## [25] cli_3.6.5         shiny_1.10.0      rlang_1.1.6       splines_4.4.1    
    ## [29] ellipsis_0.3.2    remotes_2.5.0     cachem_1.1.0      yaml_2.3.10      
    ## [33] tools_4.4.1       parallel_4.4.1    memoise_2.0.1     dplyr_1.1.4      
    ## [37] httpuv_1.6.16     MiXcan_0.1.0      broom_1.0.8       curl_6.4.0       
    ## [41] vctrs_0.6.5       R6_2.6.1          mime_0.13         lifecycle_1.0.4  
    ## [45] fs_1.6.6          htmlwidgets_1.6.4 miniUI_0.1.2      pkgconfig_2.0.3  
    ## [49] pillar_1.11.0     later_1.4.2       glue_1.8.0        profvis_0.4.0    
    ## [53] Rcpp_1.1.0        xfun_0.52         tibble_3.3.0      tidyselect_1.2.1 
    ## [57] rstudioapi_0.17.1 knitr_1.50        xtable_1.8-4      htmltools_0.5.8.1
    ## [61] rmarkdown_2.29    compiler_4.4.1

### Load a public data set

``` r
SMiXcan::safe_ACAT(c(0.5, 0.5))
```
