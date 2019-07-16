Posterior predictive check for bulk RNA sequencing data
================

The input data set is a tidy representation of a differential gene
transcript abundance analysis

``` r
library(tidyverse)
```

    ## Registered S3 methods overwritten by 'ggplot2':
    ##   method         from 
    ##   [.quosures     rlang
    ##   c.quosures     rlang
    ##   print.quosures rlang

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.1     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(ppcSeq)
```

``` r
counts # Available at ppcSeq::counts
```

You can convert a list of BAM/SAM files into a tidy data frame of
annotated counts

``` r
counts.ppc = 
    ppcSeq::counts %>%
    mutate(is_significant = FDR < 0.05) %>%
    ppc_seq(
        formula = ~ Label,
        significance_column = PValue,
        do_check_column = is_significant,
        value_column = value
    )
```

    ## [1] "2019-07-16 17:28:01 AEST"
    ## [1] "2019-07-16 17:29:47 AEST"
    ## [1] "2019-07-16 17:29:58 AEST"
    ## [1] "2019-07-16 17:33:33 AEST"

The new posterior predictive check has been added to the original data
frame

``` r
counts.ppc 
```

    ## # A tibble: 70 x 4
    ##    symbol   `sample wise data` plot   `tot deleterious outliers`
    ##    <chr>    <list>             <list>                      <int>
    ##  1 SLC16A12 <tibble [21 × 11]> <gg>                            0
    ##  2 CYP1A1   <tibble [21 × 11]> <gg>                            1
    ##  3 ART3     <tibble [21 × 11]> <gg>                            1
    ##  4 DIO2     <tibble [21 × 11]> <gg>                            0
    ##  5 OR51E2   <tibble [21 × 11]> <gg>                            0
    ##  6 MUC16    <tibble [21 × 11]> <gg>                            0
    ##  7 CCNA1    <tibble [21 × 11]> <gg>                            0
    ##  8 LYZ      <tibble [21 × 11]> <gg>                            1
    ##  9 PPM1H    <tibble [21 × 11]> <gg>                            0
    ## 10 SUSD5    <tibble [21 × 11]> <gg>                            0
    ## # … with 60 more rows

The new data frame contains plots for each gene

We can visualise the top five differentially transcribed genes

``` r
counts.ppc %>% 
    slice(1:2) %>% 
    pull(plot) %>% 
    cowplot::plot_grid(plotlist = ., align = "v", ncol = 1, axis="b", rel_widths = 1 )
```

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
