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

    ## ── Attaching packages ─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✔ ggplot2 3.1.1     ✔ purrr   0.3.2
    ## ✔ tibble  2.1.1     ✔ dplyr   0.8.1
    ## ✔ tidyr   0.8.3     ✔ stringr 1.4.0
    ## ✔ readr   1.3.1     ✔ forcats 0.4.0

    ## ── Conflicts ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
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
    counts %>% 
        mutate(is_significant = FDR < 0.05) %>% 
        ppc_seq(
            formula = ~ Label + W, 
            significance_column = "PValue", 
            do_check_column = "is_significant", 
            value_column = "value"
        )
```

    ## [1] "2019-07-09 10:16:21 AEST"
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.03 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 300 seconds.
    ## Chain 1: Adjust your expectations accordingly!
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Begin eta adaptation.
    ## Chain 1: Iteration:   1 / 250 [  0%]  (Adaptation)
    ## Chain 1: Iteration:  50 / 250 [ 20%]  (Adaptation)
    ## Chain 1: Iteration: 100 / 250 [ 40%]  (Adaptation)
    ## Chain 1: Iteration: 150 / 250 [ 60%]  (Adaptation)
    ## Chain 1: Iteration: 200 / 250 [ 80%]  (Adaptation)
    ## Chain 1: Iteration: 250 / 250 [100%]  (Adaptation)
    ## Chain 1: Success! Found best value [eta = 0.1].
    ## Chain 1: 
    ## Chain 1: Begin stochastic gradient ascent.
    ## Chain 1:   iter             ELBO   delta_ELBO_mean   delta_ELBO_med   notes 
    ## Chain 1:    100     -8267961.123             1.000            1.000
    ## Chain 1:    200     -3023712.844             1.367            1.734
    ## Chain 1:    300     -1500125.240             1.250            1.016
    ## Chain 1:    400      -803505.836             1.154            1.016
    ## Chain 1:    500      -499179.933             1.045            1.000
    ## Chain 1:    600      -356529.160             0.938            1.000
    ## Chain 1:    700      -289907.868             0.837            0.867
    ## Chain 1:    800      -254354.664             0.750            0.867
    ## Chain 1:    900      -235398.090             0.675            0.610
    ## Chain 1:   1000      -224403.662             0.613            0.610
    ## Chain 1:   1100      -216862.102             0.560            0.400   MAY BE DIVERGING... INSPECT ELBO
    ## Chain 1:   1200      -211215.770             0.516            0.400   MAY BE DIVERGING... INSPECT ELBO
    ## Chain 1:   1300      -206679.724             0.478            0.230
    ## Chain 1:   1400      -202875.636             0.445            0.230
    ## Chain 1:   1500      -199065.830             0.416            0.140
    ## Chain 1:   1600      -195912.544             0.391            0.140
    ## Chain 1:   1700      -192762.265             0.369            0.081
    ## Chain 1:   1800      -189974.874             0.350            0.081
    ## Chain 1:   1900      -187254.976             0.332            0.049
    ## Chain 1:   2000      -184644.520             0.316            0.049
    ## Chain 1:   2100      -182171.919             0.302            0.035
    ## Chain 1:   2200      -179966.137             0.289            0.035
    ## Chain 1:   2300      -177892.863             0.277            0.027
    ## Chain 1:   2400      -175870.507             0.265            0.027
    ## Chain 1:   2500      -174062.555             0.255            0.022
    ## Chain 1:   2600      -172287.641             0.246            0.022
    ## Chain 1:   2700      -170692.779             0.237            0.019
    ## Chain 1:   2800      -169203.418             0.229            0.019
    ## Chain 1:   2900      -167779.569             0.221            0.019
    ## Chain 1:   3000      -166419.294             0.214            0.019
    ## Chain 1:   3100      -165152.352             0.208            0.016
    ## Chain 1:   3200      -163954.705             0.201            0.016
    ## Chain 1:   3300      -162890.884             0.195            0.016
    ## Chain 1:   3400      -161819.849             0.190            0.016
    ## Chain 1:   3500      -160842.943             0.185            0.015
    ## Chain 1:   3600      -159915.668             0.180            0.015
    ## Chain 1:   3700      -159075.808             0.175            0.015
    ## Chain 1:   3800      -158270.676             0.170            0.015
    ## Chain 1:   3900      -157498.352             0.166            0.014
    ## Chain 1:   4000      -156775.570             0.162            0.014
    ## Chain 1:   4100      -156111.651             0.158            0.014
    ## Chain 1:   4200      -155453.624             0.155            0.014
    ## Chain 1:   4300      -154880.171             0.151            0.012
    ## Chain 1:   4400      -154298.739             0.148            0.012
    ## Chain 1:   4500      -153761.467             0.145            0.012
    ## Chain 1:   4600      -153260.332             0.142            0.012
    ## Chain 1:   4700      -152789.857             0.139            0.011
    ## Chain 1:   4800      -152333.971             0.136            0.011
    ## Chain 1:   4900      -151934.542             0.133            0.010
    ## Chain 1:   5000      -151530.257             0.130            0.010
    ## Chain 1:   5100      -151132.324             0.110            0.010
    ## Chain 1:   5200      -150794.302             0.076            0.009
    ## Chain 1:   5300      -150468.407             0.056            0.009
    ## Chain 1:   5400      -150141.596             0.038            0.008
    ## Chain 1:   5500      -149856.513             0.026            0.008
    ## Chain 1:   5600      -149596.979             0.018            0.008
    ## Chain 1:   5700      -149314.856             0.014            0.007
    ## Chain 1:   5800      -149069.269             0.011            0.007
    ## Chain 1:   5900      -148850.000             0.009            0.007
    ## Chain 1:   6000      -148619.095             0.008            0.006
    ## Chain 1:   6100      -148412.309             0.008            0.006
    ## Chain 1:   6200      -148219.508             0.007            0.005
    ## Chain 1:   6300      -148049.442             0.007            0.005
    ## Chain 1:   6400      -147860.524             0.006            0.005   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.
    ## [1] "2019-07-09 10:18:07 AEST"

The new posterior predictive check has been added to the original data
frame

``` r
counts.ppc %>% 
    distinct(symbol, FDR, `ppc samples failed`, plot)
```

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `plot`

    ## # A tibble: 18,801 x 4
    ##    symbol          FDR `ppc samples failed` plot  
    ##    <chr>         <dbl>                <int> <list>
    ##  1 SLC16A12 0.00000274                    1 <gg>  
    ##  2 CYP1A1   0.000175                      0 <gg>  
    ##  3 ART3     0.000217                      1 <gg>  
    ##  4 DIO2     0.000483                      0 <gg>  
    ##  5 OR51E2   0.000785                      0 <gg>  
    ##  6 MUC16    0.00657                       0 <gg>  
    ##  7 CCNA1    0.00734                       1 <gg>  
    ##  8 LYZ      0.00734                       1 <gg>  
    ##  9 PPM1H    0.00734                       0 <gg>  
    ## 10 SUSD5    0.00734                       0 <gg>  
    ## # … with 18,791 more rows

The new data frame contains plots for each gene

``` r
counts.ppc %>% 
    distinct(symbol, FDR, `ppc samples failed`, plot) %>%
    slice(1) %>%
    pull(plot) %>%
    `[[` (1)
```

    ## Warning: distinct() does not fully support columns of type `list`.
    ## List elements are compared by reference, see ?distinct for details.
    ## This affects the following columns:
    ## - `plot`

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->
