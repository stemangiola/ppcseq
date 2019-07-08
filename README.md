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
counts %>% 
    mutate(is_significant = FDR < 0.05) %>% 
    ppc_seq(
        formula = ~ Label + W, 
        significance_column = "PValue", 
        do_check_column = "is_significant", 
        value_column = "value"
    )
```

    ## Joining, by = "symbol"
    ## Joining, by = "symbol"

    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"
    ## Joining, by = "G"

    ## Joining, by = c("G", "idx_MPI")

    ## [1] "2019-07-08 17:56:32 AEST"
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: EXPERIMENTAL ALGORITHM:
    ## Chain 1:   This procedure has not been thoroughly tested and may be unstable
    ## Chain 1:   or buggy. The interface is subject to change.
    ## Chain 1: ------------------------------------------------------------
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: 
    ## Chain 1: Gradient evaluation took 0.02 seconds
    ## Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 200 seconds.
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
    ## Chain 1:    100     -7229661.426             1.000            1.000
    ## Chain 1:    200     -2794519.673             1.294            1.587
    ## Chain 1:    300     -1361491.168             1.213            1.053
    ## Chain 1:    400      -737507.961             1.121            1.053
    ## Chain 1:    500      -465792.971             1.014            1.000
    ## Chain 1:    600      -341679.371             0.905            1.000
    ## Chain 1:    700      -283403.496             0.805            0.846
    ## Chain 1:    800      -252923.588             0.720            0.846
    ## Chain 1:    900      -236177.601             0.648            0.583
    ## Chain 1:   1000      -226074.404             0.587            0.583
    ## Chain 1:   1100      -218925.512             0.537            0.363   MAY BE DIVERGING... INSPECT ELBO
    ## Chain 1:   1200      -213951.897             0.494            0.363
    ## Chain 1:   1300      -209094.996             0.458            0.206
    ## Chain 1:   1400      -205086.154             0.427            0.206
    ## Chain 1:   1500      -201482.907             0.399            0.121
    ## Chain 1:   1600      -198193.286             0.375            0.121
    ## Chain 1:   1700      -194527.086             0.354            0.071
    ## Chain 1:   1800      -191680.671             0.336            0.071
    ## Chain 1:   1900      -188729.688             0.319            0.045
    ## Chain 1:   2000      -186090.700             0.304            0.045
    ## Chain 1:   2100      -183618.019             0.290            0.033
    ## Chain 1:   2200      -181198.496             0.277            0.033
    ## Chain 1:   2300      -179087.329             0.266            0.023
    ## Chain 1:   2400      -176955.497             0.255            0.023
    ## Chain 1:   2500      -175033.088             0.245            0.023
    ## Chain 1:   2600      -173229.130             0.236            0.023
    ## Chain 1:   2700      -171584.452             0.228            0.020
    ## Chain 1:   2800      -169989.845             0.220            0.020
    ## Chain 1:   2900      -168545.320             0.213            0.019
    ## Chain 1:   3000      -167142.561             0.206            0.019
    ## Chain 1:   3100      -165859.806             0.200            0.018
    ## Chain 1:   3200      -164631.085             0.194            0.018
    ## Chain 1:   3300      -163455.538             0.188            0.017
    ## Chain 1:   3400      -162431.336             0.183            0.017
    ## Chain 1:   3500      -161415.008             0.178            0.016
    ## Chain 1:   3600      -160474.404             0.173            0.016
    ## Chain 1:   3700      -159595.811             0.168            0.015
    ## Chain 1:   3800      -158715.755             0.164            0.015
    ## Chain 1:   3900      -157932.142             0.160            0.014
    ## Chain 1:   4000      -157185.327             0.156            0.014
    ## Chain 1:   4100      -156466.719             0.152            0.013
    ## Chain 1:   4200      -155790.122             0.149            0.013
    ## Chain 1:   4300      -155205.487             0.145            0.013
    ## Chain 1:   4400      -154594.641             0.142            0.013
    ## Chain 1:   4500      -154053.365             0.139            0.012
    ## Chain 1:   4600      -153507.057             0.136            0.012
    ## Chain 1:   4700      -153017.845             0.133            0.012
    ## Chain 1:   4800      -152556.236             0.131            0.012
    ## Chain 1:   4900      -152125.977             0.128            0.011
    ## Chain 1:   5000      -151708.852             0.126            0.011
    ## Chain 1:   5100      -151310.496             0.106            0.010
    ## Chain 1:   5200      -150947.884             0.074            0.010
    ## Chain 1:   5300      -150599.250             0.053            0.009
    ## Chain 1:   5400      -150279.302             0.036            0.009
    ## Chain 1:   5500      -149960.494             0.024            0.008
    ## Chain 1:   5600      -149664.605             0.017            0.008
    ## Chain 1:   5700      -149413.137             0.013            0.007
    ## Chain 1:   5800      -149139.790             0.011            0.007
    ## Chain 1:   5900      -148897.458             0.009            0.006
    ## Chain 1:   6000      -148653.688             0.008            0.006
    ## Chain 1:   6100      -148464.908             0.008            0.006
    ## Chain 1:   6200      -148253.021             0.007            0.006
    ## Chain 1:   6300      -148084.245             0.007            0.006
    ## Chain 1:   6400      -147881.522             0.007            0.005   MEDIAN ELBO CONVERGED
    ## Chain 1: 
    ## Chain 1: Drawing a sample of size 1000 from the approximate posterior... 
    ## Chain 1: COMPLETED.
    ## [1] "2019-07-08 17:58:12 AEST"

    ## Warning: Expected 3 pieces. Additional pieces discarded in 22470 rows [1,
    ## 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Warning: Expected 2 pieces. Additional pieces discarded in 21 rows [1, 2,
    ## 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].

    ## Joining, by = "S"

    ## Joining, by = c("S", "G")

    ## Joining, by = "symbol"

    ## # A tibble: 394,821 x 12
    ##    sample symbol logCPM    LR   PValue     FDR value       W Label
    ##    <chr>  <chr>   <dbl> <dbl>    <dbl>   <dbl> <int>   <dbl> <chr>
    ##  1 10922… SLC16…   1.39  41.1 1.46e-10 2.74e-6   160 -0.129  High 
    ##  2 10935… SLC16…   1.39  41.1 1.46e-10 2.74e-6   150 -0.127  High 
    ##  3 10973… SLC16…   1.39  41.1 1.46e-10 2.74e-6   146 -0.426  High 
    ##  4 10976… SLC16…   1.39  41.1 1.46e-10 2.74e-6   347 -0.0164 High 
    ##  5 10985… SLC16…   1.39  41.1 1.46e-10 2.74e-6   175 -0.135  High 
    ##  6 11026… SLC16…   1.39  41.1 1.46e-10 2.74e-6   244  0.125  High 
    ##  7 11045… SLC16…   1.39  41.1 1.46e-10 2.74e-6   399 -0.0892 High 
    ##  8 11082… SLC16…   1.39  41.1 1.46e-10 2.74e-6   100  0.261  Neoa…
    ##  9 11086… SLC16…   1.39  41.1 1.46e-10 2.74e-6    37 -0.132  Neoa…
    ## 10 11103… SLC16…   1.39  41.1 1.46e-10 2.74e-6    73  0.146  Neoa…
    ## # … with 394,811 more rows, and 3 more variables: is_significant <lgl>,
    ## #   plot <list>, `ppc samples failed` <int>
