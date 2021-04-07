#' The 'ppcseq' package.
#'
#' @description Relative transcript abundance has proven to be a valuable tool for understanding the function of genes in biological systems. For the differential analysis of transcript abundance using RNA sequencing data, the negative binomial model is by far the most frequently adopted. However, common methods that are based on a negative binomial model are not robust to extreme outliers, which we found to be abundant in public datasets. So far, no rigorous and probabilistic methods for detection of outliers have been developed for RNA sequencing data, leaving the identification mostly to visual inspection. Recent advances in Bayesian computation allow large-scale comparison of observed data against its theoretical distribution given in a statistical model. Here we propose ppcseq, a key quality-control tool for identifying transcripts that include outlier data points in differential expression analysis, which do not follow a negative binomial distribution. Applying ppcseq to analyse several publicly available datasets using popular tools, we show that from 3 to 10 percent of differentially abundant transcripts across algorithms and datasets had statistics inflated by the presence of outliers.
#'
#' @docType package
#' @name ppcseq-package
#' @aliases ppcseq
#' @useDynLib ppcseq, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @import rstantools
#' @import lifecycle
#' @importFrom graphics par
#'
#' @usage data(counts)
#'
#' @return NULL
#'
#' @references
#' Stan Development Team (2019). RStan: the R interface to Stan. R package version 2.19.2. https://mc-stan.org
#' angiola S, Thomas E, Modrak M, Vehtari A, Papenfuss A (2021). “Probabilistic outlier identification for RNA sequencing generalized linear models.” NAR Genomics and Bioinformatics_, *3*(1), lqab005. <URL: https://doi.org/10.1093/nargab/lqab005>.
NULL
