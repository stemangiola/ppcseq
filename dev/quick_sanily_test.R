library(tidyverse)
library(magrittr)
library(ppcSeq)
library(furrr)
plan(multicore)

FDR_threshold = 0.01

res_1 =
	ppcSeq::counts %>%
	mutate(is_significant = FDR < FDR_threshold) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column  = is_significant,
		value_column = value,
		percent_false_positive_genes = "5%", tol_rel_obj = 0.01
	)
