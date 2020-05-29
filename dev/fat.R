library(tidyverse)
library(magrittr)
library(ppcseq)

res =
	identify_outliers(
		dplyr::mutate(ppcseq::counts,  is_significant = FDR < 0.05 ),
		formula = ~ Label,
		.sample = sample,
		.transcript = symbol,
		.abundance = value,
		.significance = PValue,
		.do_check  = is_significant,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("dev/Mangiola_et_adipo.rds", compress = "gzip")

