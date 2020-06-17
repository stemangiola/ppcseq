library(tidyverse)
library(magrittr)
library(ppcseq)

res_dt = dplyr::mutate(ppcseq::counts,  is_significant = FDR < 0.05 ) %>% select(symbol, everything())

res_dt %>% saveRDS("dev/Mangiola_et_adipo_DT.rds", compress = "gzip")

res_dt_deseq2 =
	ppcseq::counts %>%
	distinct(sample, symbol, value, W, Label) %>%
	test_differential_abundance(~ Label + W , sample, symbol, value, method="deseq2") %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 ) %>%
select(symbol, everything())

res_dt_deseq2 %>% saveRDS("dev/Mangiola_et_adipo_DT_deseq2.rds", compress = "gzip")


res =
	res_dt %>%

	identify_outliers(
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

res_deseq2 =
	res_dt_deseq2 %>%
	filter(!lowly_abundant) %>%

	identify_outliers(
		formula = ~ Label,
		.sample = sample,
		.transcript = symbol,
		.abundance = value,
		.significance = pvalue,
		.do_check  = is_significant,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_deseq2 %>% saveRDS("dev/Mangiola_et_adipo_deseq2.rds", compress = "gzip")

