library(tidyverse)
library(magrittr)
library(ppcseq)
library(tidybulk)
library(DESeq2)

res_dt = dplyr::mutate(ppcseq::counts,  is_significant = FDR < 0.05 ) %>% select(symbol, everything())

res_dt %>% saveRDS("dev/Mangiola_et_adipo_DT.rds", compress = "gzip")

res_dt_robust =
	ppcseq::counts %>%
	distinct(sample, symbol, value, W, Label) %>%
	tidybulk(sample, symbol, value) %>%
	identify_abundant(factor_of_interest = "Label") %>%
	test_differential_abundance(~ Label + W , method="edgeR_robust_likelihood_ratio") %>%
	mutate(is_significant = !is.na(FDR) & FDR < 0.05 ) %>%
	select(symbol, everything()) %>%
	filter(.abundant)

res_dt_deseq2 =
	ppcseq::counts %>%
	distinct(sample, symbol, value, W, Label) %>%
	tidybulk(sample, symbol, value) %>%
	identify_abundant(factor_of_interest = "Label") %>%
	test_differential_abundance(~ Label + W , method="deseq2") %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 ) %>%
	select(symbol, everything())

res_dt_deseq2 %>% saveRDS("dev/Mangiola_et_adipo_DT_deseq2.rds", compress = "gzip")



res =
	res_dt %>%

	identify_outliers(
		formula = ~ Label + W ,
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

res_robust =
	res_dt_robust %>%

	identify_outliers(
		formula = ~ Label + W ,
		.sample = sample,
		.transcript = symbol,
		.abundance = value,
		.significance = PValue,
		.do_check  = is_significant,
		percent_false_positive_genes = 1,
		approximate_posterior_analysis = T,
		cores=10
	)

res_robust %>% select(-plot, -`sample wise data`) %>% saveRDS("dev/Mangiola_et_adipo_robust.rds", compress = "gzip")

# Best rank outlier
res_dt_robust %>% filter(is_significant) %>% arrange(FDR) %>% pivot_transcript %>% rowid_to_column() %>% inner_join(res_robust %>% filter( `tot deleterious outliers`> 0) %>% distinct(symbol))


res_deseq2 =
	res_dt_deseq2 %>%
	filter(!lowly_abundant) %>%

	identify_outliers(
		formula = ~ Label + W ,
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

# How many outliers with Cook method
res_dt_deseq2 %>%
	attr("internals") %$%
	DESeq2 %>%
	assays() %>%
	.[["cooks"]] %>%
	apply(2, function(x) x > (4.874046)) %>%
	rowSums() %>%
	enframe(name = "symbol", value = "outliers") %>%
	filter(outliers > 0) %>%
	inner_join(
		res_dt_deseq2 %>% filter(is_significant) %>% pivot_transcript()
	)
