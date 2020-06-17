library(tidyverse)
library(magrittr)
library(ppcseq)
library(tidybulk)


df =
	read_csv("dev/GSE99374_RawCounts_CD8.csv") %>%
	gather(sample, count, -`Transcript ID`) %>%
	tidyr::extract(col = sample, into = c("sample", "type"), "(.+)_.+_([A-Z]+)") %>%
	tidybulk(sample, `Transcript ID`, count) %>%
	scale_abundance(factor_of_interest = type)


res_dt =
	df %>%
	test_differential_abundance(~ type, method = "edgeR_likelihood_ratio") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )

res_dt %>% saveRDS("dev/GSE99374_CD8_DT.rds", compress = "gzip")

res_dt_deseq2 =
	df %>%
	test_differential_abundance(~ type, method = "deseq2") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 )

res_dt_deseq2 %>% saveRDS("dev/GSE99374_CD8_DT_deseq2.rds", compress = "gzip")

res =
	res_dt %>%
	# PPCSEQ
	identify_outliers(
		formula = ~ type,
		.sample = sample,
		.transcript = `Transcript ID`,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("dev/GSE99374_CD8.rds", compress = "gzip")

res_deseq2 =
	res_dt_deseq2 %>%
	# PPCSEQ
	identify_outliers(
		formula = ~ type,
		.sample = sample,
		.transcript = `Transcript ID`,
		.significance = pvalue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_deseq2 %>% saveRDS("dev/GSE99374_CD8_deseq2.rds", compress = "gzip")


res %>%
	mutate(`Transcript ID` = factor(`Transcript ID`, levels = .$`Transcript ID` %>% unique)) %>%
	ggplot(aes(`Transcript ID`, `tot deleterious outliers`)) + geom_point() + my_theme

res %>% pull(plot) %>% cowplot::plot_grid(plotlist = ., align = "h", ncol = 3, axis="b" )
