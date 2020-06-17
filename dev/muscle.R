library(tidyverse)
library(purrr)
library(magrittr)
library(ppcseq)
library(tidybulk)


df =
	read_tsv("dev/GSE137631_raw_counts.txt")  %>%
	gather(sample, count, -X1) %>%
	rename(ens = X1) %>%
	filter(!grepl("control|PRE", sample)) %>%
	mutate(training = grepl("ET", sample)) %>%
	tidybulk(sample, ens, count) %>%
	scale_abundance(factor_of_interest = training) %>%
	# mutate(is_sp1 = grepl("SP1", sample)) %>%
	# adjust_abundance(~type + is_sp1) %>%

	{
		((.) %>%
		 	# filter(count_scaled_adjusted %>% is.na %>% `!`) %>%
			reduce_dimensions(method="MDS", action = "get") %>%
			ggplot(aes(Dim1, Dim2, color=training, label=sample)) + geom_point() + geom_text()) %>% print
		(.)
	}

res_dt =
	df %>%
	test_differential_abundance(~ training, method = "edgeR_likelihood_ratio") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )


res_dt %>% saveRDS("dev/GSE137631_muscle_DT.rds", compress = "gzip")

res_dt_deseq2 =
	df %>%
	test_differential_abundance(~ training, method = "deseq2") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 )


res_dt_deseq2 %>% saveRDS("dev/GSE137631_muscle_DT_deseq2.rds", compress = "gzip")



res =
	res_dt %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ training ,
		.sample = sample,
		.transcript = ens,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("dev/GSE137631_muscle.rds", compress = "gzip")

res_deseq2 =
	res_dt_deseq2 %>%
	filter(padj %>% is.na %>% `!`) %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ training ,
		.sample = sample,
		.transcript = ens,
		.significance = pvalue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_deseq2 %>% saveRDS("dev/GSE137631_muscle_deseq2.rds", compress = "gzip")


res %>%
	mutate(`Transcript ID` = factor(`Transcript ID`, levels = .$`Transcript ID` %>% unique)) %>%
	ggplot(aes(`Transcript ID`, `tot deleterious outliers`)) + geom_point() + my_theme

res %>% pull(plot) %>% cowplot::plot_grid(plotlist = ., align = "h", ncol = 3, axis="b" )
