library(tidyverse)
library(purrr)
library(magrittr)
library(ppcseq)
library(tidybulk)

options('future.globals.maxSize' = 10014*1024^2)

df =
	dir("dev/GSE141027_RAW/", full.names = T) %>%
	map(~ .x %>% read_tsv(col_names = F) %>% setNames(c("ens_iso", .x))) %>%
	purrr::reduce(left_join, by = "ens_iso") %>%
	gather(file, count, -ens_iso) %>%
	tidyr::extract(
		col = file,
		into = c( "sample", "type", "replicate"),
		"dev/GSE141027_RAW//([A-Z0-9]+)_([A-Za-z]+)_([0-9]+)_.+",
		remove=F,
		convert = T
	) %>%
	tidybulk(sample, ens_iso, count) %>%
	identify_abundant(factor_of_interest = type) %>%
	scale_abundance() %>%

	{
		((.) %>%
			reduce_dimensions(method="MDS", action = "get") %>%
			ggplot(aes(Dim1, Dim2, color=type)) + geom_point()) %>% print
		(.)
	}

res_dt =
	df %>%
	test_differential_abundance(~ type, method = "edgeR_likelihood_ratio") %>%
	filter(.abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )

res_dt %>% saveRDS("dev/GSE141027_lipoma_DT.rds", compress = "gzip")

res_dt_robust =
	df %>%
	test_differential_abundance(~ type, method = "edgeR_robust_likelihood_ratio") %>%
	filter(.abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )

res_dt_robust %>% saveRDS("dev/GSE141027_lipoma_DT_robust.rds", compress = "gzip")

res_dt_deseq2 =
	df %>%
	test_differential_abundance(~ type, method = "deseq2") %>%
	filter(.abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 )

res_dt_deseq2 %>% saveRDS("dev/GSE141027_lipoma_DT_deseq2.rds", compress = "gzip")


res =
	res_dt %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type,
		.sample = sample,
		.transcript = ens_iso,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("dev/GSE141027_lipoma.rds", compress = "gzip")

res_robust =
	res_dt_robust %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type,
		.sample = sample,
		.transcript = ens_iso,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		approximate_posterior_analysis = T,
		cores=30
	)

res_robust %>% select(-plot, -`sample wise data`) %>% saveRDS("dev/GSE141027_lipoma_robust.rds", compress = "gzip")

# Best rank outlier
res_dt_robust %>% filter(is_significant) %>% arrange(FDR) %>% pivot_transcript %>% rowid_to_column() %>% inner_join(res_robust2 %>% filter( `tot deleterious outliers`> 0) %>% distinct(ens_iso))

res_deseq2 =
	res_dt_deseq2 %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type,
		.sample = sample,
		.transcript = ens_iso,
		.significance = pvalue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_deseq2 %>% saveRDS("dev/GSE141027_lipoma_deseq2.rds", compress = "gzip")

# How many outliers with Cook method
res_dt_deseq2 %>%
	attr("internals") %$%
	DESeq2 %>%
	assays() %>%
	.[["cooks"]] %>%
	apply(2, function(x) x > (4.874046)) %>%
	rowSums() %>%
	enframe(name = "ens_iso", value = "outliers") %>%
	filter(outliers > 0) %>%
	inner_join(
		res_dt_deseq2 %>% filter(is_significant) %>% pivot_transcript()
	)


res %>%
	mutate(`Transcript ID` = factor(`Transcript ID`, levels = .$`Transcript ID` %>% unique)) %>%
	ggplot(aes(`Transcript ID`, `tot deleterious outliers`)) + geom_point() + my_theme

res %>% pull(plot) %>% cowplot::plot_grid(plotlist = ., align = "h", ncol = 3, axis="b" )
