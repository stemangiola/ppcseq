library(tidyverse)
library(magrittr)
library(ppcseq)
library(tidybulk)

df =
	dir("dev/GSE151005_RAW/", full.names = T) %>%
	map(~ .x %>% read_tsv) %>%
	do.call(bind_cols, .) %>%
	select(gene_id...1, contains("Col"), contains("NP")) %>%
	dplyr::rename(gene_id = gene_id...1) %>%
	gather(sample, count, -gene_id) %>%
	tidyr::extract(col = sample, into = c( "type", "batch"), "(.+)_(.+)", remove=F) %>%
	tidybulk(sample, gene_id, count) %>%
	identify_abundant(factor_of_interest = type) %>%
	scale_abundance() %>%
	mutate(batch_B = batch == "B") %>%

	{
		((.) %>%
			reduce_dimensions(method="MDS", action = "get") %>%
			ggplot(aes(Dim1, Dim2, color=type, shape=batch_B)) + geom_point()) %>% print
		(.)
	}

res_dt =
	df%>%
	test_differential_abundance(~ type + batch_B, method = "edgeR_likelihood_ratio") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )

res_dt %>% saveRDS("dev/GSE151005_arabidopsis_DT.rds", compress = "gzip")

res_dt_robust =
	df%>%
	test_differential_abundance(~ type + batch_B, method = "edgeR_robust_likelihood_ratio") %>%
	filter(.abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )

res_dt_robust %>% saveRDS("dev/GSE151005_arabidopsis_DT_robust.rds", compress = "gzip")

res_dt_deseq =
	df %>%
	test_differential_abundance(~ type + batch_B, method = "deseq2") %>%
	filter(.abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 )

res_dt_deseq %>% saveRDS("dev/GSE151005_arabidopsis_DT_deseq.rds", compress = "gzip")

res =
	res_dt %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type + batch,
		.sample = sample,
		.transcript = gene_id,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("dev/GSE151005_arabidopsis.rds", compress = "gzip")

res_robust =
	res_dt_robust %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type + batch,
		.sample = sample,
		.transcript = gene_id,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_robust %>% select(-plot, -`sample wise data`) %>% saveRDS("dev/GSE151005_arabidopsis_robust.rds", compress = "gzip")

res_deseq =
	res_dt_deseq %>%
	filter(!lowly_abundant) %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type + batch,
		.sample = sample,
		.transcript = gene_id,
		.significance = pvalue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_deseq %>% saveRDS("dev/GSE151005_arabidopsis_deseq2.rds", compress = "gzip")

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
