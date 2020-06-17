library(tidyverse)
library(purrr)
library(magrittr)
library(ppcseq)
library(tidybulk)


res = read_csv("dev/GSE150019_counts_TGFb.csv")  %>%
	gather(sample, count, -gene) %>%
	tidyr::extract(
		col = sample,
		into = c("type"),
		"[A-Z0-9]+_([A-Za-z0-9]+)",
		remove=F,
		convert = T
	) %>%
	tidybulk(sample, gene, count) %>%
	scale_abundance(factor_of_interest = type, minimum_counts = 3, minimum_proportion = 2/3) %>%
	mutate(is_sp1 = grepl("SP1", sample)) %>%
	adjust_abundance(~type + is_sp1) %>%

	{
		((.) %>%
		 	filter(count_scaled_adjusted %>% is.na %>% `!`) %>%
			reduce_dimensions(.abundance = count_scaled_adjusted, method="MDS", action = "get") %>%
			ggplot(aes(Dim1, Dim2, color=type, label=sample)) + geom_point() + geom_text()) %>% print
		(.)
	} %>%
	test_differential_abundance(~ type + is_sp1, minimum_counts = 3, minimum_proportion = 2/3) %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 ) %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ type + is_sp1,
		.sample = sample,
		.transcript = gene,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("myeloma_IFNg.rds", compress = "gzip")

res %>%
	mutate(`Transcript ID` = factor(`Transcript ID`, levels = .$`Transcript ID` %>% unique)) %>%
	ggplot(aes(`Transcript ID`, `tot deleterious outliers`)) + geom_point() + my_theme

res %>% pull(plot) %>% cowplot::plot_grid(plotlist = ., align = "h", ncol = 3, axis="b" )
