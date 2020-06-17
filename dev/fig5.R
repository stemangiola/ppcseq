library(tidyverse)
library(purrr)
library(magrittr)
library(ppcseq)
library(tidybulk)
library(patchwork)

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		legend.position="bottom",
		#axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

all_df =
	c(
		"dev/Mangiola_et_adipo.rds",
		"dev/GSE137631_muscle.rds",
		"dev/Atkins_et_brain.rds",
		"dev/GSE141027_lipoma.rds",
		"dev/GSE99374_CD8.rds",
		"dev/GSE151005_arabidopsis.rds"
	) %>%
	map_dfr(~.x %>% readRDS %>% mutate(source = .x) %>% rename_at(vars(1),function(x) "transcript")) %>%
	mutate(source = gsub("dev/|\\.rds", "", source)) %>%
	mutate(source = gsub("_", " ", source))


all_df_deseq2 =
	c(
		"dev/Mangiola_et_adipo_deseq2.rds",
		"dev/GSE137631_muscle_deseq2.rds",
		#"dev/Atkins_et_bra",
		"dev/GSE141027_lipoma_deseq2.rds",
		"dev/GSE99374_CD8_deseq2.rds",
		"dev/GSE151005_arabidopsis_deseq2.rds"
	) %>%
	map_dfr(~.x %>% readRDS %>% mutate(source = .x) %>% rename_at(vars(1),function(x) "transcript")) %>%
	mutate(source = gsub("dev/|\\.rds", "", source)) %>%
	mutate(source = gsub("_", " ", source))

# Median fraction
all_df %>%
	mutate(includes_outlier = `tot deleterious outliers` > 0) %>%
	group_by(source) %>%
	summarise(fraction_genes_includng_out =  sum(includes_outlier)/n(), tot_genes = n()) %>%
	summarise(median(fraction_genes_includng_out))

all_df_deseq2 %>%
	mutate(includes_outlier = `tot deleterious outliers` > 0) %>%
	group_by(source) %>%
	summarise(fraction_genes_includng_out =  sum(includes_outlier)/n(), tot_genes = n()) %>%
	summarise(median(fraction_genes_includng_out, na.rm=T))

p1 =
	all_df %>% mutate(algorithm="edgeR") %>%
	bind_rows(all_df_deseq2 %>% mutate(algorithm="DESeq2")) %>%

	# TEMPORARY
	bind_rows(
		all_df  %>%
			filter(source=="Atkins et brain") %>%
			mutate(source = "Atkins et brain deseq2") %>%
			mutate(algorithm = "DESeq2")
	) %>%

	mutate(source = str_replace(source, " deseq2", "")) %>%
	mutate(includes_outlier = `tot deleterious outliers` > 0) %>%
	group_by(source, algorithm) %>%
	summarise(fraction_genes_includng_out =  sum(includes_outlier)/n(), tot_genes = n()) %>%
	ungroup() %>%
	arrange(desc(fraction_genes_includng_out)) %>%
	mutate(source = factor(source, levels = unique(.$source))) %>%
	ggplot(aes(y=fraction_genes_includng_out, x=source, fill = algorithm)) +
	geom_bar(stat = "identity", position = "dodge") +
	scale_fill_brewer(palette = "Dark2") +
	my_theme +
	theme(
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=6),
		axis.title.x = element_text(size=8)
		)

p2 =
	all_df %>%

	# Filter non visual genes
	filter(transcript %in% c("401097") %>% `!`) %>%

	# Summarise
	group_by(source) %>%
	filter(`tot deleterious outliers` > 0) %>%
	slice(1) %>%
	ungroup() %>%
	mutate(
	plot =
		pmap(list(plot, source, transcript), ~ ..1 + theme(
			legend.position="none",
			axis.title.x  = element_blank(),
			axis.title.y  = element_blank(),
			axis.text.x = element_text(angle = 90, hjust = 1, size = 5)
		) +
			scale_fill_brewer(palette="Set1") +
			scale_size(range=c(1,2))+
			ggtitle(sprintf("%s\n%s",..2, ..3)) +
			theme(plot.title = element_text(size=8))
		)
) %>%
	pull(plot)


legend_p2 =
	all_df %>%
	filter(`ppc samples failed` == 1) %>%
	pull(plot) %>%
	map(~.x + theme(text = element_text(size=6))) %>%
	.[[1]] %>%
	cowplot::get_legend()

p  = (p1 + (p2[[1]] / p2[[4]]) + (p2[[2]] / p2[[5]]) + (p2[[3]] / p2[[6]]) +  plot_layout(nrow = 1))
	ggsave(
		"dev/Fig5.pdf",
		plot = p,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 130,
		limitsize = FALSE
	)

res %>% filter(`ppc samples failed` == 1) %>%
	pull(plot) %>% `[[` (1) + theme( legend.position = "bottom")

# How many outliers after first passage
res %>% filter(`tot deleterious outliers` == 0) %>% unnest(`sample wise data`) %>% rowwise() %>% filter(! between(value, .lower_2, .upper_2)) %>% ungroup()

# What is the n_eff I achieve on real data
res %>% attr("fit 2") %>%
	rstan::summary() %$% summary %>%
	as_tibble(rownames="par") %>%
	filter(grepl("counts_rng", par)) %>%
	ggplot(aes(n_eff)) + geom_histogram() +
	my_theme

# What is the top ranked transcript with outlier in the 6 data sets


# Read data edgeR
get_rank_outliers = function(file_name_list, outlier_df, pvalue_col, FC_col, FDR_col){
	pvalue_col = as.symbol(pvalue_col)
	FC_col = as.symbol(FC_col)
	FDR_col = as.symbol(FDR_col)

	file_name_list %>%
		map_dfr(~.x %>% readRDS %>% mutate(source = .x) %>% rename_at(vars(1),function(x) "transcript")) %>%
		mutate(source = gsub("dev/|\\.rds", "", source)) %>%
		mutate(source = gsub("_", " ", source)) %>%
		mutate(source = str_replace(source, " DT", "")) %>%

		# Build rank
		distinct(source, transcript, !!FC_col,!!pvalue_col, !!FDR_col, is_significant) %>%
		filter(is_significant) %>%
		nest(data = -source) %>%
		mutate(data = map(data, ~.x %>%
												arrange(!!pvalue_col) %>%
												mutate(rank = row_number())
		)) %>%
		unnest(data) %>%

		# Add inference
		right_join(outlier_df) %>%

		# select top
		filter(`tot deleterious outliers`>0) %>%
		group_by(source) %>%
		arrange(rank) %>%
		slice(1) %>%
		distinct(source, rank)
}



c(
	"dev/Mangiola_et_adipo_DT.rds",
	"dev/GSE137631_muscle_DT.rds",
	"dev/Atkins_et_brain_DT.rds",
	"dev/GSE141027_lipoma_DT.rds",
	"dev/GSE99374_CD8_DT.rds",
	"dev/GSE151005_arabidopsis_DT.rds"
) %>%
	get_rank_outliers(all_df, "PValue", "logFC", "FDR")

c(
	"dev/Mangiola_et_adipo_DT_deseq2.rds",
	"dev/GSE137631_muscle_DT_deseq2.rds",
	"dev/Atkins_et_brain_DT_deseq2.rds",
	"dev/GSE141027_lipoma_DT_deseq2.rds",
	"dev/GSE99374_CD8_DT_deseq2.rds",
	"dev/GSE151005_arabidopsis_DT_deseq.rds"
) %>%
	get_rank_outliers(all_df_deseq2, "pvalue", "log2FoldChange", "padj")


# Change in fold change
all_df %>%
	filter(`tot deleterious outliers`>0) %>%
	mutate(ratio = map_dbl(`sample wise data`, ~ .x %>%  distinct(slope_1, slope_2) %>% mutate(ratio = slope_1/slope_2) %>% pull(ratio))) %>%
	mutate(cohort_size = map_dbl(`sample wise data`, ~ .x %>%  distinct(sample) %>% nrow)) %>%
	group_by(source, cohort_size) %>% summarise(mean(ratio))


all_df_deseq2 %>%
	filter(`tot deleterious outliers`>0) %>%
	mutate(ratio = map_dbl(`sample wise data`, ~ .x %>%  distinct(slope_1, slope_2) %>% mutate(ratio = slope_1/slope_2) %>% pull(ratio))) %>%
	mutate(cohort_size = map_dbl(`sample wise data`, ~ .x %>%  distinct(sample) %>% nrow)) %>%
	group_by(source, cohort_size) %>% summarise(mean(ratio))


