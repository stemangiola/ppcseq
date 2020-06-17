library(tidyverse)
library(magrittr)
library(ppcseq)
library(tidybulk)


load("~/third_party_analyses/PUBLISHED_Atkins_et_al_2019_brain_ryan/global_analysis/counts_01_01_2016.RData")
colnames(d) = c("LN229.2neg_S3" ,   "LN229.2pos_S2"  ,  "LN229.3neg_S5"  ,  "LN229.3pos_S4" ,   "LN229.4neg_S7"   , "LN229.4pos_S6", "b", "u")
d_1 = d$counts
load("~/third_party_analyses/PUBLISHED_Atkins_et_al_2019_brain_ryan/global_analysis/counts_18_05_2016.RData")
colnames(d) = gsub("alignment_hg19\\.|\\.sam", "", colnames(d))
d_2 = d$counts
load("~/third_party_analyses/PUBLISHED_Atkins_et_al_2019_brain_ryan/global_analysis/counts_20_06_2016.RData")
colnames(d) = c("U118.A.pos_S2",  "U118.A.neg_S3",  "U118.B.pos_S4",  "U118.B.neg_S5",  "U118.pos_S6",  "U118.neg_S7")
d_3 = d$counts

annot <- read.csv("~/third_party_analyses/PUBLISHED_Atkins_et_al_2019_brain_ryan/global_analysis/annot.csv")
colnames(annot) <- c("sample", "run", "cellLine", "batch", "treatment")
annot$batch = factor(annot$batch)




df =
	as.data.frame(d_1) %>% as_tibble(rownames = "entrez") %>%
	left_join(
		as.data.frame(d_2) %>% as_tibble(rownames = "entrez")
	) %>%
	left_join(
		as.data.frame(d_3) %>% as_tibble(rownames = "entrez")
	) %>%
	gather(sample, count, -entrez) %>%
	left_join(annot)  %>%
	filter(treatment %>% is.na %>% `!`) %>%
	tidybulk(sample, entrez, count) %>%
	scale_abundance(factor_of_interest = treatment) %>%
	adjust_abundance(~  treatment + cellLine) %>%

	{
		((.) %>%
		 	filter(count_scaled_adjusted %>% is.na %>% `!`) %>%
			reduce_dimensions(.abundance = count_scaled_adjusted, method="MDS", action = "get") %>%
			ggplot(aes(Dim1, Dim2, color=treatment, shape=batch)) + geom_point()) %>% print
		(.)
	}

res_dt = df %>%
	test_differential_abundance(~ treatment + cellLine, method = "edgeR_likelihood_ratio") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = FDR < 0.05 )

res_dt %>% saveRDS("dev/Atkins_et_brain_DT.rds", compress = "gzip")

res_dt_deseq2 = df %>%
	test_differential_abundance(~ treatment + cellLine, method = "deseq2") %>%
	filter(!lowly_abundant) %>%
	mutate(count = as.integer(count)) %>%
	mutate(is_significant = !is.na(padj) & padj < 0.05 )

res_dt_deseq2 %>% saveRDS("dev/Atkins_et_brain_DT_deseq2.rds", compress = "gzip")


res =
	res_dt %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ treatment + cellLine,
		.sample = sample,
		.transcript = entrez,
		.significance = PValue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res %>% saveRDS("dev/Atkins_et_brain.rds", compress = "gzip")

res_deseq2 =
	res_dt_deseq2 %>%
	filter(padj %>% is.na %>% `!`) %>%

	# PPCSEQ
	ppcseq::identify_outliers(
		formula = ~ treatment + cellLine,
		.sample = sample,
		.transcript = entrez,
		.significance = pvalue,
		.do_check  = is_significant,
		.abundance = count,
		percent_false_positive_genes = 1,
		cores=30,
		pass_fit = T
	)

res_deseq2 %>% saveRDS("dev/Atkins_et_brain_deseq2.rds", compress = "gzip")


res %>%
	mutate(`Transcript ID` = factor(`Transcript ID`, levels = .$`Transcript ID` %>% unique)) %>%
	ggplot(aes(`Transcript ID`, `tot deleterious outliers`)) + geom_point() + my_theme

res %>% pull(plot) %>% cowplot::plot_grid(plotlist = ., align = "h", ncol = 3, axis="b" )
