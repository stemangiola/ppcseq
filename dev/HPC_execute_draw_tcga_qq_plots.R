library(tidyverse)
library(magrittr)
library(ttBulk)
library(ppcSeq)

#TCGA_tbl = readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/temp_TCGA_tbl.RData")
#load("dev/TCGA_tbl.RData")

# Interctept only

foreach(r = 1:6) %do% {
	res =
		TCGA_tbl %>%

		# Filter
		filter(`CAPRA-S` %>% is.na %>% `!`) %>%
		separate(sample, c("data base", "laboratory", "patient"), sep="-", remove = F) %>%
		inner_join(
			(.) %>% distinct(sample, laboratory) %>% count(laboratory) %>% filter(n >= 8) %>% select(-n)
		) %>%
		mutate(risk = `CAPRA-S` <= 3) %>%

		# Do check
		mutate(do_check = (!`house keeping`) & run==r) %>%

		ppc_seq(
			significance_column = PValue,
			do_check_column = do_check,
			value_column = `read count`,
			percent_false_positive_genes = "5%",
			sample_column = sample,
			gene_column = transcript,
			pass_fit = T,
			tol_rel_obj = 0.01,
			just_discovery = T, full_bayes = F,
			cores = 10,
			additional_parameters_to_save = c("sigma_raw", "alpha")
		)

	saveRDS(res, file=sprintf("dev/draw_qq_TCGA_no_covariates_%s.rds", r))

}

# Analysis of confounders

TCGA_tbl.MDS =
	TCGA_tbl %>%
	create_ttBulk(sample_column = sample, transcript_column = transcript, counts_column = `read count`) %>%
	normalise_counts() %>%
	reduce_dimensions(value_column = `read count normalised`, method = "MDS", components = 1:10)

TCGA_tbl.MDS %>%
	select(sample,  `CAPRA-S`, contains("Dimension")) %>%
	distinct() %>%
	filter(`CAPRA-S` %>% is.na %>% `!`) %>%
	mutate(risk = `CAPRA-S` <= 3) %>%
	mutate(`CAPRA-S` = `CAPRA-S` %>% as.integer %>% as.factor)  %>%
	separate(sample, c("data base", "laboratory", "patient"), sep="-", remove = F) %>%
	inner_join(
		(.) %>% distinct(sample, laboratory) %>% count(laboratory) %>% filter(n >= 8) %>% select(-n)
	) %>%
	select(contains("Dimension"), everything()) %>%
	distinct() %>%
	GGally::ggpairs(columns = 1:6, ggplot2::aes(colour=`laboratory`))

# Model with covariates

foreach(r = 1:6) %do% {
	res =
		TCGA_tbl %>%

		# Filter
		filter(`CAPRA-S` %>% is.na %>% `!`) %>%
		separate(sample, c("data base", "laboratory", "patient"), sep="-", remove = F) %>%
		inner_join(
			(.) %>% distinct(sample, laboratory) %>% count(laboratory) %>% filter(n >= 8) %>% select(-n)
		) %>%
		mutate(risk = `CAPRA-S` <= 3) %>%

		# Do check
		mutate(do_check = (!`house keeping`) & run==r) %>%

		ppc_seq(
			~ risk + laboratory,
			significance_column = PValue,
			do_check_column = do_check,
			value_column = `read count`,
			percent_false_positive_genes = "5%",
			sample_column = sample,
			gene_column = transcript,
			pass_fit = T,
			tol_rel_obj = 0.01,
			just_discovery = T, full_bayes = F,
			cores = 10,
			additional_parameters_to_save = c("sigma_raw", "alpha")
		)

	saveRDS(res, file=sprintf("dev/draw_qq_TCGA_with_covariates_%s.rds", r))

}


x = readRDS("dev/draw_qq_TCGA_with_covariates_1.rds")

x %>%
	attr("fit") %>%
	rstan::summary() %$% summary %>%
	as_tibble(rownames="par") %>%
	filter(grepl("alpha", par)) %>%
	separate(par, c("par", "C", "G"), sep="\\[|,|\\]", extra = "drop") %>%
	mutate(C = C %>% as.integer, G = G %>% as.integer) %>%
	arrange(G, C) %>%
	left_join(x %>% distinct(G, transcript), by="G") %>%
	select(transcript, G, C, mean) %>%
	spread(C, mean)
