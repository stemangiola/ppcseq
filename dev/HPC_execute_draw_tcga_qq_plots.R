library(tidyverse)
library(magrittr)
library(tidyTranscriptomics)
library(ppcSeq)

#TCGA_tbl = readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/temp_TCGA_tbl.RData")
load("dev/TCGA_tbl.RData")

res =
	TCGA_tbl %>%
	mutate(do_check = (!`house keeping`) & run==1) %>%

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
		additional_parameters_to_save = c("intercept", "sigma_raw")
	)

saveRDS(res, file="/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/temp_res_ppcSeq.RData")


#readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/temp_res_ppcSeq.RData")
