library(tidyverse)
library(ppcseq)

TCGA_tbl_subsample %>%
	mutate(do_check = (!`house keeping`) & run==1) %>%
	mutate(`read count` = `read count` %>% as.integer) %>%
	identify_outliers(
		.significance = PValue,
		.do_check = do_check,
		.abundance = `read count`,
		.sample = sample,
		.transcript = transcript,
		percent_false_positive_genes = "1%",
		tol_rel_obj = 0.01,
		cores = 2,
		approximate_posterior_inference = T,
		approximate_posterior_analysis = T,
		pass_fit = F
	)


system.time({

TCGA_tbl_subsample %>%
	inner_join(
		(.) %>% distinct(sample) %>% slice(1:20)
	) %>%
	select(-run) %>%
	left_join(
		(.) %>%
			distinct(transcript) %>%
			mutate(run = sample(1:40, replace = T, size=n()))
	) %>%

	mutate(do_check = (!`house keeping`) & run==1) %>%
	mutate(`read count` = `read count` %>% as.integer) %>%
	identify_outliers(
		.significance = PValue,
		.do_check = do_check,
		.abundance = `read count`,
		.sample = sample,
		.transcript = transcript,
		percent_false_positive_genes = "1%",
		tol_rel_obj = 0.01,
		cores = 2,
		approximate_posterior_inference = T,
		approximate_posterior_analysis = T,
		pass_fit = F
	)
})
