library(tidyverse)
library(magrittr)
library(foreach)
library(ttBulk)
library(ppcSeq)
#source("https://gist.githubusercontent.com/stemangiola/dd3573be22492fc03856cd2c53a755a9/raw/e4ec6a2348efc2f62b88f10b12e70f4c6273a10a/tidy_extensions.R")


# Calculate CAPRA-S score
capra_s = function(psa,
									 sur_mar,
									 semin_ves_invas,
									 gleason,
									 extracap_ext,
									 lymph_inv,
									 sample) {

	if(c(psa, sur_mar, semin_ves_invas, gleason, extracap_ext, lymph_inv) %>% is.na %>% any) return(NA)

	# PSA
	psa_score = function(psa) {
		psa = as.numeric(psa)
		if (psa <= 6)
			0
		else if (psa > 6 & psa <= 10)
			1
		else if (psa > 10 & psa <= 20)
			2
		else if (psa > 20)
			3
	}

	# GS
	gleason_score = function(gleason) {
		gleason = as.numeric(gleason)
		if (gleason %>% gsub("\\+", "", .) <= 33)
			0
		else if (gleason %>% gsub("\\+", "", .) == 34)
			1
		else if (gleason %>% gsub("\\+", "", .) == 43)
			2
		else if (gleason %>% gsub("\\+", "", .) > 43)
			3
	}

	# Invasive
	sum(psa_score(psa) ,
			switch(sur_mar %>% as.character, "0" = 0, "1" = 2) ,
			switch(
				semin_ves_invas %>% as.character,
				"0" = 0,
				"1" = 2
			) ,
			gleason_score(gleason),
			switch(
				extracap_ext %>% as.character,
				"0" = 0,
				"1" = 2
			) ,
			switch(lymph_inv %>% as.character, "0" = 0, "1" = 2),
			na.rm = T)
}

get_TCGA_prostate_clinical_annotaton = function() {
	mycgds = "http://www.cbioportal.org/" %>% cgdsr::CGDS()

	foreach(
		fi =
			mycgds %>%
			cgdsr::getCancerStudies() %>%
			filter(grepl("TCGA", name)) %>%
			pull(cancer_study_id) %>%
			grep("prad", ., value = T),
		.combine = bind_rows
	) %do% {
		print(fi)
		mycaselist = cgdsr::getCaseLists(mycgds, fi) %>% as_tibble() %>% filter(grepl("RNA Seq", case_list_name))
		print(ncol(mycaselist))
		cgdsr::getClinicalData(mycgds, caseList = mycaselist[1, ] %>% pull(case_list_id)) %>%
			as_tibble(rownames = "sample") %>%
			mutate_all(funs(as.character)) %>%
			mutate(cancer_study_id = fi)
	} %>%
		# Clean from absent data
		filter(!(DFS_MONTHS == "" | is.na(DFS_MONTHS))) %>%
		# Select duplicate more up to date
		group_by(sample) %>%
		arrange(desc(DFS_MONTHS)) %>%
		filter(row_number() == 1) %>%
		ungroup() %>%
		separate(sample, c("t1", "t2", "t3"), sep = "\\.") %>%
		unite(sample, c("t1", "t2", "t3"), sep = "-")  %>%

		#3 Calculate CAPRA-S score
		unite(
			GLEASON_SCORE_FORMATTED,
			c("GLEASON_PATTERN_PRIMARY", "GLEASON_PATTERN_SECONDARY"),
			sep = "",
			remove = F
		) %>%
		mutate(sur_mar = ifelse(RESIDUAL_TUMOR == "R0", 0, 1)) %>%
		mutate(LYMPH_NODES_EXAMINED_HE_COUNT = ifelse(
			is.na(LYMPH_NODES_EXAMINED_HE_COUNT),
			0,
			LYMPH_NODES_EXAMINED_HE_COUNT
		)) %>%
		#filter(!is.na(PSA_MOST_RECENT_RESULTS) &	!is.na(GLEASON_SCORE_FORMATTED)) %>%
		rowwise() %>%
		mutate(
			`CAPRA-S` = capra_s(
				psa = PSA_MOST_RECENT_RESULTS,
				sur_mar = sur_mar,
				semin_ves_invas = 0,
				gleason = GLEASON_SCORE_FORMATTED,
				extracap_ext = 0,
				lymph_inv = LYMPH_NODES_EXAMINED_HE_COUNT,
				sample = sample
			)
		) %>%
		ungroup() %>%

		# Format further
		distinct(sample,
						 DFS_MONTHS,
						 DFS_STATUS,
						 `CAPRA-S`) %>%
		mutate(is_recurred = ifelse(DFS_STATUS == "Recurred/Progressed", T, F)) %>%
		mutate(DFS_MONTHS = as.numeric(DFS_MONTHS))
}

hkg = read_csv("dev/hk_600.txt", col_names = FALSE) %>% pull(1)

TCGA_tbl = read_csv(
	#"~/unix3XX/PhD/deconvolution/TCGA_all_cancers/tibble_TCGA_files/Prostate_Adenocarcinoma_Primary_Tumor.csv"
	"~/PhD/deconvolution/TCGA_all_cancers/tibble_TCGA_files/Prostate_Adenocarcinoma_Primary_Tumor.csv"
) %>%
	separate(sample, c("t1", "t2", "t3"), sep = "-") %>%
	unite(sample, c("t1", "t2", "t3"), sep = "-") %>%
	left_join(get_TCGA_prostate_clinical_annotaton()) %>%
	mutate_if(is.character, as.factor) %>%

	# Select little samples
	#inner_join((.) %>% distinct(sample) %>% head(n = 20)) %>%

	# Prepare data frame
	tidyr::separate(ens_iso, c("ens", "iso"), sep = "\\.") %>%
	annotate_symbol(ens) %>%
	filter(transcript %>% is.na %>% `!`) %>%
	ttBulk(
		sample,
		transcript,
		 `read count`
	) %>%
	aggregate_duplicates() %>%
	mutate(`read count` = `read count` %>% as.integer) %>%
	mutate_if(is.character, as.factor) %>%

	# setup house keeping genes
	mutate(`house keeping` = transcript %in% hkg) %>%
	mutate(PValue = ifelse(!`house keeping`, 0, 1)) %>%

	# genes to check
	left_join( (.) %>% distinct(transcript) %>% mutate(run = sample( 1:6, size = n(), replace = T)))


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
		just_discovery = T,
		approximate_posterior_inference = F,
		approximate_posterior_analysis = T,
		cores = 30,
		additional_parameters_to_save = c("intercept", "sigma_raw", "sigma_intercept", "sigma_slope", "sigma_sigma")
	)

saveRDS(res, file="~/PostDoc/temp_res_ppcSeq.RData")

# res = readRDS("~/PostDoc/temp_res_ppcSeq.RData")
#
# x = res %>% distinct(transcript, G, S, `read count`, sample) %>%
#
# 	# Attach paramers NB
# 	left_join(
# 		res %>% attr("fit") %>%
# 			rstan::summary() %$% summary %>%
# 			as_tibble(rownames="par") %>%
# 			filter(grepl("intercept|sigma_raw", par)) %>%
# 			tidyr::separate(par, c("par", "G"), sep="[\\[,\\]]", extra="drop") %>%
# 			mutate(G = G %>% as.integer) %>%
# 			distinct(G, par, mean) %>%
# 			tidyr::spread(par, mean)
# 	) %>%
#
# 	# Attach paramers exposure samples
# 	left_join(
# 		res %>% attr("fit") %>%
# 			rstan::summary() %$% summary %>%
# 			as_tibble(rownames="par") %>%
# 			filter(grepl("exposure_rate", par)) %>%
# 			tidyr::separate(par, c("par", "S"), sep="[\\[,\\]]", extra="drop") %>%
# 			mutate(S = S %>% as.integer) %>%
# 			distinct(S, par, mean) %>%
# 			tidyr::spread(par, mean)
# 	) %>%
#
# 	# Normalised read counts
# 	mutate(`read count normalised bayes` = `read count` / exp(exposure_rate))	%>%
#
# 	# Calculate quantiles
# 	# QQ plots
# 	do_parallel_start(20, "transcript") %>%
# 	do({
#
# 		`%>%` = magrittr::`%>%`
# 		library(tidyverse)
# 		library(magrittr)
#
# 		(.) %>%
# 			group_by(transcript) %>%
# 			do(
# 				(.) %>%
# 					arrange(`read count normalised bayes`) %>%
# 					mutate(
# 						predicted_NB =
#
# 							qnbinom(
#
# 								# If 1 sample, just use median
# 								switch(	((.) %>% nrow>1) %>% `!` %>% `+` (1), ppoints(`read count normalised bayes`), 0.5	),
# 								size=.$sigma_raw %>% unique %>% exp %>% `^` (-1),
# 								mu=.$intercept %>% unique %>% exp
# 							)
#
# 					)
# 			) %>%
# 			ungroup()
# 	}) %>%
# 	do_parallel_end()
#
# xx = x %>%
# 	create_tt_from_tibble_bulk(sample_column = sample, transcript_column = transcript, counts_column = `read count`) %>%
# 	group_by(transcript) %>%
# 	do(
# 		(.) %>%
# 			add_rotated_dimensions(
# 				dimension_1_column = predicted_NB,
# 				dimension_2_column = `read count normalised bayes` ,
# 				rotation_degrees = -45
# 			)
# 	) %>%
# 	ungroup() %>%
# 	nest(data = c(-transcript))
#
# xx %>% slice(60) %>% unnest(data) %>% ggplot(aes(x = `predicted_NB rotated -45` +1, y = `read count normalised bayes rotated -45` +1)) + geom_point() + scale_x_log10() + scale_y_log10()
#
# mutate(`error of log` = (log(`read count normalised bayes` + 1) - log(predicted_NB + 1)) ) %>%
#


	# arrange(`read count normalised bayes`) %>%
	# tidyr::nest(data = c( S, `read count`, intercept, sigma_raw ,exposure_rate,`read count normalised bayes`)) %>%
	#
	# mutate(
	# 	theoretical_NB = map(
	# 		data,
	# 		~ qnbinom(
	# 			ppoints(.x$`read count normalised bayes`),
	# 			size= %>% unique %>% exp %>% `^` (-1),
	# 			mu=..3 %>% unique %>% exp
	# 		)
	# 	)
	# )




