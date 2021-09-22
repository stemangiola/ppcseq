
#' identify_outliers main
#'
#' @description This function runs the data modeling and statistical test for the hypothesis that a transcript includes outlier biological replicate.
#'
#' \lifecycle{maturing}
#'
#' @importFrom tibble as_tibble
#' @import dplyr
#' @importFrom tidyr spread
#' @import tidybayes
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom benchmarkme get_ram
#' @importFrom magrittr multiply_by
#' @importFrom magrittr equals
#' @import edgeR
#' @importFrom stats sd
#' @importFrom purrr map_chr
#'
#' @param .data A tibble including a transcript name column | sample name column | read counts column | covariate columns | Pvalue column | a significance column
#' @param formula A formula. The sample formula used to perform the differential transcript abundance analysis
#' @param .sample A column name as symbol. The sample identifier
#' @param .transcript A column name as symbol. The transcript identifier
#' @param .abundance A column name as symbol. The transcript abundance (read count)
#' @param .significance A column name as symbol. A column with the Pvalue, or other significance measure (preferred Pvalue over false discovery rate)
#' @param .do_check A column name as symbol. A column with a boolean indicating whether a transcript was identified as differentially abundant
#' @param percent_false_positive_genes A real between 0 and 100. It is the aimed percent of transcript being a false positive. For example, percent_false_positive_genes = 1 provide 1 percent of the calls for outlier containing transcripts that has actually not outliers.
#' @param approximate_posterior_inference A boolean. Whether the inference of the joint posterior distribution should be approximated with variational Bayes It confers execution time advantage.
#' @param approximate_posterior_analysis A boolean. Whether the calculation of the credible intervals should be done semi-analytically, rather than with pure sampling from the posterior. It confers execution time and memory advantage.
#' @param how_many_negative_controls An integer. How many transcript from the bottom non-significant should be taken for inferring the mean-overdispersion trend.
#' @param draws_after_tail An integer. How many draws should on average be after the tail, in a way to inform CI.
#' @param save_generated_quantities A boolean. Used for development and testing purposes
#' @param additional_parameters_to_save A character vector. Used for development and testing purposes
#' @param cores An integer. How many cored to be used with parallel calculations.
#' @param pass_fit A boolean. Used for development and testing purposes
#' @param do_check_only_on_detrimental A boolean. Whether to test only for detrimental outliers (same direction as the fold change). It allows to test for less transcript/sample pairs and therefore higher the probability threshold.
#' @param tol_rel_obj A real. Used for development and testing purposes
#' @param just_discovery A boolean. Used for development and testing purposes
#' @param seed An integer. Used for development and testing purposes
#' @param adj_prob_theshold_2 A boolean. Used for development and testing purposes
#' @param return_fit A boolean
#'
#' @return A nested tibble `tbl` with transcript-wise information: `sample wise data` | plot | `ppc samples failed` | `tot deleterious outliers`
#'
#' @examples
#'
#' library(dplyr)
#'
#' data("counts")
#'
#' if(Sys.info()[['sysname']] == "Linux")
#' result =
#'   counts %>%
#'   dplyr::mutate(  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), TRUE, FALSE) ) %>%
#'	 ppcseq::identify_outliers(
#'		formula = ~ Label,
#'		sample, symbol, value,
#'		.significance = PValue,
#'		.do_check  = is_significant,
#'		percent_false_positive_genes = 1,
#'		tol_rel_obj = 0.01,
#'		approximate_posterior_inference =TRUE,
#'		approximate_posterior_analysis =TRUE,
#'		how_many_negative_controls = 50,
#'		cores=1
#'	)
#'
#' @export
#'
identify_outliers = function(.data,
														 formula = ~ 1,
														 .sample,
														 .transcript,
														 .abundance,
														 .significance,
														 .do_check,
														 percent_false_positive_genes = 1,
														 how_many_negative_controls = 500,

														 approximate_posterior_inference = TRUE,
														 approximate_posterior_analysis = TRUE,
														 draws_after_tail = 10,

														 save_generated_quantities = FALSE,
														 additional_parameters_to_save = c(),  # For development purpose
														 cores = detect_cores(), # For development purpose,
														 pass_fit = FALSE,
														 do_check_only_on_detrimental = length(parse_formula(formula)) > 0,
														 tol_rel_obj = 0.01,
														 just_discovery = FALSE,
														 seed = sample(seq_len(length.out=999999), size = 1),
														 adj_prob_theshold_2 = NULL,
														 return_fit = FALSE
) {
	# Prepare column same enquo
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.significance = enquo(.significance)
	.do_check = enquo(.do_check)

	# Get factor of interest
	#factor_of_interest = ifelse(parse_formula(formula) %>% length %>% `>` (0), parse_formula(formula)[1], "")

	# Check if columns exist
	check_columns_exist(.data, !!.sample, !!.transcript, !!.abundance, !!.significance, !!.do_check)

	# Check if any column is NA or null
	check_if_any_NA(.data, !!.sample, !!.transcript, !!.abundance, !!.significance, !!.do_check, parse_formula(formula))

	# Check if I have any genes to check
	if(.data %>%	filter(!!.do_check) %>% nrow %>% equals(0)){
		warning("ppcseq says: There are not transcripts with the category .to_check. NULL is returned.")
		return(tibble(
			a = "a",
			b = list(),
			c = 1L,
			d = 1L
		) %>%
			slice(0) %>%
			setNames(c(quo_name(.transcript), "sample wise data", "ppc samples failed", "tot deleterious outliers")))
	}


	# Check is testing environment is supported
	if (approximate_posterior_inference &	save_generated_quantities)
		stop("Variational Bayes does not support tidybayes needed for save_generated_quantities, use sampling")

	# Check percent FP input
	if (percent_false_positive_genes %>% is.na |
			!(percent_false_positive_genes %>% between(0, 100)))
		stop("percent_false_positive_genes must be a string from > 0% to < 100%")

	# For reference MPI inference
	# Check if all transcripts are non NA
	if (.data %>% filter(!!.transcript %>% is.na) %>% nrow > 0)
		stop("There are NAs in the .transcript. Please filter those records")

	# Check if the counts column is an integer
	if (.data %>% pull(!!.abundance) %>% is("integer") %>% not())
		stop(
			sprintf(
				"The column %s must be of class integer. You can do as mutate(`%s` = `%s` %%>%% as.integer)",
				quo_name(.abundance),
				quo_name(.abundance),
				quo_name(.abundance)
			)
		)

	# Calculate the adj_prob_theshold
	if(adj_prob_theshold_2 %>% is.null)
		adj_prob_theshold_2 =
		percent_false_positive_genes / 100 /
		(.data %>% distinct(!!.sample) %>% nrow) *
		ifelse(do_check_only_on_detrimental, 2, 1)

	# The first pasage is at least 2 times more permissive than the second
	adj_prob_theshold_1  = 0.05 %>% max(adj_prob_theshold_2*2)

	# Calculate adj_prob_theshold
	how_many_posterior_draws_1 =  draws_after_tail %>% divide_by(adj_prob_theshold_1) %>% max(1000) # I want 5 draws in the tail
	how_many_posterior_draws_2 =  draws_after_tail %>% divide_by(adj_prob_theshold_2) %>% max(1000) # I want 5 draws in the tail

	# If too many draws required revert to approximation of CI
	if(approximate_posterior_analysis %>% is.null){
		if(how_many_posterior_draws_2 > 20000) {
			writeLines(sprintf("The number of draws needed to calculate the CI from the posterior would be larger than %s. To avoid impractical computation times, the calculation of the CI will be based on the mean, exposure and overdisperison posteriors.", how_many_posterior_draws_2))
			approximate_posterior_analysis = TRUE
		} else approximate_posterior_analysis = FALSE
	}

	# Check if enough memory for full draw
	available_memory = ifelse(
		.Platform$OS.type == "windows",
		shell('systeminfo | findstr Memory', intern = TRUE)[1] %>% gsub(",", "", .) %>% gsub(".*?([0-9]+).*", "\\1", .) %>% as.integer %>% divide_by(1000) %>% multiply_by(1e+9),
		get_ram() %>% as.numeric() %>% multiply_by(1e+9)
	)

	required_memory = ifelse(
		approximate_posterior_inference %>% `!`,
		1.044e+06 + how_many_posterior_draws_2 * 3.777e-02, # Regression taken from performances.R
		1.554e+06 + how_many_posterior_draws_2 * 7.327e-02  # Regression taken from performances.R
	)
	if(required_memory > available_memory & !approximate_posterior_analysis) {
		warning("
						You don't have enough memory to model the posterior distribution with MCMC draws.
						Therefore the parameter approximate_posterior_analysis was set to TRUE
		")
		approximate_posterior_analysis = TRUE
	}


	# distinct_at is not released yet for dplyr, thus we have to use this trick
	my_df <-
		.data %>%
		mutate(do_check___ = !!.do_check) %>%
		format_input(
		formula,
		!!.sample,
		!!.transcript,
		!!.abundance,
		do_check___,
		!!.significance,
		how_many_negative_controls
	)

	# Create design matrix
	X = create_design_matrix(my_df, formula, !!.sample)

	C = X %>% ncol

	# Prior info
	lambda_mu_mu = 5.612671

	# Scale dataset
	my_df_scaled =
		my_df %>%
		.identify_abundant(!!.sample,!!.transcript,!!.abundance) %>%
		get_scaled_counts_bulk(!!.sample,!!.transcript,!!.abundance) %>%
		left_join(my_df, by=quo_name(.sample)) %>%
		dplyr::mutate(!!as.symbol(sprintf("%s_scaled",  quo_name(.abundance))) := !!.abundance * multiplier)

	# Build better scales for the inference
	exposure_rate_multiplier =
		my_df_scaled %>%
		distinct(!!.sample, TMM, multiplier) %>%
		mutate(l = multiplier %>% log) %>%
		summarise(l %>% sd) %>%
		pull(`l %>% sd`)

	# Build better scales for the inference
	intercept_shift_scale =
		my_df_scaled %>%
		mutate(cc =
					 	!!as.symbol(sprintf(
					 		"%s_scaled",  quo_name(.abundance)
					 	)) %>%
					 	`+` (1) %>% log) %>%
		summarise(shift = cc %>% mean, scale = cc %>% sd) %>%
		as.numeric

	# Run the first discovery phase with permissive false discovery rate
	res_discovery =
		my_df %>%
		do_inference(
			formula,!!.sample ,!!.transcript ,!!.abundance ,!!.significance ,!!.do_check,
			approximate_posterior_inference,
			approximate_posterior_analysis = FALSE,
			C,
			X,
			lambda_mu_mu,
			cores,
			exposure_rate_multiplier,
			intercept_shift_scale,
			additional_parameters_to_save,
			adj_prob_theshold  = adj_prob_theshold_1,
			how_many_posterior_draws = how_many_posterior_draws_1,
			pass_fit = TRUE,
			tol_rel_obj = tol_rel_obj,
			write_on_disk = write_on_disk,
			seed = seed
		)

	# For building some figure I just need the discovery run, return prematurely
	if(just_discovery) return(res_discovery %>% filter(.variable == "counts_rng"))

	# Columns of counts to be ignored from the inference
	to_exclude =
		res_discovery %>%
		filter(`.variable` == "counts_rng") %>%
		ifelse_pipe(
			do_check_only_on_detrimental,
			~ .x %>% filter(`deleterious outliers`),
			~ .x %>% filter(!ppc)
		) %>%
		distinct(S, G, .lower, .upper)

	# Claculate how many popential non NB transcript I should check
	how_namy_to_exclude = to_exclude %>% nrow

	# Get the credible intervals for which account in the truncated NB model
	truncation_values =
		res_discovery %>%
		filter(`.variable` == "counts_rng") %>%
		distinct(S, G, .lower, .upper) %>%
		mutate(`.lower` = `.lower` %>% as.integer,
					 `.upper` = `.upper` %>% as.integer)

	# Get the inferred values from first model to possibly use them in the second model as priors
	prior_from_discovery =
		res_discovery %>%
		filter(`.variable` != "counts_rng") %>%
		select(`.variable`, S, G, mean, sd)

	# Run the second test phase with the user selected false discovery rate
	res_test =
		my_df %>%
		do_inference(
			formula,!!.sample ,!!.transcript ,!!.abundance ,!!.significance ,!!.do_check,
			approximate_posterior_inference,
			approximate_posterior_analysis,
			C,
			X,
			lambda_mu_mu,
			cores,
			exposure_rate_multiplier,
			intercept_shift_scale,
			additional_parameters_to_save,
			adj_prob_theshold = adj_prob_theshold_2, # If check only deleterious is one side test
			# * 2 because we just test one side of the distribution
			how_many_posterior_draws = how_many_posterior_draws_2,
			pass_fit = pass_fit,
			to_exclude = to_exclude,
			save_generated_quantities = save_generated_quantities,
			tol_rel_obj = tol_rel_obj,
			truncation_compensation = 0.7352941, # Taken by approximation study
			write_on_disk = write_on_disk,
			seed = seed
		)

	# Merge results and return
	merge_results(

		# Calculate CI 2 for discovery for plotting
		res_discovery,

		# %>%
		# 	left_join(
		# 		(.) %>%
		# 			attr("fit") %>%
		# 			fit_to_counts_rng_approximated(adj_prob_theshold_2, how_many_posterior_draws_2, truncation_compensation = 0.7352941, cores) %>%
		# 			select(S, G, .lower_1 = .lower, .upper_1 = .upper)
		# 	),

		res_test,
		formula,
		!!.transcript,
		!!.abundance,
		!!.sample,
		do_check_only_on_detrimental
	) %>%

		# If return_fit
		when(
			return_fit ~ (.) %>%
				add_attr(res_discovery %>% attr("fit"), "fit 1") %>%
				add_attr(res_test %>% attr("fit"), "fit 2"),
			~ (.)
		) %>%

		# Add total draws
		add_attr(res_test %>% attr("total_draws"), "total_draws") %>%
		add_attr(quo_name(.transcript), "transcript_column") %>%
		add_attr(quo_name(.abundance), "abundance_column") %>%
		add_attr(quo_name(.sample), "sample_column") %>%
		add_attr(quo_name(formula), "formula")

}

#' plot_credible interval for theoretical data distributions
#'
#' @description Plot the data along the theoretical data distribution.
#'
#' @importFrom tibble as_tibble
#' @importFrom stats as.formula
#'
#' @param .data The tibble returned by identify_outliers
#'
#' @return A tibble with an additional `plot` column
#'
#' @examples
#'
#' library(dplyr)
#'
#' data("counts")
#'
#' if(Sys.info()[['sysname']] == "Linux"){
#' result =
#'   counts %>%
#'   dplyr::mutate(  is_significant = ifelse(symbol %in% c("SLC16A12", "CYP1A1", "ART3"), TRUE, FALSE) ) %>%
#'	 ppcseq::identify_outliers(
#'		formula = ~ Label,
#'		sample, symbol, value,
#'		.significance = PValue,
#'		.do_check  = is_significant,
#'		percent_false_positive_genes = 1,
#'		tol_rel_obj = 0.01,
#'		approximate_posterior_inference =TRUE,
#'		approximate_posterior_analysis =TRUE,
#'		how_many_negative_controls = 50,
#'		cores=1
#'	)
#'
#' result_plot = result %>% plot_credible_intervals()
#' }
#'
#' @export
#'
plot_credible_intervals = function(.data){

	.transcript = .data %>% attr("transcript_column")
	.abundance = .data %>% attr("abundance_column")
	.sample = .data %>% attr("sample_column")
	formula = .data %>% attr("formula") %>% as.formula()

	.data %>%

		# Create plots for every tested transcript
		mutate(plot =
					 	pmap(
					 		list(
					 			`sample wise data`,
					 			!!as.symbol(.transcript),
					 			# nested data for plot
					 			.abundance,
					 			# name of value column
					 			.sample,
					 			# name of sample column
					 			parse_formula(formula)[1] # main covariate
					 		),
					 		~ produce_plots(..1, ..2, ..3, ..4, ..5)
					 	))
}


