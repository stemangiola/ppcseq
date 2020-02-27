
#' pcc_seq main
#'
#' @description This function calls the stan model.
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
#' @importFrom tidybulk scale_abundance
#' @importFrom benchmarkme get_ram
#' @importFrom magrittr multiply_by
#'
#' @param .data A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param .sample A column name
#' @param .transcript A column name
#' @param .abundance A column name
#' @param .significance A column name
#' @param approximate_posterior_inference A boolean
#' @param approximate_posterior_analysis A boolean
#' @param do_correct_approx A boolean
#' @param .do_check A symbol
#' @param how_many_negative_controls An integer
#' @param draws_after_tail An integer. How many draws should on average be after the tail, in a way to inform CI
#' @param save_generated_quantities A boolean
#' @param additional_parameters_to_save A character vector
#' @param cores An integer
#' @param percent_false_positive_genes A real
#' @param pass_fit A boolean
#' @param do_check_only_on_detrimental A boolean
#' @param tol_rel_obj A real
#' @param just_discovery A boolean
#' @param seed an integer
#' @param adj_prob_theshold_2 A boolean. Used for development and testing purposes
#'
#' @return A tibble with additional columns
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
														 approximate_posterior_inference = T,
														 approximate_posterior_analysis = NULL,
														 do_correct_approx = T,
														 how_many_negative_controls = 500,
														 draws_after_tail = 10,
														 save_generated_quantities = F,
														 additional_parameters_to_save = c(),  # For development purpose
														 cores = detect_cores(), # For development purpose,
														 percent_false_positive_genes = "1%",
														 pass_fit = F,
														 do_check_only_on_detrimental = length(parse_formula(formula)) > 0,
														 tol_rel_obj = 0.01,
														 just_discovery = F,
														 seed = sample(1:99999, size = 1),
														 adj_prob_theshold_2 = NULL
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

	# Check is testing environment is supported
	if (approximate_posterior_inference &	save_generated_quantities)
		stop("Variational Bayes does not support tidybayes needed for save_generated_quantities, use sampling")

	# Check percent FP input
	pfpg = percent_false_positive_genes %>% gsub("%$", "", .) %>% as.numeric
	if (pfpg %>% is.na |
			!(pfpg %>% between(0, 100)))
		stop("percent_false_positive_genes must be a string from > 0% to < 100%")

	# For reference MPI inference
	# Check if all trannscripts are non NA
	if (.data %>% filter(!!.transcript %>% is.na) %>% nrow > 0)
		stop("There are NAs in the .transcript. Please filter those records")

	# Check if the counts column is an integer
	if (.data %>% select(!!.abundance) %>% sapply(class) != "integer")
		stop(
			sprintf(
				"The column %s must be of class integer. You can do as mutate(`%s` = `%s` %%>%% as.integer)",
				quo_name(.abundance),
				quo_name(.abundance),
				quo_name(.abundance)
			)
		)

	# Calculate the adj_prob_theshold
	adj_prob_theshold_1  = 0.05
	if(adj_prob_theshold_2 %>% is.null)
		adj_prob_theshold_2 =
		pfpg / 100 /
		(.data %>% distinct(!!.sample) %>% nrow) *
		ifelse(do_check_only_on_detrimental, 2, 1)

	print(sprintf("adj_prob_theshold_2 = %s", adj_prob_theshold_2))

	# Calculate adj_prob_theshold
	how_many_posterior_draws_1 =  draws_after_tail %>% divide_by(adj_prob_theshold_1) %>% max(1000) # I want 5 draws in the tail
	how_many_posterior_draws_2 =  draws_after_tail %>% divide_by(adj_prob_theshold_2) %>% max(1000) # I want 5 draws in the tail

	print(sprintf("how_many_posterior_draws_2 = %s", how_many_posterior_draws_2))

	# If too many draws required revert to approximation of CI
	if(approximate_posterior_analysis %>% is.null){
		if(how_many_posterior_draws_2 > 20000) {
			writeLines(sprintf("The number of draws needed to calculate the CI from the posterior would be larger than %s. To avoid impractical computation times, the calculation of the CI will be based on the mean, exposure and overdisperison posteriors.", how_many_posterior_draws_2))
			approximate_posterior_analysis = T
		} else approximate_posterior_analysis = F
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
		approximate_posterior_analysis = T
	}


	# distinct_at is not released yet for dplyr, thus we have to use this trick
	my_df <- format_input(
		.data,
		formula,
		!!.sample,
		!!.transcript,
		!!.abundance,
		!!.do_check,
		!!.significance,
		how_many_negative_controls
	)

	# Create design matrix
	X = create_design_matrix(my_df, formula, !!.sample)

	C = X %>% ncol

	# Prior info
	lambda_mu_mu = 5.612671

	# Build better scales for the inference
	exposure_rate_multiplier =
		my_df %>%
		scale_abundance(!!.sample,!!.transcript,!!.abundance) %>%
		distinct(!!.sample, TMM, multiplier) %>%
		mutate(l = multiplier %>% log) %>%
		summarise(l %>% sd) %>%
		pull(`l %>% sd`)

	# Build better scales for the inference
	intercept_shift_scale =
		my_df %>%
		scale_abundance(!!.sample,!!.transcript,!!.abundance) %>%
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
			approximate_posterior_analysis = F,
			do_correct_approx = do_correct_approx,
			C,
			X,
			lambda_mu_mu,
			cores,
			exposure_rate_multiplier,
			intercept_shift_scale,
			additional_parameters_to_save,
			adj_prob_theshold  = adj_prob_theshold_1,
			how_many_posterior_draws = how_many_posterior_draws_1,
			pass_fit = pass_fit,
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
			~ .x %>% filter(ppc)
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
			do_correct_approx = do_correct_approx,
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
	merge_results(res_discovery, res_test, formula, .transcript, .abundance, .sample, do_check_only_on_detrimental) %>%

		# Add fit attribute if any
		add_attr(res_discovery %>% attr("fit"), "fit 1") %>%
		add_attr(res_test %>% attr("fit"), "fit 2") %>%

		# Add total draws
		add_attr(res_test %>% attr("total_draws"), "total_draws")

}


#' do_inference
#'
#' @description This function calls the stan model.
#'
#' @importFrom tibble tibble
#' @import rstan
#' @import dplyr
#' @importFrom tidyr spread
#' @import tidybayes
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#' @importFrom magrittr %$%
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom tidybulk scale_abundance
#'
#' @param my_df A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param .sample A column name
#' @param .transcript A column name
#' @param .abundance A column name
#' @param .significance A column name
#' @param .do_check A column name
#' @param approximate_posterior_inference A boolean
#' @param approximate_posterior_analysis A boolean
#' @param do_correct_approx A boolean
#' @param C An integer
#' @param X A tibble
#' @param lambda_mu_mu A real
#' @param cores An integer
#' @param exposure_rate_multiplier A real
#' @param intercept_shift_scale A real
#' @param additional_parameters_to_save A character vector
#' @param adj_prob_theshold A real
#' @param how_many_posterior_draws A real number of posterior draws needed
#' @param to_exclude A boolean
#' @param truncation_compensation A real
#' @param save_generated_quantities A boolean
#' @param inits_fx A function
#' @param prior_from_discovery A tibble
#' @param pass_fit A fit
#' @param tol_rel_obj A real
#' @param write_on_disk A boolean
#' @param seed an integer
#'
#' @return A tibble with additional columns
#'
do_inference = function(my_df,
												formula,
												.sample ,
												.transcript ,
												.abundance,
												.significance ,
												.do_check,
												approximate_posterior_inference = F,
												approximate_posterior_analysis = F,
												do_correct_approx = T,
												C,
												X,
												lambda_mu_mu,
												cores,
												exposure_rate_multiplier,
												intercept_shift_scale,
												additional_parameters_to_save,
												adj_prob_theshold,
												how_many_posterior_draws,
												to_exclude = tibble(S = integer(), G = integer()),
												truncation_compensation = 1,
												save_generated_quantities = F,
												inits_fx = "random",
												prior_from_discovery = tibble(`.variable` = character(),
																											mean = numeric(),
																											sd = numeric()),
												pass_fit = F,
												tol_rel_obj = 0.01,
												write_on_disk = F,
												seed) {

	writeLines(sprintf("executing %s", "do_inference"))

	# Prepare column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.significance = enquo(.significance)
	.do_check = enquo(.do_check)

	# Check that the dataset is squared
	if(my_df %>% distinct(!!.sample, !!.transcript) %>% count(!!.transcript) %>% count(n) %>% nrow %>% `>` (1))
		stop("The input data frame does not represent a rectangular structure. Each transcript must be present in all samples.")

	# Get the number of transcripts to check
	how_many_to_check =
		my_df %>%
		filter(!!.do_check) %>%
		distinct(!!.transcript) %>%
		nrow

	# if analysis approximated
	# If posterior analysis is approximated I just need enough
	how_many_posterior_draws_practical = ifelse(approximate_posterior_analysis, 1000, how_many_posterior_draws)
	if(approximate_posterior_analysis) additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique

	# Identify the optimal number of chain
	# based on how many draws we need from the posterior
	chains =
		find_optimal_number_of_chains(how_many_posterior_draws_practical) %>%
		min(cores) %>%
		max(3)

	# Find how many cores per chain, minimum 1 of course
	my_cores = cores %>% divide_by(chains) %>% floor %>% max(1)

	# Set the number of data partition = to the number of cores per chain
	shards = my_cores

	# Setup the data shards
	counts_MPI =
		my_df %>%
		select(!!.transcript,!!.sample,!!.abundance, S, G) %>%
		format_for_MPI(shards,!!.sample)

	# Setup dimensions of variables for the model
	G = counts_MPI %>% distinct(G) %>% nrow()
	S = counts_MPI %>% distinct(!!.sample) %>% nrow()
	N = counts_MPI %>% distinct(idx_MPI,!!.abundance, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_MPI %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G_per_shard = counts_MPI %>% distinct(!!.transcript, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array
	n_shards = min(shards, counts_MPI %>% distinct(idx_MPI) %>% nrow)
	G_per_shard_idx = c(
		0,
		counts_MPI %>% distinct(!!.transcript, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum
	)

	# Read count object
	counts =
		counts_MPI %>%
		distinct(idx_MPI,!!.abundance, `read count MPI row`)  %>%
		spread(idx_MPI,!!.abundance) %>%
		select(-`read count MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	# Indexes of the samples
	sample_idx =
		counts_MPI %>%
		distinct(idx_MPI, S, `read count MPI row`)  %>%
		spread(idx_MPI, S) %>%
		select(-`read count MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	# In the data structure when data for a transcript starts and ends
	symbol_end =
		counts_MPI %>%
		distinct(idx_MPI, end, `symbol MPI row`)  %>%
		spread(idx_MPI, end) %>%
		bind_rows((.) %>% head(n = 1) %>%  mutate_all(function(x) {
			0
		})) %>%
		arrange(`symbol MPI row`) %>%
		select(-`symbol MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	# Indexes of the transcripts
	G_ind =
		counts_MPI %>%
		distinct(idx_MPI, G, `symbol MPI row`)  %>%
		spread(idx_MPI, G) %>%
		arrange(`symbol MPI row`) %>%
		select(-`symbol MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	# Get the matrix of the idexes of the outlier data points
	# to explude from the model if it is the second passage
	to_exclude_MPI = get_outlier_data_to_exlude(counts_MPI, to_exclude, shards)

	# Package data
	counts_package =

		# Dimensions data sets
		rep(c(M, N, S), shards) %>%
		matrix(nrow = shards, byrow = T) %>%
		cbind(G_per_shard) %>%
		cbind(symbol_end) %>%
		cbind(sample_idx) %>%
		cbind(counts) %>%
		cbind(to_exclude_MPI)

	# Dimension od the data package to pass to Stan
	CP = ncol(counts_package)

	# Run model
	#writeLines(sprintf("- Roughly the memory allocation for the fit object is %s Gb", object.size(1:(S * how_many_to_check * how_many_posterior_draws))/1e9))

	# Set up environmental variable for threading
	Sys.setenv("STAN_NUM_THREADS" = my_cores)

	# Run Stan
	fit =
		switch(
			approximate_posterior_inference %>% `!` %>% as.integer %>% sum(1),

			# VB Repeat strategy for failures of vb
			vb_iterative(
				stanmodels$negBinomial_MPI,
				#pcc_seq_model, #
				output_samples = how_many_posterior_draws_practical,
				iter = 50000,
				tol_rel_obj = 0.005,
				additional_parameters_to_save = additional_parameters_to_save
			),

			# MCMC
			sampling(
				stanmodels$negBinomial_MPI,
				#pcc_seq_model, #
				chains = chains,
				cores = cores,
				iter = (how_many_posterior_draws_practical / chains) %>% ceiling %>% sum(150),
				warmup = 150,
				save_warmup = FALSE,
				seed = seed,
				init = inits_fx,
				pars = c(
					"counts_rng",
					"exposure_rate",
					additional_parameters_to_save
				)
			)
		)

	# Parse and return
	fit %>%

		ifelse_pipe(
			approximate_posterior_analysis,
			~ .x %>% fit_to_counts_rng_approximated(adj_prob_theshold, how_many_posterior_draws * 10, truncation_compensation, do_correct_approx, cores),
			~ .x %>% fit_to_counts_rng(adj_prob_theshold)
		) %>%

		# If generated quantities are saved
		save_generated_quantities_in_case(fit, save_generated_quantities) %>%

		# Add exposure rate
		#add_exposure_rate(fit) %>%

		# Check if data is within posterior
		check_if_within_posterior(my_df, .do_check, .abundance) %>%

		# Add annotation if sample belongs to high or low group
		add_deleterious_if_covariate_exists(X) %>%

		# Add position in MPI package for next inference
		left_join(counts_MPI %>% distinct(S, G, idx_MPI, `read count MPI row`),
							by = c("S", "G")) %>%

		# Add exposure rate
		add_exposure_rate(fit) %>%

		# needed for the figure article
		ifelse_pipe(pass_fit,	~ .x %>% add_attr(fit, "fit")	) %>%

		# Passing the amout of sampled data
		add_attr(S * how_many_to_check * how_many_posterior_draws, "total_draws")
}


