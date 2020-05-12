#' Add attribute
#'
#' @param var A character
#' @param attribute An object
#' @param name A character
#' @export
add_attr = function(var, attribute, name){
	attr(var, name) <- attribute
	var
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @import dplyr
#' @importFrom purrr as_mapper
#'
#' @param .x A tibble
#' @param .p A boolean
#' @param .f1 A function
#' @param .f2 A function
#'
#'
#' @return A tibble
ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}

#' format_for_MPI
#'
#' @description Format reference data frame for MPI
#'
#' @param df A tibble
#' @param shards A integer
#' @param .sample A symbol
#'
format_for_MPI = function(df, shards, .sample) {
	.sample = enquo(.sample)

	df %>%

		left_join((.) %>%
								distinct(G) %>%
								arrange(G) %>%
								mutate(idx_MPI = head(
									rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
								)),
							by = "G") %>%
		arrange(idx_MPI, G) %>%

		# Decide start - end location
		group_by(idx_MPI) %>%
		do(
			(.) %>%
				left_join(
					(.) %>%
						distinct(!!.sample, G) %>%
						arrange(G) %>%
						count(G) %>%
						mutate(end = cumsum(n)) %>%
						mutate(start = c(
							1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)
						)),
					by = "G"
				)
		) %>%
		ungroup() %>%

		# Add symbol MPI rows indexes - otherwise spread below gives error
		left_join(
			(.) %>%
				group_by(idx_MPI) %>%
				distinct(G) %>%
				arrange(G) %>%
				mutate(`symbol MPI row` = 1:n()) %>%
				ungroup,
			by = c("G", "idx_MPI")
		) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(G) %>%
		mutate(`read count MPI row` = 1:n()) %>%
		ungroup

}

#' add_partition
#'
#' @description Add partition column dto data frame
#'
#' @param df.input A tibble
#' @param partition_by A symbol. Column we want to partition by
#' @param n_partitions An integer number of partition
add_partition = function(df.input, partition_by, n_partitions) {
	df.input %>%
		left_join(
			(.) %>%
				select(!!partition_by) %>%
				distinct %>%
				mutate(
					partition = 1:n() %>%
						divide_by(length((.))) %>%
						#	multiply_by(min(n_partitions, df.input %>% distinct(symbol) %>% nrow)) %>%
						multiply_by(n_partitions) %>%
						ceiling
				)
		)
}

#' Formula parser
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
	if (attr(terms(fm), "response") == 1)
		stop("The formula must be of the kind \"~ covariates\" ")
	else
		as.character(attr(terms(fm), "variables"))[-1]
}

#' Get matrix from tibble
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#'
#' @return A matrix
as_matrix <- function(tbl, rownames = NULL) {
	tbl %>%

		ifelse_pipe(
			tbl %>%
				ifelse_pipe(!is.null(rownames),		~ .x %>% dplyr::select(-contains(rownames))) %>%
				summarise_all(class) %>%
				gather(variable, class) %>%
				pull(class) %>%
				unique() %>%
				`%in%`(c("numeric", "integer")) %>% `!`() %>% any(),
			~ {
				warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
				.x
			}
		) %>%
		as.data.frame() %>%

		# Deal with rownames column if present
		ifelse_pipe(!is.null(rownames),
								~ .x  %>%
									set_rownames(tbl %>% pull(!!rownames)) %>%
									select(-!!rownames)) %>%

		# Convert to matrix
		as.matrix()
}

#' vb_iterative
#'
#' @description Runs iteratively variational bayes until it suceeds
#'
#' @importFrom rstan vb
#'
#' @param model A Stan model
#' @param output_samples An integer of how many samples from posteriors
#' @param iter An integer of how many max iterations
#' @param tol_rel_obj A real
#' @param additional_parameters_to_save A character vector
#' @param ... List of paramaters for vb function of Stan
#'
#' @return A Stan fit object
#'
vb_iterative = function(model,
												output_samples,
												iter,
												tol_rel_obj,
												additional_parameters_to_save,
												...) {
	res = NULL
	i = 0
	while (res %>% is.null | i > 5) {
		res = tryCatch({
			my_res = vb(
				model,
				output_samples = output_samples,
				iter = iter,
				tol_rel_obj = tol_rel_obj,
				#seed = 654321,
				pars=c("counts_rng", "exposure_rate", additional_parameters_to_save),
				...
			)
			boolFalse <- T
			return(my_res)
		},
		error = function(e) {
			i = i + 1
			writeLines(sprintf("Further attempt with Variational Bayes: %s", e))
			return(NULL)
		},
		finally = {
		})
	}

	return(res)
}

#' Choose the number of chains baed on how many draws we need from the posterior distribution
#' Because there is a fix cost (warmup) to starting a new chain,
#' we need to use the minimum amount that we can parallelise
#' @param how_many_posterior_draws A real number of posterior draws needed
#' @param max_number_to_check A sane upper plateau
#'
#' @return A Stan fit object
find_optimal_number_of_chains = function(how_many_posterior_draws,
																				 max_number_to_check = 100) {
	foreach(cc = 2:max_number_to_check, .combine = bind_rows) %do%
		{
			tibble(chains = cc, tot = how_many_posterior_draws / cc + 150 * cc)
		} %>%
		filter(tot == tot %>% min) %>%
		pull(chains)

}

#' Identify the optimal number of chain
#' based on how many draws we need from the posterior
#'
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map
#'
#' @param counts_MPI A matrix of read count information
#' @param to_exclude A vector of oulier data points to exclude
#' @param shards An integer
#'
#' @return A matrix
get_outlier_data_to_exlude = function(counts_MPI, to_exclude, shards) {
	# If there are genes to exclude
	switch(
		to_exclude %>% nrow %>% `>` (0) %>% `!` %>% sum(1),
		foreach(s = 1:shards, .combine = full_join) %do% {
			counts_MPI %>%
				inner_join(to_exclude, by = c("S", "G")) %>%
				filter(idx_MPI == s) %>%
				distinct(idx_MPI, `read count MPI row`) %>%
				rowid_to_column %>%
				spread(idx_MPI, `read count MPI row`) %>%

				# If a shard is empty create a dummy data set to avoid error
				ifelse_pipe((.) %>% nrow == 0, ~ tibble(rowid = 1,!!as.symbol(s) := NA))

		} %>%

			# Anonymous function - Add length array to the first row for indexing in MPI
			# Input: tibble
			# Output: tibble
			{
				bind_rows((.) %>% map(function(x)
					x %>% is.na %>% `!` %>% as.numeric %>% sum) %>% unlist,
					(.))
			} %>%

			select(-rowid) %>%
			replace(is.na(.), 0 %>% as.integer) %>%
			as_matrix() %>% t,

		# Otherwise
		matrix(rep(0, shards))
	)
}

#' function to pass initialisation values
#'
#' @return A list
inits_fx =
	function () {
		pars =
			res_discovery %>%
			filter(`.variable` != "counts_rng") %>%
			distinct(`.variable`) %>%
			pull(1)

		foreach(
			par = pars,
			.final = function(x)
				setNames(x, pars)
		) %do% {
			res_discovery %>%
				filter(`.variable` == par) %>%
				mutate(init = rnorm(n(), mean, sd)) %>%
				mutate(init = 0) %>%
				select(`.variable`, S, G, init) %>%
				pull(init)
		}
	}

#' Produce generated quantities plots with marked uotliers
#'
#' @importFrom purrr pmap
#' @importFrom purrr map_int
#' @import ggplot2
#'
#' @param .x A tibble
#' @param symbol A symbol object
#' @param .abundance A symbol object
#' @param .sample A symbol object
#' @param covariate A character string
#'
#' @return A ggplot
produce_plots = function(.x,
												 symbol,
												 .abundance,
												 .sample,
												 covariate) {
	# Set plot theme
	my_theme =
		theme_bw() +
		theme(
			panel.border = element_blank(),
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size = 12),
			#aspect.ratio = 1,
			axis.text.x = element_text(
				angle = 90,
				hjust = 1,
				vjust = 0.5
			),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(
				t = 10,
				r = 10,
				b = 10,
				l = 10
			)),
			axis.title.y  = element_text(margin = margin(
				t = 10,
				r = 10,
				b = 10,
				l = 10
			))
		)

	{
		ggplot(data = .x, aes(
			y = !!as.symbol(.abundance),
			x = !!as.symbol(.sample)
		)) +
			geom_errorbar(
				aes(ymin = `.lower`,
						ymax = `.upper`),
				width = 0,
				linetype = "dashed",
				color = "#D3D3D3"
			) +
			geom_errorbar(aes(
				ymin = `.lower_2`,
				ymax = `.upper_2`,
				color = `deleterious outliers`
			),
			width = 0) +
			scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
			my_theme
		} %>%

		ifelse_pipe(
			covariate %>% is.null %>% `!`,
			~ .x + geom_point(aes(
				size = `exposure rate`, fill = !!as.symbol(covariate)
			), shape = 21),
			~ .x + geom_point(
				aes(size = `exposure rate`),
				shape = 21,
				fill = "black"
			)
		) +
		ggtitle(symbol)
}

# Add annotation if sample belongs to high or low group
add_deleterious_if_covariate_exists = function(.data, X){
	.data %>%
		ifelse_pipe(
			X %>% ncol %>% `>` (1),
			~ .x %>%
				left_join(
					X %>%
						as_tibble %>%
						select(2) %>%
						setNames("factor or interest") %>%
						mutate(S = 1:n()) %>%
						mutate(`is group high` = `factor or interest` > mean(`factor or interest`)),
					by = "S"
				) %>%

				# Check if outlier might be deleterious for the statistics
				mutate(`deleterious outliers` = (!ppc) &
							 	(`is higher than mean` == `is group high`))
		)
}

#' merge_results
#'
#' @importFrom tidyr nest
#'
#' @param res_discovery A tibble
#' @param res_test A tibble
#' @param formula A formula
#' @param .sample A column name
#' @param .transcript A column name
#' @param .abundance A column name
#' @param do_check_only_on_detrimental A boolean
#'
#' @export
merge_results = function(res_discovery, res_test, formula, .transcript, .abundance, .sample, do_check_only_on_detrimental){

	res_discovery %>%
		filter(`.variable` == "counts_rng") %>%
		select(
			S,
			G,
			!!.transcript,
			!!.abundance,
			!!.sample,
			mean,
			`.lower`,
			`.upper`,
			`exposure rate`,
			one_of(parse_formula(formula))
		) %>%

		# Attach results of tests
		left_join(
			res_test %>% 		filter(`.variable` == "counts_rng") %>%
				select(
					S,
					G,
					mean,
					`.lower`,
					`.upper`,
					ppc,
					one_of(c("generated quantities", "deleterious outliers"))
				) %>%
				rename(mean_2 = mean, `.lower_2` = `.lower`, `.upper_2` = `.upper`),
			by = c("S", "G")
		) %>%

		# Check if new package is installed with different sintax
		ifelse_pipe(
			packageVersion("tidyr") >= "0.8.3.9000",
			~ .x %>% nest(`sample wise data` = c(-!!.transcript)),
			~ .x %>%
				group_by(!!.transcript) %>%
				nest(`sample wise data` = -!!.transcript)
		) %>%

		# Create plots for every tested transcript
		mutate(plot =
					 	pmap(
					 		list(
					 			`sample wise data`,
					 			!!.transcript,
					 			# nested data for plot
					 			quo_name(.abundance),
					 			# name of value column
					 			quo_name(.sample),
					 			# name of sample column
					 			parse_formula(formula)[1] # main covariate
					 		),
					 		~ produce_plots(..1, ..2, ..3, ..4, ..5)
					 	)) %>%

		# Add summary statistics
		mutate(`ppc samples failed` = map_int(`sample wise data`, ~ .x %>% pull(ppc) %>% `!` %>% sum)) %>%

		# If deleterious detection add summary as well
		ifelse_pipe(
			do_check_only_on_detrimental,
			~ .x %>%
				mutate(
					`tot deleterious outliers` =
						map_int(`sample wise data`, ~ .x %>% pull(`deleterious outliers`) %>% sum)
				)
		)
}

#' Select only significant genes plus background for efficient normalisation
#'
#' @importFrom rstan sampling
#' @importFrom rstan vb
#'
#' @param .data A tibble
#' @param .do_check A boolean
#' @param .significance A symbol
#' @param .transcript A column name
#' @param how_many_negative_controls An integer
#'
select_to_check_and_house_keeping = function(.data, .do_check, .significance, .transcript, how_many_negative_controls  = 500){
	.data %>%
		{
			bind_rows(
				# Genes to check
				(.) %>%
					filter((!!.do_check)),

				# Least changing genes, negative controls
				(.) %>%
					filter((!!.do_check) %>% `!`) %>%
					inner_join(
						(.) %>%
							arrange(!!.significance) %>%
							select(!!.transcript) %>%
							distinct() %>%
							tail(how_many_negative_controls),
						by = quo_name(.transcript)
					)
			)
		}
}


#' add_exposure_rate
#'
#' @importFrom tidyr separate
#'
#' @param .data A data frame
#' @param fit A fit object
#'
add_exposure_rate = function(.data, fit){

	writeLines(sprintf("executing %s", "add_exposure_rate"))

	.data %>%
		left_join(
			fit %>%
				summary("exposure_rate") %$%
				summary %>%
				as_tibble(rownames = ".variable") %>%
				separate(
					.variable,
					c(".variable", "S"),
					sep = "[\\[,\\]]",
					extra = "drop"
				) %>%
				mutate(S = S %>% as.integer) %>%
				rename(`exposure rate` = mean) %>%
				select(S, `exposure rate`),
			by = "S"
		)
}

check_if_within_posterior = function(.data, my_df, .do_check, .abundance){

	writeLines(sprintf("executing %s", "check_if_within_posterior"))

	.data %>%
		left_join(my_df, by = c("S", "G")) %>%
		filter((!!.do_check)) %>% # Filter only DE genes
		rowwise() %>%
		mutate(`ppc` = !!.abundance %>% between(`.lower`, `.upper`)) %>%
		mutate(`is higher than mean` = (!`ppc`) &
					 	(!!.abundance > mean)) %>%
		ungroup
}

#' fit_to_counts_rng
#'
#' @importFrom tidyr separate
#' @importFrom tidyr nest
#' @importFrom rstan summary
#'
#' @param fit A fit object
#' @param adj_prob_theshold fit real
#'
fit_to_counts_rng = function(fit, adj_prob_theshold){

	writeLines(sprintf("executing %s", "fit_to_counts_rng"))

	fit %>%
		rstan::summary("counts_rng",
									 prob = c(adj_prob_theshold, 1 - adj_prob_theshold)) %$%
		summary %>%
		as_tibble(rownames = ".variable") %>%
		separate(.variable,
						 c(".variable", "S", "G"),
						 sep = "[\\[,\\]]",
						 extra = "drop") %>%
		mutate(S = S %>% as.integer, G = G %>% as.integer) %>%
		select(-one_of(c("n_eff", "Rhat", "khat"))) %>%
		rename(`.lower` = (.) %>% ncol - 1,
					 `.upper` = (.) %>% ncol)
}

#' fit_to_counts_rng_approximated
#'
#' @importFrom tidyr separate
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom rstan summary
#' @importFrom furrr future_map
#' @importFrom future plan
#' @importFrom future multiprocess
#'
#' @param fit A fit object
#' @param adj_prob_theshold fit real
#' @param how_many_posterior_draws An integer
#' @param truncation_compensation A real
#' @param cores An integer
#'
#' @export

fit_to_counts_rng_approximated = function(fit, adj_prob_theshold, how_many_posterior_draws, truncation_compensation, cores){

	writeLines(sprintf("executing %s", "fit_to_counts_rng_approximated"))

	draws_mu =
		fit %>% extract("lambda_log_param") %>% `[[` (1) %>% as.data.frame() %>% setNames(sprintf("mu.%s", colnames(.))) %>%
		as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, mu, -.draw) %>% separate(par, c("par", "S", "G"), sep="\\.") %>% select(-par)
	draws_sigma =
		fit %>% extract("sigma_raw") %>% `[[` (1) %>% as.data.frame() %>% setNames(sprintf("sigma.%s", colnames(.) %>% gsub("V", "", .))) %>%
		as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, sigma, -.draw) %>% separate(par, c("par", "G"), sep="\\.") %>% select(-par)
	draws_exposure =
		fit %>% extract("exposure_rate") %>% `[[` (1) %>% as.data.frame() %>% setNames(sprintf("exposure.%s", colnames(.) %>% gsub("V", "", .))) %>%
		as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, exposure, -.draw) %>% separate(par, c("par", "S"), sep="\\.") %>% select(-par)

	draws_mu %>%
		left_join(draws_sigma) %>%
		left_join(draws_exposure) %>%
		nest(data = -c(S, G)) %>%
		mutate(
			CI = map(
				data,
				~ {
					.x_supersampled = .x %>%	sample_n(how_many_posterior_draws, replace = T)

					draws = rnbinom(n =how_many_posterior_draws,	mu = exp(.x_supersampled$mu + .x_supersampled$exposure),	size = 1/exp(.x_supersampled$sigma) * truncation_compensation	)

					draws %>%

						# Process quantile
						quantile(c(adj_prob_theshold, 1 - adj_prob_theshold)) %>%

					tibble::as_tibble(rownames="prop") %>%
						tidyr::spread(prop, value) %>%
						setNames(c(".lower", ".upper")) %>%

						# Add mean and sd
						dplyr::mutate(mean = mean(draws), sd = sd(draws))
				}
			)
		) %>%
		select(-data) %>%
		unnest(CI) %>%

		# Adapt to old dataset
		mutate(.variable = "counts_rng") %>%
		mutate(S = as.integer(S), G = as.integer(G))


}

save_generated_quantities_in_case = function(.data, fit, save_generated_quantities){

	writeLines(sprintf("executing %s", "save_generated_quantities_in_case"))

	.data %>%
		ifelse_pipe(
			save_generated_quantities,
			~ .x %>%

				# Add generated quantities
				left_join(fit %>% tidybayes::gather_draws(counts_rng[S, G])) %>%

				# Nest them in the data frame
				nest(`generated quantities` = c(.chain, .iteration, .draw, .value ))

		)
}

check_columns_exist = function(.data, .sample, .transcript, .abundance, .significance, .do_check){

	# Prepare column same enquo
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.significance = enquo(.significance)
	.do_check = enquo(.do_check)

	columns = c(quo_name(.sample), quo_name(.transcript), quo_name(.abundance), quo_name(.significance), quo_name(.do_check))
	if((!columns %in% (.data %>% colnames)) %>% any)
		stop(
			sprintf(
				"The columns %s are not present in your tibble",
				paste(columns[(!columns %in% (.data %>% colnames))], collapse=" ")
			)
		)
}


#' Check if NA
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr enquo
#'
#' @param .data A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param .sample A column name
#' @param .transcript A column name
#' @param .abundance A column name
#' @param .significance A column name
#' @param .do_check A column name
#' @param formula_columns A symbol vector
#'
check_if_any_NA = function(.data, .sample, .transcript, .abundance, .significance, .do_check, formula_columns){

	# Prepare column same enquo
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.significance = enquo(.significance)
	.do_check = enquo(.do_check)

	columns = c(quo_name(.sample), quo_name(.transcript), quo_name(.abundance), quo_name(.significance), quo_name(.do_check), formula_columns)

	if(
		.data %>%
		drop_na(columns) %>%
		nrow %>% `<`
		(
			.data %>% nrow
		)
	)
		stop(sprintf("There are NA values in you tibble for any of the column %s", paste(columns, collapse=", ")))
}

detect_cores = function(){

	if(.Platform$OS.type == "unix")
		system("nproc", intern = TRUE) %>% as.integer %>% sum(-1)
	else if(.Platform$OS.type == "windows")
		parallel::detectCores()  %>% as.integer %>% sum(-1)
	else stop("Your platform type is not recognised")

}

#' Create the design matrix
#'
#' @param .data A tibble
#' @param formula A formula
#' @param .sample A symbol
#' @export
create_design_matrix = function(.data, formula, .sample){

	.sample = enquo(.sample)

	model.matrix(
		object = formula,
		data =
			.data %>%
			select(!!.sample, one_of(parse_formula(formula))) %>%
			distinct %>% arrange(!!.sample)

	)

}

#' Format the input
#'
#' @param .data A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param .sample A column name
#' @param .transcript A column name
#' @param .abundance A column name
#' @param .do_check A symbol
#' @param .significance A column name
#' @param how_many_negative_controls An integer
#'
#' @export
format_input = function(.data, formula, .sample, .transcript, .abundance, .do_check, .significance, how_many_negative_controls = 500){

	# Prepare column same enquo
	.sample =       enquo(.sample)
	.transcript =         enquo(.transcript)
	.abundance =        enquo(.abundance)
	.do_check =     enquo(.do_check)
	.significance = enquo(.significance)

	.data %>%

		# Select only significant genes plus background for efficient normalisation
		select_to_check_and_house_keeping(.do_check, .significance, .transcript, how_many_negative_controls) %>%

		# Prepare the data frame
		select(
			!!.transcript,
			!!.sample,
			!!.abundance,
			one_of(parse_formula(formula)),
			!!.do_check
		) %>%
		distinct() %>%

		# Add symbol idx
		left_join((.) %>%
								distinct(!!.transcript) %>%
								mutate(G = 1:n()),
							by = quo_name(.transcript)) %>%

		# Add sample indeces
		mutate(S = factor(
			!!.sample,
			levels = (.) %>% pull(!!.sample) %>% unique
		) %>% as.integer)
}


run_model = function(model, approximate_posterior_inference, chains, how_many_posterior_draws, inits_fx, tol_rel_obj, additional_parameters_to_save, seed){


	writeLines(sprintf("executing %s", "run_model"))

	switch(
		approximate_posterior_inference %>% `!` %>% as.integer %>% sum(1),

		# VB Repeat strategy for failures of vb
		vb_iterative(
			model,
			#pcc_seq_model, #
			output_samples = how_many_posterior_draws,
			iter = 50000,
			tol_rel_obj = tol_rel_obj,
			pars = c(
				"counts_rng",
				"exposure_rate",
				additional_parameters_to_save
			)
			#,
			#sample_file = "temp_stan_sampling.txt"
		),

		# MCMC
		sampling(
			model,
			#pcc_seq_model, #
			chains = chains,
			cores = chains,
			iter = (how_many_posterior_draws / chains) %>% ceiling %>% sum(150),
			warmup = 300,
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

}
