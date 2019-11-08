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
#' @param sample_column A symbol
#'
format_for_MPI = function(df, shards, sample_column) {
	sample_column = enquo(sample_column)

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
						distinct(!!sample_column, G) %>%
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
				seed = 654321,
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
#' @param value_column A symbol object
#' @param sample_column A symbol object
#' @param covariate A character string
#'
#' @return A ggplot
produce_plots = function(.x,
												 symbol,
												 value_column,
												 sample_column,
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
			aspect.ratio = 1,
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
			y = !!as.symbol(value_column),
			x = !!as.symbol(sample_column)
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
add_deleterious_if_covariate_exists = function(input.df, X){
	input.df %>%
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
#' @param sample_column A column name
#' @param gene_column A column name
#' @param value_column A column name
#' @param do_check_only_on_detrimental A boolean
#'
#' @export
merge_results = function(res_discovery, res_test, formula, gene_column, value_column, sample_column, do_check_only_on_detrimental){

	res_discovery %>%
		filter(`.variable` == "counts_rng") %>%
		select(
			S,
			G,
			!!gene_column,
			!!value_column,
			!!sample_column,
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
			packageVersion("tidyr") == "0.8.3.9000",
			~ .x %>% nest(`sample wise data` = c(-!!gene_column)),
			~ .x %>%
				group_by(!!gene_column) %>%
				nest(-!!gene_column, .key = `sample wise data`)
		) %>%

		# Create plots for every tested transcript
		mutate(plot =
					 	pmap(
					 		list(
					 			`sample wise data`,
					 			symbol,
					 			# nested data for plot
					 			quo_name(value_column),
					 			# name of value column
					 			quo_name(sample_column),
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
#' @param input.df A tibble
#' @param do_check_column A boolean
#' @param significance_column A symbol
#' @param gene_column A column name
#' @param how_many_negative_controls An integer
#'
select_to_check_and_house_keeping = function(input.df, do_check_column, significance_column, gene_column, how_many_negative_controls  = 500){
	input.df %>%
		{
			bind_rows(
				# Genes to check
				(.) %>%
					filter((!!do_check_column)),

				# Least changing genes, negative controls
				(.) %>%
					filter((!!do_check_column) %>% `!`) %>%
					inner_join(
						(.) %>%
							arrange(!!significance_column) %>%
							select(!!gene_column) %>%
							distinct() %>%
							tail(how_many_negative_controls),
						by = quo_name(gene_column)
					)
			)
		}
}

run_model = function(model, approximate_posterior_inference, chains, how_many_posterior_draws, inits_fx, tol_rel_obj, additional_parameters_to_save){

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
			warmup = 150,
			save_warmup = FALSE,
			seed = 654321,
			init = inits_fx,
			pars = c(
				"counts_rng",
				"exposure_rate",
				additional_parameters_to_save
			)
		)
	)
}

#' add_exposure_rate
#'
#' @importFrom tidyr separate
#'
#' @param input.df A data frame
#' @param fit A fit object
#'
add_exposure_rate = function(input.df, fit){

	writeLines(sprintf("executing %s", "add_exposure_rate"))

	input.df %>%
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

check_if_within_posterior = function(input.df, my_df, do_check_column, value_column){

	writeLines(sprintf("executing %s", "check_if_within_posterior"))

	input.df %>%
		left_join(my_df, by = c("S", "G")) %>%
		filter((!!do_check_column)) %>% # Filter only DE genes
		rowwise() %>%
		mutate(`ppc` = !!value_column %>% between(`.lower`, `.upper`)) %>%
		mutate(`is higher than mean` = (!`ppc`) &
					 	(!!value_column > mean)) %>%
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
#' @param cores An integer
#'
#' @export

fit_to_counts_rng_approximated = function(fit, adj_prob_theshold, how_many_posterior_draws, cores){

	writeLines(sprintf("executing %s", "fit_to_counts_rng_approximated"))

	fit_summary  = fit %>% rstan::summary() %$% summary %>%
		as_tibble(rownames="par")

	fit_mu =
		fit_summary %>%
		filter(grepl("lambda_log_param", par)) %>%
		rename(mu_mean = mean, mu_sd = sd) %>%
		separate(par, c(".variable", "S", "G"), sep="[\\[\\]\\,]", extra = "drop") %>%
		mutate(S = S %>% as.integer, G = G %>% as.integer)

	fit_sigma =
		fit_summary %>%
		filter(grepl("sigma_raw", par))	%>%
		rename(sigma_mean = mean, sigma_sd = sd) %>%
		separate(par, c(".variable", "G"), sep="[\\[\\]\\,]", extra = "drop") %>%
		mutate( G = G %>% as.integer)

	# Calculate quantiles
	fit_mu %>%
		select(S, G, mu_mean, mu_sd) %>%
		left_join(
			fit_sigma %>%
				select(G, sigma_mean, sigma_sd),
			by="G"
		) %>%

		# Pass variables
		mutate(how_many_posterior_draws = how_many_posterior_draws, adj_prob_theshold = adj_prob_theshold) %>%

		do_parallel_start(cores, "G") %>%
		do({

			`%>%` = magrittr::`%>%`

			dplyr::do(
				dplyr::group_by((.), S, G),
				{
					.x = (.)
					how_many_posterior_draws = unique(.x$how_many_posterior_draws)
					adj_prob_theshold = unique(.x$adj_prob_theshold)

					x = rnbinom(
						how_many_posterior_draws,
						mu = exp(rnorm(how_many_posterior_draws, .x$mu_mean, .x$mu_sd)),
						size = 1/exp(rnorm(how_many_posterior_draws, .x$sigma_mean, .x$sigma_sd))
					)

					x %>%
						quantile(c(adj_prob_theshold, 1 - adj_prob_theshold)) %>%
						tibble::as_tibble(rownames="prop") %>%
						tidyr::spread(prop, value) %>%
						setNames(c(".lower", ".upper")) %>%

						# Add mean and sd
						dplyr::mutate(mean = mean(x), sd = sd(x))
				})
		}) %>%
		do_parallel_end() %>%


		mutate(.variable = "counts_rng") %>%
		select(.variable, S, G, mean, sd, .lower, .upper)

}

save_generated_quantities_in_case = function(input.df, fit, save_generated_quantities){

	writeLines(sprintf("executing %s", "save_generated_quantities_in_case"))

	input.df %>%
		ifelse_pipe(
			save_generated_quantities,
			~ .x %>%

				# Add generated quantities
				left_join(fit %>% tidybayes::gather_draws(counts_rng[S, G])) %>%

				# Nest them in the data frame
				group_by(`.variable`, S , G, mean, se_mean , sd, `.lower`, `.upper`) %>%
				nest(.key = "generated quantities")
		)
}

check_columns_exist = function(input.df, sample_column, gene_column, value_column, significance_column, do_check_column){

	# Prepare column same enquo
	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)

	columns = c(quo_name(sample_column), quo_name(gene_column), quo_name(value_column), quo_name(significance_column), quo_name(do_check_column))
	if((!columns %in% (input.df %>% colnames)) %>% any)
		stop(
			sprintf(
				"The columns %s are not present in your tibble",
				paste(columns[(!columns %in% (input.df %>% colnames))], collapse=" ")
			)
		)
}


#' Check if NA
#'
#' @importFrom tidyr drop_na
#' @importFrom dplyr enquo
#'
#' @param input.df A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param sample_column A column name
#' @param gene_column A column name
#' @param value_column A column name
#' @param significance_column A column name
#' @param do_check_column A column name
#' @param formula_columns A symbol vector
#'
check_if_any_NA = function(input.df, sample_column, gene_column, value_column, significance_column, do_check_column, formula_columns){

	# Prepare column same enquo
	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)

	columns = c(quo_name(sample_column), quo_name(gene_column), quo_name(value_column), quo_name(significance_column), quo_name(do_check_column), formula_columns)

	if(
		input.df %>%
		drop_na(columns) %>%
		nrow %>% `<`
		(
			input.df %>% nrow
		)
	)
		stop(sprintf("There are NA values in you tibble for any of the column %s", paste(columns, collapse=", ")))
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
#' @importFrom ttBulk scale_abundance
#'
#' @param my_df A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param sample_column A column name
#' @param gene_column A column name
#' @param value_column A column name
#' @param significance_column A column name
#' @param do_check_column A column name
#' @param approximate_posterior_inference A boolean
#' @param approximate_posterior_analysis A boolean
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
#'
#' @return A tibble with additional columns
#'
#' @export
do_inference = function(my_df,
												formula,
												sample_column ,
												gene_column ,
												value_column,
												significance_column ,
												do_check_column,
												approximate_posterior_inference = F,
												approximate_posterior_analysis = F,
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
												write_on_disk = F) {

	writeLines(sprintf("executing %s", "do_inference"))

	# Prepare column names
	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)

	# Get the number of transcripts to check
	how_many_to_check =
		my_df %>%
		filter(!!do_check_column) %>%
		distinct(!!gene_column) %>%
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
		select(!!gene_column,!!sample_column,!!value_column, S, G) %>%
		format_for_MPI(shards,!!sample_column)

	# Setup dimensions of variables for the model
	G = counts_MPI %>% distinct(G) %>% nrow()
	S = counts_MPI %>% distinct(!!sample_column) %>% nrow()
	N = counts_MPI %>% distinct(idx_MPI,!!value_column, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_MPI %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G_per_shard = counts_MPI %>% distinct(!!gene_column, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array
	n_shards = min(shards, counts_MPI %>% distinct(idx_MPI) %>% nrow)
	G_per_shard_idx = c(
		0,
		counts_MPI %>% distinct(!!gene_column, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum
	)

	# Read count object
	counts =
		counts_MPI %>%
		distinct(idx_MPI,!!value_column, `read count MPI row`)  %>%
		spread(idx_MPI,!!value_column) %>%
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
				seed = 654321,
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
			~ .x %>% fit_to_counts_rng_approximated(adj_prob_theshold, how_many_posterior_draws * 10, cores),
			~ .x %>% fit_to_counts_rng(adj_prob_theshold)
		) %>%

		# If generated quantities are saved
		save_generated_quantities_in_case(fit, save_generated_quantities) %>%

		# Add exposure rate
		#add_exposure_rate(fit) %>%

		# Check if data is within posterior
		check_if_within_posterior(my_df, do_check_column, value_column) %>%

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

detect_cores = function(){

	if(.Platform$OS.type == "unix")
		system("nproc", intern = TRUE) %>% as.integer %>% sum(-1)
	else if(.Platform$OS.type == "windows")
		parallel::detectCores()  %>% as.integer %>% sum(-1)
	else stop("Your platform type is not recognised")

}

#' Create the design matrix
#'
#' @param input.df A tibble
#' @param formula A formula
#' @param sample_column A symbol
#' @export
create_design_matrix = function(input.df, formula, sample_column){

	sample_column = enquo(sample_column)

	model.matrix(
		object = formula,
		data =
			input.df %>%
			select(!!sample_column, one_of(parse_formula(formula))) %>%
			distinct %>% arrange(!!sample_column)

	)

}

#' Format the input
#'
#' @param input.df A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param sample_column A column name
#' @param gene_column A column name
#' @param value_column A column name
#' @param do_check_column A symbol
#' @param significance_column A column name
#' @param how_many_negative_controls An integer
#'
#' @export
format_input = function(input.df, formula, sample_column, gene_column, value_column, do_check_column, significance_column, how_many_negative_controls = 500){

	# Prepare column same enquo
	sample_column =       enquo(sample_column)
	gene_column =         enquo(gene_column)
	value_column =        enquo(value_column)
	do_check_column =     enquo(do_check_column)
	significance_column = enquo(significance_column)

	input.df %>%

		# Select only significant genes plus background for efficient normalisation
		select_to_check_and_house_keeping(do_check_column, significance_column, gene_column, how_many_negative_controls) %>%

		# Prepare the data frame
		select(
			!!gene_column,
			!!sample_column,
			!!value_column,
			one_of(parse_formula(formula)),
			!!do_check_column
		) %>%
		distinct() %>%

		# Add symbol idx
		left_join((.) %>%
								distinct(!!gene_column) %>%
								mutate(G = 1:n()),
							by = quo_name(gene_column)) %>%

		# Add sample indeces
		mutate(S = factor(
			!!sample_column,
			levels = (.) %>% pull(!!sample_column) %>% unique
		) %>% as.integer)
}

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
#' @importFrom ttBulk scale_abundance
#' @importFrom benchmarkme get_ram
#' @importFrom magrittr multiply_by
#'
#' @param input.df A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param sample_column A column name
#' @param gene_column A column name
#' @param value_column A column name
#' @param significance_column A column name
#' @param approximate_posterior_inference A boolean
#' @param approximate_posterior_analysis A boolean
#' @param do_check_column A symbol
#' @param how_many_negative_controls An integer
#' @param save_generated_quantities A boolean
#' @param additional_parameters_to_save A character vector
#' @param cores An integer
#' @param percent_false_positive_genes A real
#' @param pass_fit A boolean
#' @param do_check_only_on_detrimental A boolean
#' @param tol_rel_obj A real
#' @param just_discovery A boolean
#' @param write_on_disk A boolean
#'
#' @return A tibble with additional columns
#'
#' @export
#'
ppc_seq = function(input.df,
									 formula = ~ 1,
									 sample_column = `sample`,
									 gene_column = `symbol`,
									 value_column = `read count`,
									 significance_column = `p-value`,
									 do_check_column,
									 approximate_posterior_inference = T,
									 approximate_posterior_analysis = F,
									 how_many_negative_controls = 500,
									 save_generated_quantities = F,
									 additional_parameters_to_save = c(),  # For development purpose
									 cores = detect_cores(), # For development purpose,
									 percent_false_positive_genes = "1%",
									 pass_fit = F,
									 do_check_only_on_detrimental = length(parse_formula(formula)) > 0,
									 tol_rel_obj = 0.01,
									 just_discovery = F,
									 write_on_disk = F
) {
	# Prepare column same enquo
	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)

	# Get factor of interest
	#factor_of_interest = ifelse(parse_formula(formula) %>% length %>% `>` (0), parse_formula(formula)[1], "")

	# Check if columns exist
	check_columns_exist(input.df, !!sample_column, !!gene_column, !!value_column, !!significance_column, !!do_check_column)

	# Check if any column is NA or null
	check_if_any_NA(input.df, !!sample_column, !!gene_column, !!value_column, !!significance_column, !!do_check_column, parse_formula(formula))

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
	if (input.df %>% filter(!!gene_column %>% is.na) %>% nrow > 0)
		stop("There are NAs in the gene_column. Please filter those records")

	# Check if the counts column is an integer
	if (input.df %>% select(!!value_column) %>% sapply(class) != "integer")
		stop(
			sprintf(
				"The column %s must be of class integer. You can do as mutate(`%s` = `%s` %%>%% as.integer)",
				quo_name(value_column),
				quo_name(value_column),
				quo_name(value_column)
			)
		)

	# Calculate the adj_prob_theshold
	adj_prob_theshold_1  = 0.05
	adj_prob_theshold_2 =
		pfpg / 100 /
		(input.df %>% distinct(!!sample_column) %>% nrow) *
		ifelse(do_check_only_on_detrimental, 2, 1)

	# Calculate adj_prob_theshold
	how_many_posterior_draws_1 =  5 %>% divide_by(adj_prob_theshold_1) %>% max(1000) # I want 5 draws in the tail
	how_many_posterior_draws_2 =  5 %>% divide_by(adj_prob_theshold_2) %>% max(1000) # I want 5 draws in the tail

	# Check if enough memory for full draw
	available_memory = get_ram() %>% as.numeric() %>% multiply_by(1e+9)
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
		input.df,
		formula,
		!!sample_column,
		!!gene_column,
		!!value_column,
		!!do_check_column,
		!!significance_column,
		how_many_negative_controls
	)

	# Create design matrix
	X = create_design_matrix(my_df, formula, !!sample_column)

	C = X %>% ncol

	# Prior info
	lambda_mu_mu = 5.612671

	# Build better scales for the inference
	exposure_rate_multiplier =
		my_df %>%
		scale_abundance(!!sample_column,!!gene_column,!!value_column) %>%
		distinct(!!sample_column, TMM, multiplier) %>%
		mutate(l = multiplier %>% log) %>%
		summarise(l %>% sd) %>%
		pull(`l %>% sd`)

	# Build better scales for the inference
	intercept_shift_scale =
		my_df %>%
		scale_abundance(!!sample_column,!!gene_column,!!value_column) %>%
		mutate(cc =
					 	!!as.symbol(sprintf(
					 		"%s normalised",  quo_name(value_column)
					 	)) %>%
					 	`+` (1) %>% log) %>%
		summarise(shift = cc %>% mean, scale = cc %>% sd) %>%
		as.numeric

	# Run the first discovery phase with permissive false discovery rate
	res_discovery =
		my_df %>%
		do_inference(
			formula,!!sample_column ,!!gene_column ,!!value_column ,!!significance_column ,!!do_check_column,
			approximate_posterior_inference,
			approximate_posterior_analysis = F,
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
			write_on_disk = write_on_disk
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
			formula,!!sample_column ,!!gene_column ,!!value_column ,!!significance_column ,!!do_check_column,
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
			write_on_disk = write_on_disk
		)

	# Merge results and return
	merge_results(res_discovery, res_test, formula, gene_column, value_column, sample_column, do_check_only_on_detrimental) %>%

		# Add fit attribute if any
		add_attr(res_discovery %>% attr("fit"), "fit 1") %>%
		add_attr(res_test %>% attr("fit"), "fit 2") %>%

		# Add total draws
		add_attr(res_test %>% attr("total_draws"), "total_draws")

}
