# Greater than
gt = function(a, b){	a > b }

# Smaller than
st = function(a, b){	a < b }

# Negation
not = function(is){	!is }

#' Get matrix from tibble
#'
#' @keywords internal
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#' @importFrom rlang quo_is_null
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @param do_check A boolean
#'
#' @return A matrix
#'
#' @noRd
as_matrix <- function(tbl,
											rownames = NULL,
											do_check = TRUE) {
	rownames = enquo(rownames)
	tbl %>%

		# Through warning if data frame is not numerical beside the rownames column (if present)
		when(
			do_check &&
				(.) %>%
				# If rownames defined eliminate it from the data frame
				when(!quo_is_null(rownames) ~ (.)[,-1], ~ (.)) %>%
				dplyr::summarise_all(class) %>%
				tidyr::gather(variable, class) %>%
				pull(class) %>%
				unique() %>%
				`%in%`(c("numeric", "integer")) %>% not() %>% any() ~ {
					warning("ppcseq says: there are NON-numerical columns, the matrix will NOT be numerical")
					(.)
				},
			~ (.)
		) %>%
		as.data.frame() %>%

		# Deal with rownames column if present
		when(
			!quo_is_null(rownames) ~ (.) %>%
				magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
				select(-1),
			~ (.)
		) %>%

		# Convert to matrix
		as.matrix()
}

#' Add attribute
#'
#' @keywords internal
#'
#' @param var A character
#' @param attribute An object
#' @param name A character
#'
#' @examples
#'
#' # Not needed because internal
#'
#' TRUE
#'
#' @return Same object
#'
#' @noRd
add_attr = function(var, attribute, name){
	attr(var, name) <- attribute
	var
}

#' This is a generalisation of ifelse that acceots an object and return an objects
#'
#' @keywords internal
#'
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
#' @noRd
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
#' @keywords internal
#'
#' @importFrom utils head
#'
#' @description Format reference data frame for MPI
#'
#' @param df A tibble
#' @param shards A integer
#' @param .sample A symbol
#'
#' @return A `tbl`
#'
#' @noRd
format_for_MPI = function(df, shards, .sample) {
	.sample = enquo(.sample)

	df %>%

		left_join((.) %>%
								distinct(G) %>%
								arrange(G) %>%
								mutate(idx_MPI = head(
									rep(seq_len(length.out=shards), (.) %>% nrow %>% `/` (shards) %>% ceiling), n = (.) %>% nrow
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
						mutate(`end` = cumsum(n)) %>%
						mutate(`start` = c(
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
				mutate(`symbol MPI row` = seq_len(length.out=n())) %>%
				ungroup,
			by = c("G", "idx_MPI")
		) %>%

		# Add counts MPI rows indexes
		group_by(idx_MPI) %>%
		arrange(G) %>%
		mutate(`read count MPI row` = seq_len(length.out=n())) %>%
		ungroup

}

#' add_partition
#'
#' @keywords internal
#'
#'
#' @description Add partition column dto data frame
#'
#' @param df.input A tibble
#' @param partition_by A symbol. Column we want to partition by
#' @param n_partitions An integer number of partition
#'
#' @return A `tbl`
#'
#' @noRd
add_partition = function(df.input, partition_by, n_partitions) {
	df.input %>%
		left_join(
			(.) %>%
				select(!!partition_by) %>%
				distinct %>%
				mutate(
					partition =
						seq_len(length.out=n() ) %>%
						divide_by(length((.))) %>%
						#	multiply_by(min(n_partitions, df.input %>% distinct(symbol) %>% nrow)) %>%
						multiply_by(n_partitions) %>%
						ceiling
				)
		)
}

#' Formula parser
#'
#' @importFrom stats terms
#'
#' @keywords internal
#'
#'
#' @param fm A formula
#'
#' @return A character vector
#'
#'
#' @noRd
parse_formula <- function(fm) {
	if (attr(terms(fm), "response") == 1)
		stop("The formula must be of the kind \"~ covariates\" ")
	else
		as.character(attr(terms(fm), "variables"))[-1]
}

#' vb_iterative
#'
#' @keywords internal
#'
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
#' @noRd
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
				pars=c("counts_rng", "exposure_rate", "alpha_sub_1", additional_parameters_to_save),
				...
			)
			boolFalse <- TRUE
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
#' @keywords internal
#'
#'
#' @return A Stan fit object
#' @noRd
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
#' @keywords internal
#'
#'
#' @importFrom tibble rowid_to_column
#' @importFrom purrr map
#'
#' @param counts_MPI A matrix of read count information
#' @param to_exclude A vector of oulier data points to exclude
#' @param shards An integer
#'
#' @return A matrix
#' @noRd
get_outlier_data_to_exlude = function(counts_MPI, to_exclude, shards) {
	# If there are genes to exclude
	switch(
		to_exclude %>% nrow %>% gt(0) %>% `!` %>% sum(1),
		foreach(s =  seq_len(length.out=shards), .combine = full_join) %do% {
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
#' @keywords internal
#' @importFrom stats rnorm
#'
#'
#' @return A list
#' @noRd
inits_fx =	function () {
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
#' @keywords internal
#'
#'
#' @importFrom purrr pmap
#' @importFrom purrr map_int
#' @importFrom purrr when
#' @import ggplot2
#'
#' @param .x A tibble
#' @param symbol A symbol object
#' @param .abundance A symbol object
#' @param .sample A symbol object
#' @param covariate A character string
#'
#' @return A ggplot
#' @noRd
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

	max_y  = .x %>% summarise(a = max(!!as.symbol(.abundance)), b = max(.upper_2)) %>% as.numeric %>% max

	{
		ggplot(data = .x, aes(
			y = !!as.symbol(.abundance),
			x = !!as.symbol(.sample)
		))
		# +
		# 	geom_errorbar(
		# 		aes(ymin = `.lower_1`,
		# 				ymax = `.upper_1`),
		# 		width = 0,
		# 		linetype = "dashed",
		# 		color = "#D3D3D3"
		# 	)
	} %>%
		when(
			".lower_2" %in% colnames(.x) ~ (.) +
				geom_errorbar(aes(
					ymin = `.lower_2`,
					ymax = `.upper_2`,
					color = `deleterious outliers`
				),
				width = 0),
			~ (.)

		) %>%
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
		scale_colour_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
		coord_cartesian(ylim = c(NA, max_y)) +
		my_theme +
		ggtitle(symbol)
}

# Add annotation if sample belongs to high or low group
add_deleterious_if_covariate_exists = function(.data, X){
	.data %>%
		when(
			X %>% ncol %>% gt(1) ~ left_join(.,
					X %>%
						as_tibble %>%
						select(2) %>%
						setNames("factor or interest") %>%
						mutate(S = seq_len(length.out=n())) %>%
						mutate(`is group right` = `factor or interest` > mean(`factor or interest`)) ,
					by = "S"
				) %>%

				mutate(`is group high` = (slope > 0 & `is group right`) |  (slope < 0 & !`is group right`)) %>%

				# Check if outlier might be deleterious for the statistics
				mutate(`deleterious outliers` = (!ppc) &
							 	(`is higher than mean` == `is group high`)),
			~ (.)
		)
}

#' merge_results
#'
#' @keywords internal
#'
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
#' @examples
#'
#' # Not needed because internal
#'
#' TRUE
#'
#' @return A `tbl`
#'
#' @noRd
merge_results = function(res_discovery, res_test, formula, .transcript, .abundance, .sample, do_check_only_on_detrimental){

	# Prepare column same enquo
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	res_discovery %>%
		filter(`.variable` == "counts_rng") %>%
		select(
			S,
			G,
			!!.transcript,
			!!.abundance,
			!!.sample,
			mean,
			# `.lower_1`,
			# `.upper_1`,
			`exposure rate`,
			slope_1 = slope,
			one_of(parse_formula(formula))
		) %>%

		# Attach results of tests
		left_join(
			res_test %>%
				filter(`.variable` == "counts_rng") %>%
					select(
						S,
						G,
						mean_2 = mean,
						.lower_2 = `.lower`,
						.upper_2 = `.upper`,
						slope_2 = slope,
						ppc,
						one_of(c("generated quantities", "deleterious outliers"))
					),
				by = c("S", "G")
		) %>%

		# format results
		format_results(formula, !!.transcript, !!.abundance, !!.sample, do_check_only_on_detrimental)


}

#' @importFrom utils packageVersion
format_results = function(.data, formula, .transcript, .abundance, .sample, do_check_only_on_detrimental){

	# Prepare column same enquo
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	.data %>%

		# Check if new package is installed with different sintax
		ifelse_pipe(
			packageVersion("tidyr") >= "0.8.3.9000",
			~ .x %>% nest(`sample wise data` = c(-!!.transcript)),
			~ .x %>%
				group_by(!!.transcript) %>%
				nest(`sample wise data` = -!!.transcript)
		) %>%

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
#' @keywords internal
#'
#'
#' @importFrom rstan sampling
#' @importFrom rstan vb
#' @importFrom utils tail
#'
#' @param .data A tibble
#' @param .do_check A boolean
#' @param .significance A symbol
#' @param .transcript A column name
#' @param how_many_negative_controls An integer
#'
#' @return A `tbl`
#'
#' @noRd
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
#' @keywords internal
#'
#'
#' @importFrom tidyr separate
#'
#' @param .data A data frame
#' @param fit A fit object
#'
#' @return A `tbl`
#'
#' @noRd
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
#' @keywords internal
#'
#'
#' @importFrom tidyr separate
#' @importFrom tidyr nest
#' @importFrom rstan summary
#'
#' @param fit A fit object
#' @param adj_prob_theshold fit real
#'
#' @examples
#'
#' # Not needed because internal
#'
#' TRUE
#'
#' @return A `tbl`
#' @noRd
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
#' @keywords internal
#'
#'
#' @importFrom tidyr separate
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @importFrom tidyr expand_grid
#' @importFrom rstan summary
#' @importFrom furrr future_map
#' @importFrom future multiprocess
#' @importFrom rstan extract
#' @importFrom stats sd
#'
#' @param fit A fit object
#' @param adj_prob_theshold fit real
#' @param how_many_posterior_draws An integer
#' @param truncation_compensation A real
#' @param cores An integer
#'
#' @examples
#'
#' # Not needed because internal
#'
#' TRUE
#'
#' @return A `tbl`
#'
#' @noRd
fit_to_counts_rng_approximated = function(fit, adj_prob_theshold, how_many_posterior_draws, truncation_compensation, cores, how_many_to_check){


	writeLines(sprintf("executing %s", "fit_to_counts_rng_approximated"))

	draws_mu = fit %>% rstan::extract("lambda_log_param") %>% .[[1]] %>% .[,,seq_len(length.out=how_many_to_check), drop=FALSE]

	# %>% as.data.frame() %>% setNames(sprintf("mu.%s", colnames(.))) %>%
	# 	as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, mu, -.draw) %>% separate(par, c("par", "S", "G"), sep="\\.") %>% select(-par)

	draws_sigma = fit %>% rstan::extract("sigma_raw") %>% .[[1]] %>% .[,seq_len(length.out=how_many_to_check), drop=FALSE]

	# %>% as.data.frame() %>% setNames(sprintf("sigma.%s", colnames(.) %>% gsub("V", "", .))) %>%
	# 	as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, sigma, -.draw) %>% separate(par, c("par", "G"), sep="\\.") %>% select(-par)

	draws_exposure = 	fit %>% rstan::extract("exposure_rate") %>% .[[1]]

	# %>% as.data.frame() %>% setNames(sprintf("exposure.%s", colnames(.) %>% gsub("V", "", .))) %>%
	# 	as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, exposure, -.draw) %>% separate(par, c("par", "S"), sep="\\.") %>% select(-par)

	expand_grid(S =seq_len(length.out=dim(draws_mu)[2]), G = seq_len(length.out=dim(draws_mu)[3])) %>%
		mutate(truncation_compensation = !!truncation_compensation) %>%
		mutate(
			CI = pmap(
				list(  S,G, truncation_compensation),
				~ {

					i_supersampled =	seq_len(length.out=sample(length(draws_mu[,..1, ..2]), how_many_posterior_draws, replace = TRUE ))
					draws = rnbinom(
						n = how_many_posterior_draws,
						mu = exp(draws_mu[,..1, ..2][i_supersampled] + draws_exposure[,..1][i_supersampled]),
						size = 1/exp(draws_sigma[,..2][i_supersampled]) * ..3
					)
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
		mutate(.variable = "counts_rng") %>%
		unnest(CI)

	# draws_mu =
	# 	fit %>% extract("lambda_log_param") %>% `[[` (1) %>% as.data.frame() %>% setNames(sprintf("mu.%s", colnames(.))) %>%
	# 	as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, mu, -.draw) %>% separate(par, c("par", "S", "G"), sep="\\.") %>% select(-par)
	# draws_sigma =
	# 	fit %>% extract("sigma_raw") %>% `[[` (1) %>% as.data.frame() %>% setNames(sprintf("sigma.%s", colnames(.) %>% gsub("V", "", .))) %>%
	# 	as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, sigma, -.draw) %>% separate(par, c("par", "G"), sep="\\.") %>% select(-par)
	# draws_exposure =
	# 	fit %>% extract("exposure_rate") %>% `[[` (1) %>% as.data.frame() %>% setNames(sprintf("exposure.%s", colnames(.) %>% gsub("V", "", .))) %>%
	# 	as_tibble() %>% mutate(.draw = 1:n()) %>% gather(par, exposure, -.draw) %>% separate(par, c("par", "S"), sep="\\.") %>% select(-par)
	#
	# draws_mu %>%
	# 	left_join(draws_sigma, by = c(".draw", "G")) %>%
	# 	left_join(draws_exposure, by = c(".draw", "S")) %>%
	# 	nest(data = -c(S, G)) %>%
	# 	mutate(
	# 		CI = map(
	# 			data,
	# 			~ {
	# 				.x_supersampled = .x %>%	sample_n(how_many_posterior_draws, replace = TRUE)
	# 				draws = rnbinom(n =how_many_posterior_draws,	mu = exp(.x_supersampled$mu + .x_supersampled$exposure),	size = 1/exp(.x_supersampled$sigma) * truncation_compensation	)
	# 				draws %>%
	# 					# Process quantile
	# 					quantile(c(adj_prob_theshold, 1 - adj_prob_theshold)) %>%
	# 					tibble::as_tibble(rownames="prop") %>%
	# 					tidyr::spread(prop, value) %>%
	# 					setNames(c(".lower", ".upper")) %>%
	# 					# Add mean and sd
	# 					dplyr::mutate(mean = mean(draws), sd = sd(draws))
	# 			}
	# 		)
	# 	) %>%
	# 	select(-data) %>%
	# 	unnest(CI) %>%
	#
	# 	# Adapt to old dataset
	# 	mutate(.variable = "counts_rng") %>%
	# 	mutate(S = as.integer(S), G = as.integer(G))


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
#' @keywords internal
#'
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
#' @return A `tbl`
#'
#' @noRd
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
		nrow %>% st(	.data %>% nrow	)
	)
		stop(sprintf("There are NA values in you tibble for any of the column %s", paste(columns, collapse=", ")))
}

#' @importFrom parallel detectCores
detect_cores = function(){

	detectCores()  %>% as.integer

	# if(.Platform$OS.type == "unix")
	# 	system("nproc", intern = TRUE) %>% as.integer %>% sum(-1)
	# else if(.Platform$OS.type == "windows")
	# 	parallel::detectCores()  %>% as.integer %>% sum(-1)
	# else stop("Your platform type is not recognised")

}

#' Create the design matrix
#'
#' @importFrom stats model.matrix
#'
#' @keywords internal
#'
#'
#' @param .data A tibble
#' @param formula A formula
#' @param .sample A symbol
#'
#' @examples
#'
#' # Not needed because internal
#'
#' TRUE
#'
#' @return A `matrix`
#'
#' @noRd
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
#' @keywords internal
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
#' @examples
#'
#' # Not needed because internal
#'
#' TRUE
#'
#' @return A `tbl`
#'
#' @noRd
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
								mutate(G = seq_len(length.out=n())),
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

#' draws_to_tibble_x_y
#'
#' @keywords internal
#'
#'
#' @return A `tbl`
#' @importFrom tidyr pivot_longer
#' @noRd
draws_to_tibble_x_y = function(fit, par, x, y) {

	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

	fit %>%
		rstan::extract(par_names, permuted=FALSE) %>%
		as.data.frame %>%
		as_tibble() %>%
		mutate(.iteration = seq_len(length.out=n())) %>%
		pivot_longer(names_to = c("dummy", ".chain", ".variable", x, y),  cols = contains(par), names_sep = "\\.|\\[|,|\\]|:", names_ptypes = list(".chain" = integer(), ".variable" = character(), "A" = integer(), "C" = integer()), values_to = ".value") %>%
		select(-dummy) %>%
		arrange(.variable, !!as.symbol(x), !!as.symbol(y), .chain) %>%
		group_by(.variable, !!as.symbol(x), !!as.symbol(y)) %>%
		mutate(.draw = seq_len(length.out=n())) %>%
		ungroup() %>%
		select(!!as.symbol(x), !!as.symbol(y), .chain, .iteration, .draw ,.variable ,     .value)

}

draws_to_tibble_x = function(fit, par, x) {

	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

	fit %>%
		rstan::extract(par_names, permuted=FALSE) %>%
		as.data.frame %>%
		as_tibble() %>%
		mutate(.iteration = seq_len(length.out=n())) %>%
		pivot_longer(names_to = c("dummy", ".chain", ".variable", x),  cols = contains(par), names_sep = "\\.|\\[|,|\\]|:", names_ptypes = list(".chain" = integer(), ".variable" = character(), "A" = integer(), "C" = integer()), values_to = ".value") %>%
		select(-dummy) %>%
		arrange(.variable, !!as.symbol(x), .chain) %>%
		group_by(.variable, !!as.symbol(x)) %>%
		mutate(.draw = seq_len(length.out=n())) %>%
		ungroup() %>%
		select(!!as.symbol(x), .chain, .iteration, .draw ,.variable ,     .value)

}

#' @importFrom stats sd
#' @importFrom purrr map_chr
identify_outliers_1_step = function(.data,
																		formula = ~ 1,
																		.sample,
																		.transcript,
																		.abundance,
																		.significance,
																		.do_check,
																		percent_false_positive_genes = 1,
																		how_many_negative_controls = 500,

																		approximate_posterior_inference = TRUE,
																		approximate_posterior_analysis = NULL,
																		draws_after_tail = 10,

																		save_generated_quantities = FALSE,
																		additional_parameters_to_save = c(),  # For development purpose
																		cores = detect_cores(), # For development purpose,
																		pass_fit = FALSE,
																		do_check_only_on_detrimental = length(parse_formula(formula)) > 0,
																		tol_rel_obj = 0.01,
																		just_discovery = FALSE,
																		seed = sample(seq_len(length.out=999999), size = 1),
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
	if (percent_false_positive_genes %>% is.na |
			!(percent_false_positive_genes %>% between(0, 100)))
		stop("percent_false_positive_genes must be a string from > 0% to < 100%")

	# For reference MPI inference
	# Check if all trannscripts are non NA
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
	adj_prob_theshold_1  = 0.05
	if(adj_prob_theshold_2 %>% is.null)
		adj_prob_theshold_2 =
		percent_false_positive_genes / 100 /
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

	# Scale dataset
	my_df_scaled = scale_abundance(my_df, !!.sample,!!.transcript,!!.abundance)

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
			adj_prob_theshold  = adj_prob_theshold_2,
			how_many_posterior_draws = how_many_posterior_draws_2,
			pass_fit = pass_fit,
			tol_rel_obj = tol_rel_obj,
			write_on_disk = write_on_disk,
			seed = seed
		) %>%

		# For building some figure I just need the discovery run, return prematurely
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
			ppc,
			one_of(c("generated quantities", "deleterious outliers")),
			`exposure rate`,
			one_of(parse_formula(formula))
		) %>%

		# format results
		format_results(formula, !!.transcript, !!.abundance, !!.sample, do_check_only_on_detrimental)

}

summary_to_tibble = function(fit, par, x, y = NULL) {

	par_names = names(fit) %>% grep(sprintf("%s", par), ., value = TRUE)

	fit %>%
		rstan::summary(par_names) %$%
		summary %>%
		as_tibble(rownames = ".variable") %>%
		when(
			is.null(y) ~ (.) %>% tidyr::separate(.variable, c(".variable", x), sep="\\[|\\]", extra = "drop", convert=TRUE),
			~ (.) %>% tidyr::separate(.variable, c(".variable", x, y), sep="\\[|\\]|,", extra = "drop", convert=TRUE)
		)

}

#' do_inference
#'
#' @keywords internal
#'
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
#' @importFrom tidybayes gather_draws
#' @importFrom stats sd
#' @importFrom stats start
#' @importFrom stats end
#' @importFrom utils head
#'
#' @param my_df A tibble including a transcript name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param .sample A column name as symbol
#' @param .transcript A column name as symbol
#' @param .abundance A column name as symbol
#' @param .significance A column name as symbol
#' @param .do_check A column name as symbol
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
#' @param seed an integer
#'
#' @return A tibble with additional columns
#'
#' @noRd
do_inference = function(my_df,
												formula,
												.sample ,
												.transcript ,
												.abundance,
												.significance ,
												.do_check,
												approximate_posterior_inference = FALSE,
												approximate_posterior_analysis = FALSE,
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
												save_generated_quantities = FALSE,
												inits_fx = "random",
												prior_from_discovery = tibble(`.variable` = character(),
																											mean = numeric(),
																											sd = numeric()),
												pass_fit = FALSE,
												tol_rel_obj = 0.01,
												write_on_disk = FALSE,
												seed) {

	writeLines(sprintf("executing %s", "do_inference"))

	# Prepare column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	.significance = enquo(.significance)
	.do_check = enquo(.do_check)

	# Check that the dataset is squared
	if(my_df %>% distinct(!!.sample, !!.transcript) %>% count(!!.transcript) %>% count(n) %>% nrow %>% gt(1))
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
	additional_parameters_to_save = additional_parameters_to_save %>% c("lambda_log_param", "sigma_raw") %>% unique

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
	M = counts_MPI %>% distinct(`start`, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
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
		distinct(idx_MPI, `end`, `symbol MPI row`)  %>%
		spread(idx_MPI, `end`) %>%
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
		matrix(nrow = shards, byrow = TRUE) %>%
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
					"alpha_sub_1",
					additional_parameters_to_save
				)
			)
		)

	# Parse and return
	fit %>%

		ifelse_pipe(
			approximate_posterior_analysis,
			~ .x %>% fit_to_counts_rng_approximated(adj_prob_theshold, how_many_posterior_draws, truncation_compensation, cores, how_many_to_check),
			~ .x %>% fit_to_counts_rng(adj_prob_theshold)
		) %>%

		# If generated quantities are saved
		save_generated_quantities_in_case(fit, save_generated_quantities) %>%

		# Add exposure rate
		#add_exposure_rate(fit) %>%

		# Check if data is within posterior
		check_if_within_posterior(my_df, .do_check, .abundance) %>%

		# Attach slope
		left_join(summary_to_tibble(fit, "alpha_sub_1", "G") %>% select(G, slope = mean), by = "G") %>%

		# Add annotation if sample belongs to high or low group
		add_deleterious_if_covariate_exists(X) %>%

		# Add position in MPI package for next inference
		left_join(counts_MPI %>% distinct(S, G, idx_MPI, `read count MPI row`),
							by = c("S", "G")) %>%

		# Add exposure rate
		add_exposure_rate(fit) %>%

		# needed for the figure article
		ifelse_pipe(pass_fit,	~ .x %>% add_attr(fit, "fit")	) %>%

		# Passing the amount of sampled data
		add_attr(S * how_many_to_check * how_many_posterior_draws, "total_draws")


}


