#' format_for_MPI
#'
#' @description Format reference data frame for MPI
format_for_MPI = function(df, shards, sample_column){

	sample_column = enquo(sample_column)

	df %>%

		left_join(
			(.) %>%
				distinct(G) %>%
				arrange(G) %>%
				mutate( idx_MPI = head( rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling ), n=(.) %>% nrow) ),
			by = "G"
		) %>%
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
						mutate(start = c(1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1))),
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
add_partition = function(df.input, partition_by, n_partitions){
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
#' @return A character vector
#'
#'
parse_formula <- function(fm) {
	if(attr(terms(fm), "response") == 1) stop("The formula must be of the kind \"~ covariates\" ")
	else as.character(attr(terms(fm), "variables"))[-1]
}

#' Get matrix from tibble
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom magrittr set_rownames
#'
#' @param tbl A tibble
#' @param rownames A character string of the rownames
#' @return A matrix
as_matrix <- function(tbl, rownames = NULL) {
	tbl %>%

		# Check if data frame is not numerical beside the rownames column (if present)
		{
			if (
				tbl %>% {
					if (!is.null(rownames)) (.) %>% dplyr::select(-contains(rownames)) else (.)
				} %>%
				summarise_all(class) %>%
				gather(variable, class) %>%
				pull(class) %>%
				unique() %>%
				`%in%`(c("numeric", "integer")) %>% `!`() %>% any()
				# identical("numeric")
			) {
				warning("to_matrix says: there are NON-numerical columns, the matrix will NOT be numerical")
			}

			(.)
		} %>%
		as.data.frame() %>%

		# Deal with rownames column if present
		{
			if (!is.null(rownames)) {
				(.) %>%
					set_rownames(tbl %>% pull(!!rownames)) %>%
					select(-!!rownames)
			} else {
				(.)
			}
		} %>%

		# Convert to matrix
		as.matrix()
}

other_code = function(){

	# Relation expected value, variance
	# fit_draws = fit@sim$samples[[1]] %>% bind_rows() %>% mutate(.draw = 1:n()) %>% gather(.variable, .value, -.draw)
	#
	# # Plot relation lambda sigma
	# fit_draws %>% filter(grepl("^alpha", .variable)) %>% separate(.variable, c(".variable", "C", "G"), sep="\\.") %>%
	# 	mutate(C = C %>% as.integer, G = G %>% as.integer) %>%
	# 	filter(C == 1) %>%
	# 	select(-C) %>%
	# 	bind_rows(
	# 		fit_draws %>% filter(grepl("^sigma_raw_param", .variable)) %>% separate(.variable, c(".variable", "G"), sep="\\.") %>%
	# 			mutate( G = G %>% as.integer)
	# 	) %>%
	# 	spread(.variable, .value) %>%
	# 	ggplot(aes(x=alpha, y=sigma_raw_param, group=G)) +
	# 	stat_ellipse( alpha=0.2) +
	# 	my_theme

	# fit %>%
	# 	tidybayes::spread_draws(counts_rng[S,G]) %>%
	# 	filter(G==1) %>%
	# 	ggplot(aes(x=counts_rng+1, group=S)) +
	# 	geom_density(fill="grey") +
	# 	geom_vline(data = my_df %>% filter(G==1), aes(xintercept = !!value_column, color=ct),
	# 		linetype="dotted",
	# 		size=1.5
	# 	) +
	# 	facet_wrap(~ S) +
	# 	my_theme
}

#' pcc_seq main
#'
#' @description This function calls the stan model.
#'
#' @importFrom tibble as_tibble
#' @importFrom rstan sampling
#' @importFrom rstan vb
#' @importFrom rstan summary
#' @import dplyr
#' @importFrom tidyr spread
#' @import tidybayes
#' @importFrom magrittr %$%
#' @importFrom purrr map2
#' @importFrom purrr map_int
#' @importFrom tidyTranscriptomics add_normalised_counts
#'
#' @param input.df A tibble including a gene name column | sample name column | read counts column | covariates column
#' @param formula A formula
#' @param sample_column A column name
#' @param gene_column A column name
#' @param value_column A column name
#' @param significance_column A column name
#' @param full_bayes A boolean
#' @param how_many_negative_controls An integer
#' @param how_many_posterior_draws An integer

#'
#' @return A tibble with additional columns
#'
#' @export
#'
ppc_seq = function(
	input.df,
	formula = ~ 1,
	sample_column = `sample`,
	gene_column = `symbol`,
	value_column = `read count`,
	significance_column = `p-value`,
	do_check_column,
	full_bayes = T,
	how_many_negative_controls = 500,
	how_many_posterior_draws = 500
){

	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)

	#input = c(as.list(environment()))
	cores = 30/3 %>% floor
	shards = cores * 2

	my_theme =
		theme_bw() +
		theme(
			panel.border = element_blank(),
			axis.line = element_line(),
			panel.grid.major = element_line(size = 0.2),
			panel.grid.minor = element_line(size = 0.1),
			text = element_text(size=12),
			legend.position="bottom",
			#aspect.ratio=1,
			axis.text.x = element_text(angle = 90, hjust = 0.5, vjust=1),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		)

	#########################################
	# For  reference MPI inference

	if(input.df %>% filter(!!gene_column %>% is.na) %>% nrow > 0) stop("There are NAs in the gene_column. Please filter those records")

	if(input.df %>% select(!!value_column) %>% sapply(class) != "integer") stop("The algorithm takes raw (un-normalised) integer read counts only")

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	my_df <- input.df %>%

		# Select only significant genes plus background
		{
			bind_rows(

				# Genes to check
				(.) %>%
					filter((!!do_check_column)),

				# Least changing genes
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
		} %>%

		select(!!gene_column, !!sample_column, !!value_column, one_of(parse_formula(formula)), !!do_check_column) %>%
		#setNames(c("symbol", "sample", "read count", parse_formula(formula), do_check_column)) %>%
		distinct() %>%

		# Add symbol idx
		left_join(
			(.) %>%
				distinct(!!gene_column) %>%
				mutate(G = 1:n()),
			by = quo_name(gene_column)
		) %>%

		# Add sample indeces
		mutate(S = factor(!!sample_column, levels = (.) %>% pull(!!sample_column) %>% unique) %>% as.integer)

	how_many_to_check =
		ifelse(
			parse_formula(formula) %>% length > 0,
			input.df %>% filter(!!do_check_column) %>% select(!!gene_column) %>% distinct() %>% nrow,
			0
		)

	# Create design matrix
	X =
		model.matrix(
			object = formula,
			data = my_df %>% select(!!sample_column, one_of(parse_formula(formula))) %>% distinct %>% arrange(!!sample_column)
		)
	C = X %>% ncol
	#%>%
	#	magrittr::set_colnames(c("(Intercept)", . %>% gsub parse_formula(formula)))

	counts_MPI =
		my_df %>%
		select(!!gene_column, !!sample_column, !!value_column, S, G) %>%
		format_for_MPI(shards, !!sample_column)

	G = counts_MPI %>% distinct(G) %>% nrow()
	S = counts_MPI %>% distinct(!!sample_column) %>% nrow()
	N = counts_MPI %>% distinct(idx_MPI, !!value_column, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_MPI %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G_per_shard = counts_MPI %>% distinct(!!gene_column, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array
	n_shards = min(shards, counts_MPI %>% distinct(idx_MPI) %>% nrow)
	G_per_shard_idx = c(0, counts_MPI %>% distinct(!!gene_column, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum)

	counts =
		counts_MPI %>%
		distinct(idx_MPI, !!value_column, `read count MPI row`)  %>%
		spread(idx_MPI,  !!value_column) %>%
		select(-`read count MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	sample_idx =
		counts_MPI %>%
		distinct(idx_MPI, S, `read count MPI row`)  %>%
		spread(idx_MPI, S) %>%
		select(-`read count MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	symbol_end =
		counts_MPI %>%
		distinct(idx_MPI, end, `symbol MPI row`)  %>%
		spread(idx_MPI, end) %>%
		bind_rows( (.) %>% head(n=1) %>%  mutate_all(function(x) {0}) ) %>%
		arrange(`symbol MPI row`) %>%
		select(-`symbol MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	G_ind =
		counts_MPI %>%
		distinct(idx_MPI, G, `symbol MPI row`)  %>%
		spread(idx_MPI, G) %>%
		arrange(`symbol MPI row`) %>%
		select(-`symbol MPI row`) %>%
		replace(is.na(.), 0 %>% as.integer) %>%
		as_matrix() %>% t

	counts_package =
		# Dimensions data sets
		rep(c(M, N, S), shards) %>%
		matrix(nrow = shards, byrow = T) %>%
		cbind(G_per_shard) %>%
		cbind(symbol_end) %>%
		cbind(sample_idx) %>%
		cbind(counts)
	CP = ncol(counts_package)

	########################################
	# Prior info

	lambda_mu_mu = 5.612671

	########################################
	# Build better scales

	exposure_rate_multiplier =
		my_df %>%
		add_normalised_counts(!!sample_column, !!gene_column, !!value_column) %>%
		distinct(!!sample_column, TMM, multiplier) %>%
		mutate(l = multiplier %>% log) %>%
		summarise(l %>% sd) %>%
		pull(`l %>% sd`)

	intercept_shift_scale =
		my_df %>%
		add_normalised_counts(!!sample_column, !!gene_column, !!value_column) %>%
		mutate(
			cc =
				!!as.symbol(sprintf("%s normalised",  quo_name(value_column))) %>%
				`+` (1) %>% log
					 ) %>%
		summarise(shift = cc %>% mean, scale = cc %>% sd) %>%
		as.numeric

	# sigma_shift_scale =
	# 	my_df %>%
	# 	add_normalised_counts() %>%
	# 	mutate(cc = `read count normalised` %>% `+` (1) %>% log) %>%
	# 	group_by(symbol) %>%
	# 	summarise(sigma = cc %>% sd %>% log) %>%
	# 	ungroup() %>%
	# 	summarise(shift = sigma %>% mean, scale = sigma %>% sd)

	########################################
	# MODEL
	Sys.setenv("STAN_NUM_THREADS" = cores)

	# fileConn<-file("~/.R/Makevars")
	# writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	# close(fileConn)
	# pcc_seq_model = stan_model("inst/stan/negBinomial_MPI.stan")
	# fit = vb(
	# 	pcc_seq_model, #
	# 	output_samples=1000,
	# 	iter = 50000,
	# 	tol_rel_obj=0.001
	# )

	Sys.time() %>% print
	fit =
		switch(
			full_bayes %>% `!` %>% as.integer %>% sum(1),
			sampling(
				stanmodels$negBinomial_MPI, #pcc_seq_model, #
				chains=3, cores=3,
				iter=how_many_posterior_draws + 150, warmup=150,   save_warmup = FALSE, pars=c("counts_rng", "exposure_rate")
			),
			vb(
				stanmodels$negBinomial_MPI, #pcc_seq_model, #
				output_samples=1000,
				iter = 50000,
				tol_rel_obj=0.005, pars=c("counts_rng", "exposure_rate")
			)
		)
	Sys.time() %>% print

	########################################
	# Parse results
	# Return
	input.df %>%
		left_join(

			# Parse fit
			fit %>%
				summary("counts_rng") %$%
				summary %>%
				as_tibble(rownames = ".variable") %>%
				separate(.variable, c(".variable", "S", "G"), sep="[\\[,\\]]") %>%
				mutate(S = S %>% as.integer, G = G %>% as.integer) %>%

				# Add exposure rate
				left_join(
					fit %>%
						summary("exposure_rate") %$%
						summary %>%
						as_tibble(rownames = ".variable") %>%
						separate(.variable, c(".variable", "S"), sep="[\\[,\\]]") %>%
						mutate(S = S %>% as.integer) %>%
						rename(`exposure rate` = mean) %>%
						select(S, `exposure rate`)
				) %>%

				# Check if data is within posterior
				left_join(my_df) %>%
				filter((!!do_check_column)) %>% # Filter only DE genes
				rowwise() %>%
				mutate(`ppc` = !!value_column %>% between(`2.5%`, `97.5%`)) %>%
				ungroup %>%

				# Add plots
				group_by(!!gene_column) %>%
				nest %>%
				mutate(plot = map2(data, !!gene_column, ~
													 	{
													 		ggplot(data = .x, aes(y=!!value_column, x=!!sample_column)) +
													 			geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`, color =ppc), width=0) +
													 			scale_colour_manual(values = c( "FALSE"= "red", "TRUE"= "black")) +
													 			my_theme
													 	} %>%
													 	{
													 		if(parse_formula(formula)[1] %>% is.null %>% `!`)
													 			(.) + geom_point(aes(size=`exposure rate`, fill = !!as.symbol(parse_formula(formula)[1])), shape = 21)
													 		else
													 			(.) + geom_point(aes(size=`exposure rate`), shape=21, fill="black")
													 	}
				)) %>%

				# Add summary statistics
				mutate(
					`ppc samples failed` = map_int(data, ~ .x %>% pull(ppc) %>% `!` %>% sum)
				) %>%
				#rename(!!gene_column := !!gene_column) %>%
				select(-data)
		)
}

