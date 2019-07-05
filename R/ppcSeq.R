#' format_for_MPI
#'
#' @description Format reference data frame for MPI
format_for_MPI = function(df, shards){
	df %>%

		left_join(
			(.) %>%
				distinct(G) %>%
				arrange(G) %>%
				mutate( idx_MPI = head( rep(1:shards, (.) %>% nrow %>% `/` (shards) %>% ceiling ), n=(.) %>% nrow) )
		) %>%
		arrange(idx_MPI, G) %>%

		# Decide start - end location
		group_by(idx_MPI) %>%
		do(
			(.) %>%
				left_join(
					(.) %>%
						distinct(sample, G) %>%
						arrange(G) %>%
						count(G) %>%
						mutate(end = cumsum(n)) %>%
						mutate(start = c(1, .$end %>% rev() %>% `[` (-1) %>% rev %>% `+` (1)))
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
				ungroup
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
#' @import tidyr
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
				dplyr::summarise_all(class) %>%
				tidyr::gather(variable, class) %>%
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
					magrittr::set_rownames(tbl %>% pull(!!rownames)) %>%
					dplyr::select(-!!rownames)
			} else {
				(.)
			}
		} %>%

		# Convert to matrix
		as.matrix()
}

#' pcc_seq main
#'
#' @description This function calls the stan model.
#'
#'
#' @importFrom tibble tibble
#'
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr mutate
#' @importFrom dplyr filter
#' @importFrom dplyr mutate_if
#'
#' @importFrom tidyr spread
#' @importFrom tidyr gather
#' @importFrom tidyr drop_na
#'
#' @importFrom tidybayes gather_samples
#' @importFrom tidybayes median_qi
#'
#' @importFrom foreach foreach
#' @importFrom foreach %do%
#'
#' @param mix A matrix
#' @param my_design A matrix
#' @param cov_to_test A character string
#' @param fully_bayesian A boolean
#' @param is_mix_microarray A boolean
#' @param ct_to_omit A character string
#' @param verbose A boolean
#' @param save_report A boolean
#' @param custom_ref A matrix
#' @param multithread A boolean
#' @param do_debug A boolean
#' @param cell_type_root A character string
#' @param choose_internal_ref A design matrix
#' @param omit_regression A boolean
#' @param save_fit A boolean
#' @param seed An integer
#'
#' @return An ARMET object
#'
#' @export
#'
pcc_seq = function(
	input.df,
	formula,
	sample_column = "sample",
	gene_column = "symbol",
	value_column = "read count"
){

	library(rstan)

	input = c(as.list(environment()))
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
			aspect.ratio=1,
			axis.text.x = element_text(angle = 90, hjust = 1),
			strip.background = element_blank(),
			axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
			axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
		)

	#########################################
	# For  reference MPI inference

	if(input.df %>% filter(!!as.symbol(gene_column) %>% is.na) %>% nrow > 0) stop("There are NAs in the gene_column. Please filter those records")

	# distinct_at is not released yet for dplyr, thus we have to use this trick
	my_df <- input.df %>%
		select(!!gene_column, !!sample_column, !!value_column, one_of(parse_formula(formula))) %>%
		setNames(c("symbol", "sample", "read count", parse_formula(formula))) %>%
		distinct() %>%

		# Add symbol idx
		left_join(
			(.) %>%
				distinct(symbol) %>%
				mutate(G = 1:n())
		) %>%

		# Add sample indeces
		mutate(S = factor(sample, levels = .$sample %>% unique) %>% as.integer)


	# Create design matrix
	X =
		model.matrix(
			object = formula,
			data = my_df %>% select(sample, one_of(parse_formula(formula))) %>% distinct %>% arrange(sample)
		)
	C = X %>% ncol
	#%>%
	#	magrittr::set_colnames(c("(Intercept)", . %>% gsub parse_formula(formula)))

	counts_MPI =
		my_df %>%
		select(symbol, sample, `read count`, S, G) %>%
		format_for_MPI(shards)

	G = counts_MPI %>% distinct(G) %>% nrow()
	S = counts_MPI %>% distinct(sample) %>% nrow()
	N = counts_MPI %>% distinct(idx_MPI, `read count`, `read count MPI row`) %>%  count(idx_MPI) %>% summarise(max(n)) %>% pull(1)
	M = counts_MPI %>% distinct(start, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% max
	G_per_shard = counts_MPI %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% as.array
	n_shards = min(shards, counts_MPI %>% distinct(idx_MPI) %>% nrow)
	G_per_shard_idx = c(0, counts_MPI %>% distinct(symbol, idx_MPI) %>% count(idx_MPI) %>% pull(n) %>% cumsum)

	counts =
		counts_MPI %>%
		distinct(idx_MPI, `read count`, `read count MPI row`)  %>%
		spread(idx_MPI,  `read count`) %>%
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
browser()
	########################################
	# MODEL

	# fileConn<-file("~/.R/Makevars")
	# writeLines(c( "CXX14FLAGS += -O3","CXX14FLAGS += -DSTAN_THREADS", "CXX14FLAGS += -pthread"), fileConn)
	# close(fileConn)
	# Sys.setenv("STAN_NUM_THREADS" = cores)
	# pcc_seq_model = stan_model("inst/stan/negBinomial_MPI.stan")

	Sys.time() %>% print
	fit =
		sampling(
			stanmodels$negBinomial_MPI, #pcc_seq_model, #
			chains=3, cores=3,
			iter=600, warmup=500,   save_warmup = FALSE
		)
	Sys.time() %>% print

	Sys.time() %>% print
	fit_vb =
		vb(
			stanmodels$negBinomial_MPI, #pcc_seq_model, #
			output_samples=2000,
			iter = 50000,
			tol_rel_obj=0.001

		)
	Sys.time() %>% print



	########################################
	# Parse results

	# Relation expected value, variance
	fit %>%
		tidybayes::spread_draws(lambda_log_param[G], sigma_raw_param[G]) %>%
		ggplot(aes(x=lambda_log_param, y=sigma_raw_param, group=G, alpha=)) +
		stat_ellipse( alpha=0.2) +
		my_theme

	fit %>%
		tidybayes::spread_draws(counts_rng[S,G]) %>%
		filter(G==1) %>%
		ggplot(aes(x=counts_rng+1, group=S)) +
		geom_density(fill="grey") +
		geom_vline(data = my_df %>% filter(G==1), aes(xintercept = `read count`, color=ct),
			linetype="dotted",
			size=1.5
		) +
		facet_wrap(~ S) +
		my_theme

	fit %>%
		tidybayes::spread_draws(counts_rng[S,G]) %>%
		tidybayes::median_qi() %>%
		left_join(my_df) %>%
		rowwise() %>%
		mutate(`ppc` = `read count` %>% between(`.lower`, `.upper`)) %>%
		filter(G==30) %>%
		ggplot(aes(y=`read count`, x=sample)) +
		geom_errorbar(aes(ymin=`.lower`, ymax=`.upper`, color=ppc)) +
		geom_point(aes(color=ct)) +
		my_theme

	# Return
	list(

		# Return the input itself
		input = input,

		# Return the fitted object
		fit = fit
	)

}
