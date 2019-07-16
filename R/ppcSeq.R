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


do_inference = function(
	my_df,
	formula,
	sample_column ,
	gene_column ,
	value_column,
	significance_column ,
	do_check_column,
	full_bayes = T,
	C,
	X,
	lambda_mu_mu,
	cores,
	exposure_rate_multiplier, intercept_shift_scale,
	additional_parameters_to_save,
	adj_prob_theshold,
	to_exclude = tibble(S = integer(), G = integer()),
	save_generated_quantities = F
){

	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)


	how_many_to_check =
		my_df %>%
		filter(!!do_check_column) %>%
		distinct(!!gene_column) %>%
		nrow

	# Calculate the needed posterior draws
	how_many_posterior_draws =  5 %>% divide_by(adj_prob_theshold) %>% max(500)

	chains =
		foreach(cc = 2:min(cores, 6), .combine = bind_rows) %do%
		{ tibble(chains = cc, tot = how_many_posterior_draws / cc + 150 * cc )  } %>%
		filter(tot == tot %>% min) %>%
		pull(chains)

	my_cores = cores %>% divide_by(chains) %>% floor
	shards = my_cores

	counts_MPI =
		my_df %>%
		select(!!gene_column, !!sample_column, !!value_column, S, G) %>%
		format_for_MPI(shards, !!sample_column)

	to_exclude_MPI =
		switch(
			to_exclude %>% nrow %>% `>` (0) %>% `!` %>% sum(1), # If there are genes to exclude
			foreach(s = 1:shards, .combine=full_join) %do% {
				counts_MPI %>%
					inner_join(to_exclude, by=c("S", "G")) %>%
					filter(idx_MPI == s) %>%
					distinct(idx_MPI, `read count MPI row`) %>%
					rowid_to_column %>%
					spread(idx_MPI, `read count MPI row`)

			} %>%

				# Add length array to the first row for indexing in MPI
				{
					bind_rows(
						(.) %>% map(function(x) x %>% is.na %>% `!` %>% as.numeric %>% sum) %>% unlist,
						(.)
					)
				} %>%
				select(-rowid) %>%
				replace(is.na(.), 0 %>% as.integer) %>%
				as_matrix() %>% t,
			matrix(rep(0,shards))
		)


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
		cbind(counts) %>%
		cbind(to_exclude_MPI)

	CP = ncol(counts_package)

	######################################
	# Run model

	Sys.setenv("STAN_NUM_THREADS" = my_cores)


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
				chains=chains, cores=chains,
				iter=(how_many_posterior_draws/chains) %>% ceiling %>% sum(150), warmup=150,   save_warmup = FALSE, pars=c("counts_rng", "exposure_rate", additional_parameters_to_save)
			),
			vb(
				stanmodels$negBinomial_MPI, #pcc_seq_model, #
				output_samples=how_many_posterior_draws,
				iter = 50000,
				tol_rel_obj=0.005, pars=c("counts_rng", "exposure_rate", additional_parameters_to_save)
			)
		)
	Sys.time() %>% print

	########################################
	# Parse results
	# Return

	# Parse fit
	fit %>%
		summary("counts_rng", prob=c(adj_prob_theshold, 1-adj_prob_theshold)) %$%
		summary %>%
		as_tibble(rownames = ".variable") %>%
		separate(.variable, c(".variable", "S", "G"), sep="[\\[,\\]]", extra="drop") %>%
		mutate(S = S %>% as.integer, G = G %>% as.integer) %>%
		rename(`.lower` = 7, `.upper` = 8) %>%

		{
			if(save_generated_quantities)
				(.) %>%
				left_join(
					fit %>% tidybayes::gather_draws(counts_rng[S,G])
				) %>%
				group_by(`.variable`, S ,G, mean,se_mean , sd, `.lower`, `.upper`, n_eff,  Rhat) %>%
				nest(.key = "generated quantities")
			else
				(.)
		} %>%

		# Add exposure rate
		left_join(
			fit %>%
				summary("exposure_rate") %$%
				summary %>%
				as_tibble(rownames = ".variable") %>%
				separate(.variable, c(".variable", "S"), sep="[\\[,\\]]", extra="drop") %>%
				mutate(S = S %>% as.integer) %>%
				rename(`exposure rate` = mean) %>%
				select(S, `exposure rate`),
			by="S"
		) %>%

		# Check if data is within posterior
		left_join(my_df, by = c("S", "G")) %>%
		filter((!!do_check_column)) %>% # Filter only DE genes
		rowwise() %>%
		mutate(`ppc` = !!value_column %>% between(`.lower`, `.upper`)) %>%
		mutate(`is higher than mean` = (!`ppc`) & (!!value_column > mean)) %>%
		ungroup %>%

		# Add annotation if sample belongs to high or low group
		left_join(
			X %>%
				as_tibble %>%
				select(2) %>%
				setNames("factor or interest") %>%
				mutate(S=1:n()) %>%
				mutate(`is group high` = `factor or interest` > mean(`factor or interest`)),
			by = "S"
		) %>%

		# Check if outlier might be deleterious for the statistics
		mutate(`deleterious outliers` = (!ppc) & (`is higher than mean` == `is group high`)) %>%

		# Add position in MPI package for next inference
		left_join(	counts_MPI %>% distinct(S, G, idx_MPI, `read count MPI row`)	)
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
#' @importFrom magrittr divide_by
#' @importFrom magrittr multiply_by
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
	save_generated_quantities = F,       # For development purpose
	additional_parameters_to_save = c(), # For development purpose,
	cores = system("nproc", intern = TRUE) %>% as.integer %>% sum(-1),
	error_rate = 0.05
){

	sample_column = enquo(sample_column)
	gene_column = enquo(gene_column)
	value_column = enquo(value_column)
	significance_column = enquo(significance_column)
	do_check_column = enquo(do_check_column)

	#input = c(as.list(environment()))


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

	if(input.df %>% select(!!value_column) %>% sapply(class) != "integer")
		stop(sprintf("The column %s must be of class integer. You can do as mutate(`%s` = `%s` %>% as.integer)", quo_name(value_column), quo_name(value_column), quo_name(value_column)))

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

	# Create design matrix
	X =
		model.matrix(
			object = formula,
			data = my_df %>% select(!!sample_column, one_of(parse_formula(formula))) %>% distinct %>% arrange(!!sample_column)
		)
	C = X %>% ncol
	#%>%
	#	magrittr::set_colnames(c("(Intercept)", . %>% gsub parse_formula(formula)))

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

	########################################
	# MODEL

	res_discovery =
		my_df %>%
		do_inference(
			formula,
			!!sample_column ,
			!!gene_column ,
			!!value_column ,
			!!significance_column ,
			!!do_check_column,
			full_bayes,
			C,
			X,
			lambda_mu_mu,
			cores,
			exposure_rate_multiplier,
			intercept_shift_scale,
			additional_parameters_to_save,
			adj_prob_theshold  = error_rate
		)

	# Columns of counts to be ignored from the inference
	to_exclude =
		res_discovery %>%
		filter(`deleterious outliers`) %>%
		distinct(S, G)

	how_namy_to_exclude = to_exclude %>% nrow

	res_test =
		my_df %>%
		do_inference(
			formula,
			!!sample_column ,
			!!gene_column ,
			!!value_column ,
			!!significance_column ,
			!!do_check_column,
			full_bayes,
			C,
			X,
			lambda_mu_mu,
			cores,
			exposure_rate_multiplier,
			intercept_shift_scale,
			additional_parameters_to_save,
			adj_prob_theshold = error_rate / how_namy_to_exclude * 2, # * 2 because we just test one side of the distribution
			to_exclude = to_exclude,
			save_generated_quantities = save_generated_quantities
		)

		# Merge results
		res_discovery %>%
			select(
				S, G, !!gene_column, !!value_column, !!sample_column,
				`.lower`, `.upper`, `exposure rate`, !!as.symbol(parse_formula(formula)[1])
			) %>%
			left_join(
				res_test %>%
					select(S, G, `.lower`, `.upper`, `deleterious outliers`) %>%
					rename(`.lower_2` = `.lower`, `.upper_2` = `.upper`),
				by = c("S", "G")
			) %>%

		# Rpoduce summary results
			# Add plots
			group_by(!!gene_column) %>%
			nest(.key = "sample wise data") %>%
			mutate(plot = map2(`sample wise data`, !!gene_column, ~
												 	{
												 		ggplot(data = .x, aes(y=!!value_column, x=!!sample_column)) +
												 			geom_errorbar(aes(
												 				ymin=`.lower`,
												 				ymax=`.upper`
												 			), width=0, linetype="dashed", color="#D3D3D3") +
												 			geom_errorbar(aes(
												 				ymin=`.lower_2`,
												 				ymax=`.upper_2`,
												 				color =`deleterious outliers`
												 			), width=0) +
												 			scale_colour_manual(values = c( "TRUE"= "red", "FALSE"= "black")) +
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
				# `ppc samples failed` = map_int(data, ~ .x %>% pull(ppc) %>% `!` %>% sum),
				`tot deleterious outliers` = map_int(`sample wise data`, ~ .x %>% pull(`deleterious outliers`) %>% sum)
			)
			#%>%
			#rename(!!gene_column := !!gene_column) %>%
			#select(-data)


}

