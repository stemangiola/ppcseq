library(tidyverse)
library(magrittr)
library(ttBulk)
library(ppcSeq)

#TCGA_tbl = readRDS("/wehisan/bioinf/bioinf-data/Papenfuss_lab/projects/mangiola.s/PostDoc/temp_TCGA_tbl.RData")
# load("dev/TCGA_tbl.RData")

# Interctept only

ifelse_pipe = function(.x, .p, .f1, .f2 = NULL) {
	switch(.p %>% `!` %>% sum(1),
				 as_mapper(.f1)(.x),
				 if (.f2 %>% is.null %>% `!`)
				 	as_mapper(.f2)(.x)
				 else
				 	.x)

}

create_input_df = function(df, my_run){
	df %>%
		# Filter
		filter(`CAPRA-S` %>% is.na %>% `!`) %>%
		separate(sample, c("data base", "laboratory", "patient"), sep="-", remove = F) %>%
		inner_join(
			(.) %>% distinct(sample, laboratory) %>% count(laboratory) %>% filter(n >= 8) %>% select(-n)
		) %>%
		mutate(risk = `CAPRA-S` <= 3) %>%

		# Do check
		mutate(do_check = (!`house keeping`) & run==my_run)
}

get_NB_qq_values = function(input.df) {
	input.df = input.df %>% arrange(`read count normalised adjusted`)

	predicted_NB =
		qnbinom(
			# If 1 sample, just use median
			switch(
				input.df %>% nrow %>% `>` (1) %>% `!` %>% sum(1),
				ppoints(input.df$`read count normalised adjusted`),
				0.5
			),
			size = input.df$sigma_raw %>% unique %>% exp %>% `^` (-1),
			mu = input.df$lambda %>% unique %>% exp
		)

	input.df %>%	mutate(predicted_NB = predicted_NB)
}

plot_trends = function(input.df, symbols) {
	input.df %>%
		filter(transcript %in% symbols) %>%
		ggplot(
			aes(
				x = `predicted_NB` + 1,
				y = `read count normalised adjusted` + 1,
				label = transcript,
				color = `regression`
			)
		) +
		geom_abline(intercept = 0,
								slope = 1,
								color = "grey") +
		geom_jitter() +
		geom_smooth(method = "lm", formula = y ~ x + I(x ^ 2)) +
		facet_wrap( ~ transcript, scales = "free")  +
		#expand_limits(y = 1, x = 1) +
		scale_y_log10() +
		scale_x_log10()

}

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

foreach(r = 1:6) %do% {
	res =
		TCGA_tbl %>%

		create_input_df(r)  %>%

		ppc_seq(
			significance_column = PValue,
			do_check_column = do_check,
			value_column = `read count`,
			percent_false_positive_genes = "5%",
			sample_column = sample,
			gene_column = transcript,
			pass_fit = T,
			tol_rel_obj = 0.01,
			just_discovery = T, full_bayes = F,
			cores = 10,
			additional_parameters_to_save = c("sigma_raw", "alpha")
		)

	saveRDS(res, file=sprintf("dev/draw_qq_TCGA_no_covariates_%s.rds", r))

}

# Analysis of confounders

TCGA_tbl.MDS =
	TCGA_tbl %>%
	create_ttBulk(sample_column = sample, transcript_column = transcript, counts_column = `read count`) %>%
	normalise_counts() %>%
	reduce_dimensions(value_column = `read count normalised`, method = "MDS", components = 1:10)

TCGA_tbl.MDS %>%
	select(sample,  `CAPRA-S`, contains("Dimension")) %>%
	distinct() %>%
	filter(`CAPRA-S` %>% is.na %>% `!`) %>%
	mutate(risk = `CAPRA-S` <= 3) %>%
	mutate(`CAPRA-S` = `CAPRA-S` %>% as.integer %>% as.factor)  %>%
	separate(sample, c("data base", "laboratory", "patient"), sep="-", remove = F) %>%
	inner_join(
		(.) %>% distinct(sample, laboratory) %>% count(laboratory) %>% filter(n >= 8) %>% select(-n)
	) %>%
	select(contains("Dimension"), everything()) %>%
	distinct() %>%
	GGally::ggpairs(columns = 1:6, ggplot2::aes(colour=`laboratory`))

# Model with covariates

foreach(r = 1:6) %do% {
	res =
		TCGA_tbl %>%

		create_input_df(r)  %>%

		ppc_seq(
			~ risk + laboratory,
			significance_column = PValue,
			do_check_column = do_check,
			value_column = `read count`,
			percent_false_positive_genes = "5%",
			sample_column = sample,
			gene_column = transcript,
			pass_fit = T,
			tol_rel_obj = 0.01,
			just_discovery = T, full_bayes = F,
			cores = 10,
			additional_parameters_to_save = c("sigma_raw", "alpha")
		)

	saveRDS(res, file=sprintf("dev/draw_qq_TCGA_with_covariates_%s.rds", r))

}



get_regression_adjusted = function(){
	x = readRDS("dev/draw_qq_TCGA_with_covariates_1.rds")


X =
	TCGA_tbl %>%
	create_input_df(1) %>%
	format_input(
		~ risk + laboratory,
		sample_column = sample,
		gene_column = transcript,
		value_column = `read count`,
		do_check_column = do_check,
		significance_column= PValue,
		how_many_negative_controls = 500
	) %>%
	create_design_matrix(~ risk + laboratory, sample)


alpha =
	x %>%
	attr("fit") %>%
	rstan::summary() %$% summary %>%
	as_tibble(rownames="par") %>%
	filter(grepl("alpha", par)) %>%
	separate(par, c("par", "C", "G"), sep="\\[|,|\\]", extra = "drop") %>%
	mutate(C = C %>% as.integer, G = G %>% as.integer) %>%
	arrange(G, C) %>%
	left_join(x %>% distinct(G, transcript), by="G") %>%
	select(transcript, G, C, mean) %>%
	spread(C, mean) %>% filter(transcript %>% is.na %>% `!`) %>%
	select(-G) %>% ttBulk::as_matrix(rownames="transcript")

X[,-1] %*% t(alpha[,-1]) %>%
	as_tibble(rownames="S") %>%
	gather(transcript, `log adjustment`, -S) %>%
	mutate(S = S %>% as.integer) %>%
	left_join(
		x %>% distinct(transcript, S, sample, G, `exposure rate`, `read count`)
	) %>%

	# normalise
	mutate(`read count normalised` = `read count` / exp(`exposure rate`))	%>%

	# Adjust
	mutate(`read count normalised adjusted` = `read count normalised` / exp(`log adjustment`))	%>%

	# Add lambda
	left_join(
		alpha %>% as_tibble(rownames="transcript") %>% select(1:2) %>% rename(lambda = `1`)
	) %>%

	# Add sigma
	left_join(
		x %>%
			attr("fit") %>%
			rstan::summary() %$% summary %>%
			as_tibble(rownames="par") %>%
			filter(grepl("sigma_raw", par)) %>%
			separate(par, c("par", "G"), sep="\\[|,|\\]", extra = "drop") %>%
			mutate( G = G %>% as.integer) %>%
			rename(sigma_raw = mean) %>%
			select(G, sigma_raw)
	) %>%

	filter(`read count normalised adjusted` > 0) %>%
	group_by(transcript) %>%
	do({
		mdf = (.) %>%
			get_NB_qq_values %>%
			mutate(a = `read count normalised adjusted` %>% `+` (1) %>% log,
						 b = predicted_NB %>% `+` (1) %>% log)

		mdf %>%
			mutate(regression = lm(a ~  b + I(b ^ 2), data = (.))$coefficients[3])
	}) %>%
	ungroup()

}

get_regression = function(){
	x = readRDS("dev/draw_qq_TCGA_no_covariates_1.rds")



	alpha =
		x %>%
		attr("fit") %>%
		rstan::summary() %$% summary %>%
		as_tibble(rownames="par") %>%
		filter(grepl("alpha", par)) %>%
		separate(par, c("par", "C", "G"), sep="\\[|,|\\]", extra = "drop") %>%
		mutate(C = C %>% as.integer, G = G %>% as.integer) %>%
		arrange(G, C) %>%
		left_join(x %>% distinct(G, transcript), by="G") %>%
		select(transcript, G, C, mean) %>%
		spread(C, mean) %>% filter(transcript %>% is.na %>% `!`) %>%
		select(-G) %>% ttBulk::as_matrix(rownames="transcript")


		x %>%
			distinct(transcript, S, sample, G, `exposure rate`, `read count`) %>%

		# normalise
		mutate(`read count normalised` = `read count` / exp(`exposure rate`))	%>%

		# Adjust
		mutate(`read count normalised adjusted` = `read count normalised` )	%>%

		# Add lambda
		left_join(
			alpha %>% as_tibble(rownames="transcript") %>% select(1:2) %>% rename(lambda = `1`)
		) %>%

		# Add sigma
		left_join(
			x %>%
				attr("fit") %>%
				rstan::summary() %$% summary %>%
				as_tibble(rownames="par") %>%
				filter(grepl("sigma_raw", par)) %>%
				separate(par, c("par", "G"), sep="\\[|,|\\]", extra = "drop") %>%
				mutate( G = G %>% as.integer) %>%
				rename(sigma_raw = mean) %>%
				select(G, sigma_raw)
		) %>%

		filter(`read count normalised adjusted` > 0) %>%
		group_by(transcript) %>%
		do({
			mdf = (.) %>%
				get_NB_qq_values %>%
				mutate(a = `read count normalised adjusted` %>% `+` (1) %>% log,
							 b = predicted_NB %>% `+` (1) %>% log)

			mdf %>%
				mutate(regression = lm(a ~  b + I(b ^ 2), data = (.))$coefficients[3])
		}) %>%
		ungroup()

}

get_regression_adjusted_theoretical  = function(){

	yes %>%
		distinct(transcript, lambda, sigma_raw, sample) %>%
		group_by(transcript) %>%
		do(
			(.) %>% mutate(`read count normalised adjusted` = rnbinom((.) %>% nrow, mu = (.) %>% slice(1) %>% pull(lambda) %>% exp, size = (.) %>% slice(1) %>% pull(sigma_raw) %>% exp %>% `^` (-1) ))
		) %>%
		ungroup() %>%

	filter(`read count normalised adjusted` > 0) %>%
		group_by(transcript) %>%
		do({
			mdf = (.) %>%
				get_NB_qq_values %>%
				mutate(a = `read count normalised adjusted` %>% `+` (1) %>% log,
							 b = predicted_NB %>% `+` (1) %>% log)

			mdf %>%
				mutate(regression = lm(a ~  b + I(b ^ 2), data = (.))$coefficients[3])
		}) %>%
		ungroup()

}


yes = get_regression_adjusted()
no = get_regression()
yes_theoretical = get_regression_adjusted_theoretical()
lowly_transcribed =
	TCGA_tbl %>%
	ttBulk::normalise_counts(sample_column = sample, transcript_column=transcript, counts_column = `read count`, action="get") %>%
	distinct(transcript, `filtered_out_low_counts`)

(
	no %>%
		drop_na() %>%
		# Prepare for plot
		mutate(transcript = factor(transcript, unique(transcript))) %>%
		distinct(transcript, regression, lambda) %>%
		arrange(`regression` %>% desc) %>%
		mutate(tt = 1:n(), adjusted = F, theoretical = F) %>%

		bind_rows({

			my_n = (.) %>% nrow
			yes %>%

				# Prepare for plot
				mutate(transcript = factor(transcript, unique(transcript))) %>%
				distinct(transcript, regression, lambda) %>%
				drop_na() %>%
				sample_n( my_n) %>%
				arrange(`regression` %>% desc) %>%
				mutate(tt = 1:n(), adjusted = T, theoretical = F)
		}) %>%

		bind_rows({

			my_n = (.) %>% nrow
			yes_theoretical %>%

				# Prepare for plot
				mutate(transcript = factor(transcript, unique(transcript))) %>%
				distinct(transcript, regression, lambda) %>%
				drop_na() %>%
				ifelse_pipe( (.) %>% nrow %>% `>` (my_n), ~ .x %>% sample_n( my_n)) %>%
				arrange(`regression` %>% desc) %>%
				mutate(tt = 1:n(), adjusted = T, theoretical = T)
		}) %>%

		# Filter lowly transcribed
		left_join( lowly_transcribed ) %>%
		filter(!filtered_out_low_counts) %>%

		ggplot( aes(
			x = `tt`,
			y = `regression` ,
			label = transcript,
			#alpha=lambda,
			color = interaction(adjusted, theoretical),
			group=adjusted,
			#size=lambda,
			lmbda = lambda
		)) +
		geom_point() +
		scale_color_brewer(palette="Set1") +
		theme(
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank()
		) +
		my_theme +
		ylab("Quadratic regression curvature coefficient") +
		xlab("Gene rank")
)

%>% plotly::ggplotly()


yes %>% plot_trends("WDR5")



no %>% plot_trends("MRPS6")

