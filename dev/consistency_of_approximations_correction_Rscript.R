library(tidyverse)
library(magrittr)
library(ppcSeq)
# plan(multicore)
library(foreach)

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
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

# Produce input
FDR_threshold = 0.2
percent_false_positive_genes = "5%"
additional_parameters_to_save = c("intercept", "alpha", "sigma_raw")
adj_prob_theshold_2 = 0.00476190476190476

make_df_plot = function(.data){
	.data %>%
		ungroup() %>%
		mutate(
			.upper_diff = .upper_2.x - .upper_2.y,
			.lower_diff = .lower_2.x - .lower_2.y,
			mean_diff = mean_2.x - mean_2.y
		) %>%
		# mutate(
		# 	.upper_diff = .upper_diff / mean_2.x,
		# 	.lower_diff = .lower_diff  / mean_2.x,
		# 	mean_diff = mean_diff  / mean_2.x
		# ) %>%
		distinct(sample, symbol, .upper_diff, .lower_diff, mean_diff, mean_2.x, intercept, sigma_raw, title) %>%
		drop_na() %>%
		gather(Comparison, Difference, c(.upper_diff, .lower_diff, mean_diff)) %>%
		#mutate(Difference = ifelse(Difference == 0, 1, Difference)) %>%
		mutate(sign = Difference > 0) %>%
		mutate(title = factor(title,  levels=c("full", "Approximate posterior inference", "Approximate posterior analysis", "Approximate posterior inference and analysis"))) %>%
		mutate(Comparison = factor(Comparison,  levels=c("mean_diff", ".lower_diff", ".upper_diff")))

}

wrapper = function(adj_prob_theshold_2, do_correct_approx = F){

	res_0 =
		ppcSeq::counts %>%
		mutate(is_significant = FDR < FDR_threshold) %>%
		ppc_seq(
			formula = ~ Label,
			significance_column = PValue,
			do_check_column = is_significant,
			value_column = value,
			percent_false_positive_genes = percent_false_positive_genes,
			approximate_posterior_inference = F,
			approximate_posterior_analysis = F,
			pass_fit = T,
			additional_parameters_to_save = additional_parameters_to_save,
			save_generated_quantities = T,
			adj_prob_theshold_2 = adj_prob_theshold_2, do_correct_approx = do_correct_approx
		)

	fit_to_lambda_sigma = function(fit){
		fit %>% attr("fit 2") %>%
			rstan::summary(pars=c("intercept", "sigma_raw")) %>%
			as.data.frame() %>% as_tibble(rownames = "par") %>%
			separate(par, c("par", "G"), sep="\\[|\\]") %>%
			select(par, G, summary.mean) %>% mutate(G = G %>% as.integer) %>%
			spread(par, summary.mean)
	}

	input_0 =
		res_0 %>%
		ungroup %>%
		select(symbol, `sample wise data`) %>%
		unnest(cols = `sample wise data`) %>%
		select(symbol, sample, mean_2, .lower_2, .upper_2 , `generated quantities`) %>%

		# Filter first draw
		mutate(`generated quantities` = map(`generated quantities`, ~ .x %>% filter(`.draw` %in% c(1, 10, 30, 50, 200) & .chain == 1))) %>%

		# Attach annotation
		left_join( ppcSeq::counts %>% select(-value) %>% distinct()	) %>%

		# unpack
		unnest(cols =  `generated quantities`) %>%
		unite(sample_draw, c(sample, .draw), sep="_", remove = F) %>%
		select(-`.chain`, -`.iteration`, -`.draw`) %>%
		rename(value = `.value`) %>%


		# Add negative controls
		bind_rows(
			(.) %>% distinct(sample, sample_draw) %>%
				left_join(
					ppcSeq::counts %>%
						inner_join(
							(.) %>%
								arrange(PValue) %>%
								distinct(symbol) %>%
								tail(n=2000)
						)
				)

		) %>%
		select(-sample) %>%
		rename(sample = sample_draw) %>%
		mutate(value = value %>% as.integer) %>%
		mutate(is_significant = FDR < FDR_threshold)

	# saveRDS(input_0, file="dev/input_0.rds")
	#
	# input_0 = readRDS("dev/input_0.rds")

	res_1 =
		input_0 %>%
		ppc_seq(
			formula = ~ Label,
			significance_column = PValue,
			do_check_column = is_significant,
			value_column = value,
			approximate_posterior_inference = F,
			approximate_posterior_analysis = F,
			percent_false_positive_genes = percent_false_positive_genes,
			cores = 20,
			seed = 654321,
			pass_fit = T,
			additional_parameters_to_save = additional_parameters_to_save,
			adj_prob_theshold_2 = adj_prob_theshold_2, do_correct_approx = do_correct_approx
		)

	res_2 =
		input_0 %>%
		ppc_seq(
			formula = ~ Label,
			significance_column = PValue,
			do_check_column = is_significant,
			value_column = value,
			approximate_posterior_inference = T,
			approximate_posterior_analysis = T,
			percent_false_positive_genes = percent_false_positive_genes,
			cores = 20,
			additional_parameters_to_save = additional_parameters_to_save,
			adj_prob_theshold_2 = adj_prob_theshold_2, pass_fit = T, do_correct_approx = do_correct_approx
		)

	res_3 =
		input_0 %>%
		ppc_seq(
			formula = ~ Label,
			significance_column = PValue,
			do_check_column = is_significant,
			value_column = value,
			approximate_posterior_inference = F,
			approximate_posterior_analysis = T,
			percent_false_positive_genes = percent_false_positive_genes,
			cores = 20,
			adj_prob_theshold_2 = adj_prob_theshold_2,
			additional_parameters_to_save = additional_parameters_to_save,
			pass_fit = T, do_correct_approx = do_correct_approx
		)

	res_4 =
		input_0 %>%
		ppc_seq(
			formula = ~ Label,
			significance_column = PValue,
			do_check_column = is_significant,
			value_column = value,
			approximate_posterior_inference = T,
			approximate_posterior_analysis = F,
			percent_false_positive_genes = percent_false_positive_genes,
			cores = 20,
			adj_prob_theshold_2 = adj_prob_theshold_2,
			additional_parameters_to_save = additional_parameters_to_save,
			pass_fit = T, do_correct_approx = do_correct_approx
		)
	res_1_parsed =
		res_1 %>%
		unnest(col = "sample wise data") %>%
		distinct(symbol, G, sample, value, mean_2, .lower_2, .upper_2)%>%
		left_join(
			fit_to_lambda_sigma(res_1)
		)

	res_2_parsed =
		res_2 %>%
		unnest(col = "sample wise data") %>%
		distinct(symbol, G, sample, value, mean_2, .lower_2, .upper_2) %>%
		left_join(
			fit_to_lambda_sigma(res_2)
		)

	res_3_parsed =
		res_3 %>%
		unnest(col = "sample wise data") %>%
		distinct(symbol, G, sample, value, mean_2, .lower_2, .upper_2) %>%
		left_join(
			fit_to_lambda_sigma(res_3)
		)

	res_4_parsed =
		res_4 %>%
		unnest(col = "sample wise data") %>%
		distinct(symbol, G, sample, value, mean_2, .lower_2, .upper_2) %>%
		left_join(
			fit_to_lambda_sigma(res_4)
		)

	res_parsed_1 =
		input_0 %>%
		left_join(res_1_parsed,	by=c("sample", "symbol")) %>%
		mutate(title = "full") %>%
		bind_rows(
			input_0 %>%
				left_join(res_2_parsed,	by=c("sample", "symbol")) %>%
				mutate(title = "Approximate posterior inference and analysis")
		) %>%
		bind_rows(
			input_0 %>%
				left_join(res_3_parsed,	by=c("sample", "symbol")) %>%
				mutate(title = "Approximate posterior analysis")
		) %>%
		bind_rows(
			input_0 %>%
				left_join(res_4_parsed,	by=c("sample", "symbol")) %>%
				mutate(title = "Approximate posterior inference")
		)

	res_parsed_2 =
		res_1_parsed %>%
		left_join(res_2_parsed,	by=c("sample", "symbol")) %>%
		mutate(title = "Approximate posterior inference and analysis") %>%
		bind_rows(
			res_1_parsed %>%
				left_join(res_3_parsed,	by=c("sample", "symbol")) %>%
				mutate(title = "Approximate posterior analysis")
		) %>%
		bind_rows(
			res_1_parsed %>%
				left_join(res_4_parsed,	by=c("sample", "symbol")) %>%
				mutate(title = "Approximate posterior inference")
		)


	df_plot_1 =
		res_parsed_1 %>%
		make_df_plot %>%
		mutate(p = 1)

	df_plot_2 =
		res_parsed_2 %>%
		make_df_plot %>%
		mutate(p = 2)


	df_plot = bind_rows(df_plot_1, df_plot_2) %>%
		filter(Difference != 0) %>%
		mutate(adj_prob_theshold_2 = adj_prob_theshold_2)
}

df_plot_all =
	foreach(p = c(0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001), .combine=bind_rows ) %do% {
		w = wrapper(p, do_correct_approx = T)
		save(w, file=sprintf("dev/w_corrected_%s_exp_model.rda", p))
		w
	}
