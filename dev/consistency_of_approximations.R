library(tidyverse)
library(magrittr)
library(ppcSeq)
plan(multicore)
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

wrapper = function(adj_prob_theshold_2){

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
		adj_prob_theshold_2 = adj_prob_theshold_2
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
		adj_prob_theshold_2 = adj_prob_theshold_2
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
		adj_prob_theshold_2 = adj_prob_theshold_2, pass_fit = T
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
		pass_fit = T
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
		pass_fit = T
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
		w = wrapper(p)
		save(w, file=sprintf("dev/w_%s.rda", p))
		w
	}

# Plot bias trends
df_plot_all =
	dir(path="dev/", pattern = "^w_", full.names = T) %>%
	map_dfr(~ {
		load(.x)
		w
	})

df_plot_all %>%
	filter(Difference != 0) %>%
	group_by(Comparison, adj_prob_theshold_2, title, p) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%
	summarise(m = Difference %>% mean) %>%
	ggplot(aes(y=m, x=adj_prob_theshold_2)) + geom_point() + facet_wrap(title + p ~ Comparison, scales = "free_y")

# Plot biases
df_plot_all%>%
	filter(Difference != 0) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%
	left_join(
		ppcSeq::counts %>% group_by( symbol) %>% summarise(value = value %>% `+` (1) %>% log %>% mean, PValue = PValue %>% log %>% mean)
	) %>%
	#sample_frac(0.05) %>%
	{
		my_df = (.)
		ggplot(my_df, aes(x = interaction(sample, symbol), y=Difference, group=sign, sample = sample, symbol=symbol, color = value)) +

			geom_point(alpha=0.5, size=0.1) +
			scale_color_viridis(option="magma") +
			# geom_hline(
			# 	data = my_df  %>% group_by(Comparison, adj_prob_theshold_2, sign) %>% summarise(Difference %>% mean),
			# 	aes(yintercept=`Difference %>% mean`, color=sign),
			# 	size=1
			# ) +
			# geom_hline(
			# 	data =my_df %>% group_by(Comparison, adj_prob_theshold_2) %>% summarise(Difference %>% median),
			# 	aes(yintercept=`Difference %>% median`),
			# 	color = "yellow",
			# 	size=1
			# ) +
			facet_grid(Comparison ~ adj_prob_theshold_2, scales = "free_y") +
			theme(axis.line = element_line(),
						axis.title.x=element_blank(),
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank(),
						#legend.position = "none",
						text = element_text(size=12),
						strip.background = element_blank(),
						axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
			)
	}

# Plot lambda error
df_plot_all%>%
	filter(Difference != 0) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%
	left_join(
		ppcSeq::counts %>% group_by( symbol) %>% summarise(value = value %>% `+` (1) %>% log %>% mean, PValue = PValue %>% log %>% mean)
	) %>%
	ggplot( aes(x = mean_2.x, y=Difference, group=sign, sample = sample, symbol=symbol, color = value)) +

	geom_point(alpha=0.5, size=0.1) 	+		facet_grid(Comparison ~ adj_prob_theshold_2, scales = "free_y") +
	theme(axis.line = element_line(),
				legend.position = "none",
				text = element_text(size=12),
				strip.background = element_blank(),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)



# Scatter plot intercept sigma colored by error - UPPER diff
df_plot_all%>%
	filter(Difference != 0) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%

	# Correct bug
	select(-intercept, -sigma_raw) %>%
	left_join(
		df_plot_all %>% filter(p==1) %>% distinct(title, symbol, intercept, sigma_raw, adj_prob_theshold_2) %>% distinct()
	) %>%

	filter(Comparison == ".upper_diff") %>%
	sample_frac(0.5) %>%
	ggplot( aes(x = intercept, y=sigma_raw, group=sign, sample = sample, symbol=symbol, color = Difference)) +

	geom_jitter() 	+		facet_grid(Comparison ~ adj_prob_theshold_2, scales = "free_y") +
	scale_color_viridis(option="magma", trans = "log") +
	theme(axis.line = element_line(),
				#legend.position = "none",
				text = element_text(size=12),
				strip.background = element_blank(),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

# Scatter plot intercept sigma colored by error - LOWER diff
df_plot_all%>%
	filter(Difference != 0) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%

	# Correct bug
	select(-intercept, -sigma_raw) %>%
	left_join(
		df_plot_all %>% filter(p==1) %>% distinct(title, symbol, intercept, sigma_raw, adj_prob_theshold_2) %>% distinct()
	) %>%

	filter(Comparison == ".lower_diff") %>%
	sample_frac(0.5) %>%
	ggplot( aes(x = intercept, y=sigma_raw, group=sign, sample = sample, symbol=symbol, color = -Difference)) +

	geom_jitter() 	+		facet_grid(Comparison ~ adj_prob_theshold_2, scales = "free_y") +
	scale_color_gradient( trans = "log") +
	theme(axis.line = element_line(),
				#legend.position = "none",
				text = element_text(size=12),
				strip.background = element_blank(),
				axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

# Build linear model - UPPER diff
df_plot_all%>%
	filter(Difference != 0) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%

	# Correct bug
	select(-intercept, -sigma_raw) %>%
	left_join(
		df_plot_all %>% filter(p==1) %>% distinct(title, symbol, intercept, sigma_raw, adj_prob_theshold_2) %>% distinct()
	) %>%
	filter(Comparison == ".upper_diff") %>%
	group_by(symbol, intercept, sigma_raw, adj_prob_theshold_2) %>%
	summarise(Difference = Difference %>% median) %>%
	ungroup() %>%
	lm(log(Difference) ~ intercept + sigma_raw + log(adj_prob_theshold_2), data = .) %>%
	{
		lm_approx_bias_upper = (.)
		save(lm_approx_bias_upper, file = "dev/lm_approx_bias_upper.rda")
		(.) %>% summary
	}

# Build linear model - LOWER diff
df_plot_all%>%
	filter(Difference != 0) %>%
	filter(p==2 & title=="Approximate posterior analysis") %>%

	# Correct bug
	select(-intercept, -sigma_raw) %>%
	left_join(
		df_plot_all %>% filter(p==1) %>% distinct(title, symbol, intercept, sigma_raw, adj_prob_theshold_2) %>% distinct()
	) %>%
	filter(Comparison == ".lower_diff") %>%
	group_by(symbol, intercept, sigma_raw, adj_prob_theshold_2) %>%
	summarise(Difference = Difference %>% median) %>%
	ungroup() %>%
	lm(log(-Difference) ~ intercept + sigma_raw + log(adj_prob_theshold_2), data = .) %>%
	{
		lm_approx_bias_upper = (.)
		save(lm_approx_bias_upper, file = "dev/lm_approx_bias_lower.rda")
		(.) %>% summary
	}

save(df_plot_all, file="dev/df_plot_all.rda")

save(list=c("res_1", "res_2", "res_3", "res_4"), file="dev/consistency_or_approximation_fits.rda")

# Scatter plot




# res_parsed %>%
# 	ggplot(aes(x=mean_2.x + 1, y=mean_2.y + 1)) +
# 	geom_errorbar(aes(ymin = .lower_2.y + 1, ymax = .upper_2.y + 1), alpha = 0.05) +
# 	geom_errorbarh(aes(xmin = .lower_2.x + 1, xmax = .upper_2.x + 1), alpha = 0.05) +
# 	geom_point(color="#1B9CFC") +
# 	geom_abline(color="#e74c3c", linetype="dashed") +
# 	scale_y_log10() +
# 	scale_x_log10() +
# 	facet_wrap(~title) +
# 	xlab("Credible interval (0.95) of the posterior distribution using the full MCMC sampling") +
# 	ylab("Credible interval (0.95)") +
# 	my_theme

# Plot differences

signed_log <- scales::trans_new("signed_log",
																transform=function(x) sign(x)*log(abs(x)),
																inverse=function(x) sign(x)*exp(abs(x)))
library(scales)


save(df_plot, file="dev/df_plot.rda")

df_plot%>%
	# mutate(Difference = Difference %>% abs) %>%
	sample_frac(0.05) %>%
	{
		ggplot((.), aes(x = interaction(sample, symbol), y=Difference, group=sign, sample = sample, symbol=symbol)) +

			geom_point(alpha=0.5, size=0.1) +
			geom_hline(
				data = df_plot  %>% group_by(p, title, Comparison, sign) %>% summarise(Difference %>% mean),
				aes(yintercept=`Difference %>% mean`, color=sign),
				size=1
			) +
			geom_hline(
				data =df_plot %>% group_by(p, title, Comparison) %>% summarise(Difference %>% median),
				aes(yintercept=`Difference %>% median`),
				color = "yellow",
				size=1
			) +
			facet_grid(Comparison ~ p + title, scales = "free_y") +
			theme(axis.line = element_line(),
						axis.title.x=element_blank(),
						axis.text.x=element_blank(),
						axis.ticks.x=element_blank(),
						legend.position = "none",
						text = element_text(size=12),
						strip.background = element_blank(),
						axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
					)
	} %>%
	ggsave(
		"dev/consistency_of_approximation.pdf",
		plot = .,
		device = "pdf",
		useDingbats=FALSE,
		width=183,
		height = 120,
		units = "mm"
	)


# Model after correction

df_plot %>% group_by(p, title, Comparison) %>% summarise(Difference %>% median) %>% filter(title=="Approximate posterior analysis" & p ==2)
