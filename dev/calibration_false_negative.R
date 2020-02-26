library(tidyverse)
library(magrittr)
library(ppcSeq)
library(furrr)
plan(multiprocess)


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

FDR_threshold = 0.2
#
# res_1 =
# 	ppcSeq::counts %>%
# 	mutate(is_significant = FDR < FDR_threshold) %>%
# 	ppc_seq(
# 		formula = ~ Label,
# 		significance_column = PValue,
# 		do_check_column = is_significant,
# 		value_column = value,
# 		save_generated_quantities = T,
# 		percent_false_positive_genes = "5%",
# 		approximate_posterior_inference = F,
# 		approximate_posterior_analysis = F,
# 		additional_parameters_to_save = c("exposure_rate", "lambda_log_param",	"sigma_raw"),
# 		pass_fit = T,
# 		cores = 30
# 	)
#
# res_1 %>% saveRDS("dev/res_1_for_calibration_false_negative.rds")

res_1 = readRDS("dev/res_1_for_calibration_false_negative.rds")

outlier_prop = 1e-10
how_many_outliers = 0.2

# outlier_df =
	# res_1 %>%
	# attr("fit 2") %>%
	# tidybayes::spread_draws(lambda_log_param[S, G], sigma_raw[G], exposure_rate[S]) %>%
	# tidybayes::median_qi() %>%
	# ungroup()  %>%
	# nest(data = -c(S, G)) %>%
	# mutate(
	# 	CI = map(
	# 		data,
	#
	# 		# Calculate quantiles
	# 		~ qnbinom(
	# 			c(outlier_prop, 1-outlier_prop),
	# 			mu=exp(.x$lambda_log_param + .x$exposure_rate),
	# 			size = 1/exp(.x$sigma_raw)
	# 		) %>%
	# 			enframe(name=NULL) %>%
	# 			mutate(id = c("outlier_low", "outlier_high")) %>%
	# 			spread(id, value)
	# 	)
	# ) %>%
	# unnest(CI) %>%
	# select(S, G, contains("outlier"))
#
# outlier_df %>% saveRDS("dev/outlier_for_calibration_false_negative.rds")

outlier_df = readRDS("dev/outlier_for_calibration_false_negative.rds")
set.seed(534231)
input_2 =
	res_1 %>%
	ungroup %>%
	select(symbol, `sample wise data`) %>%
	unnest(cols = `sample wise data`) %>%
	select(S, G, symbol, sample, mean_2, .lower_2, .upper_2 , `generated quantities`) %>%

	# Filter first draw
	mutate(`generated quantities` = map(`generated quantities`, ~ .x %>% filter(`.draw` == 1))) %>%

	# unpack
	unnest(cols =  `generated quantities`) %>%
	select(-`.chain`, -`.iteration`, -`.draw`) %>%
	rename(value = `.value`) %>%
	left_join( ppcSeq::counts %>% select(-value) %>% distinct()	) %>%

	# Add outliers
	left_join( outlier_df	) %>%
	left_join(
		(.) %>%
			distinct(symbol) %>%
			mutate(is_symbol_outlier = sample(1:3,size = n(), prob = c(1-how_many_outliers, rep(how_many_outliers, 2)/2), replace = T ))
	)	%>%
	group_by(symbol) %>%
	mutate(is_outlier = ifelse(is_symbol_outlier & row_number() == sample(size=1, 1:n()), is_symbol_outlier, 1)) %>%
	ungroup() %>%

	# Outliers all up
	mutate(value = ifelse(is_outlier == 1, value, outlier_high)) %>%
	#mutate(value = recode(is_outlier, value, outlier_low, outlier_high)) %>%

	# Add negative controls
	bind_rows(
		ppcSeq::counts %>%
			inner_join(
				(.) %>%
					arrange(PValue) %>%
					distinct(symbol) %>%
					tail(n=2000)
			)
	) %>%
	mutate(value = value %>% as.integer) %>%
	mutate(is_significant = FDR < FDR_threshold)

es =
	expand.grid(
		fp = c(seq(0.2, 0.9, 0.1), seq(1, 10, 1)),
		run = 1:3
	) %>%
		as_tibble() %>%
		#slice(1) %>%   # <<<<<<<<< TESING
		mutate(`data source` = list(input_2)) %>%
		mutate(
			`inference` =
				map2(fp, `data source`, ~
						.y %>%
						ppc_seq(
							formula = ~ Label,
							significance_column = PValue,
							do_check_column = is_significant,
							value_column = value,
							percent_false_positive_genes = sprintf("%s%%", .x),
							cores = 30,
							approximate_posterior_inference = F, approximate_posterior_analysis = F
						)
				 	)
		)

es %>% saveRDS("dev/es_calibration_false_negative.rds")

# es = readRDS("dev/es_calibration_false_negative.rds")


stats =
	es %>%
	mutate(
		inference = map(inference, ~ .x %>% select(symbol, `ppc samples failed`)) # unnest(`sample wise data`) %>% select(symbol, sample, ppc))
	) %>%
	mutate(`data source` = map(`data source`, ~ .x  %>% filter(is_significant) %>% distinct(symbol, is_symbol_outlier))) %>%
	mutate(integrated = map2(inference, `data source`, ~ .x %>% left_join(.y))) %>%
	select(-`data source` , -  inference) %>%
	unnest(integrated)

# Where the data does not have outliers
fp_stat =
	stats %>%
	filter(is_symbol_outlier == 1) %>%
	mutate(false_positive =   `ppc samples failed` > 0 ) %>%
	mutate(true_negative =  `ppc samples failed` == 0) %>%
	group_by(fp, run) %>%
	summarise(sum(false_positive), sum(true_negative)) %>%
	mutate(`false positive predicted` = `sum(false_positive)` / (`sum(false_positive)` + `sum(true_negative)`))

p1 =
	fp_stat %>%
	ggplot(aes(x = fp/100, y = `false positive predicted` )) +
	geom_smooth(method = "lm", color = "#4b68b1") +
	geom_point() +
	#geom_abline(intercept = 0, slope = 1, linetype ="dashed") +
	xlab("False positive aimed") +
	ylab("False positive predicted") +
	my_theme

fn_stat =
	stats %>%
	filter(is_symbol_outlier > 1) %>%
	mutate(true_positive =   `ppc samples failed` > 0 ) %>%
	mutate(false_negative =  `ppc samples failed` == 0) %>%
	group_by(fp, run) %>%
	summarise(sum(true_positive), sum(false_negative)) %>%
	mutate(`false negative predicted` = `sum(false_negative)` / (`sum(false_negative)` + `sum(true_positive)`))

p2 =
	fn_stat %>%
	ggplot(aes(x = fp / 100, y = `false negative predicted` )) + geom_point()  +
	geom_abline(intercept = 0, slope = 1) +
	xlab("False positive aimed") +
	ylab("False negative predicted") +
	my_theme

# ROC
p3 =
	tibble(
	`false positive predicted` = fp_stat %>% pull(`false positive predicted`) %>% sort,
	`false negative predicted` = fn_stat %>% pull(`false negative predicted`) %>% sort,
	) %>%
		ggplot(aes(y = `false negative predicted`, x = `false positive predicted`)) +
		geom_line() +
		xlab("False positive predicted") +
		ylab("False negative predicted") +
		my_theme

# Compose plots
cowplot::plot_grid(plotlist = list(p1, p2, p3), align = "v", ncol = 3, axis="b", rel_widths = 1 ) %>%

	# Save plots
	ggsave(
		"dev/fp_fn_ROC.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183
	)


#
# es %>% slice(1) %>% unnest(`false positive predicted`)
# (
# 	es %>%
# 		mutate(fp = fp %>% divide_by(100)) %>%
# 		unnest(`false positive predicted`) %>%
# 		ggplot(aes(x=fp, y= `false positive predicted`)) +
# 		geom_jitter(width = 0.01) +
# 		geom_smooth(method = "lm") +
# 		xlab("false positive") +
# 		my_theme
# ) %>%
# 	ggsave(plot=., "dev/false_positive_study.pdf", device = cairo_pdf, width=89, height = 90, units = "mm" )
#
#
# save(list=c("es", "input_2", "res_1"), file="dev/false_positive_study_3_runs.RData")
#
#
# res_2 %>% unnest(`sample wise data`)
