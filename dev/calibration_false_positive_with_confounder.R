library(tidyverse)
library(magrittr)
library(ppcseq)
# library(furrr)
# plan(multisession, workers = 15)


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

# res_1 =
# 	ppcseq::counts %>%
# 	mutate(is_significant = FDR < FDR_threshold) %>%
# 	identify_outliers(
# 		formula = ~ Label + W,
# 		sample, symbol, value,
# 		.significance = PValue,
# 		.do_check = is_significant,
# 		save_generated_quantities = T,
# 		percent_false_positive_genes = 5,
# 		approximate_posterior_inference = F,
# 		approximate_posterior_analysis = F
# 	)
#
# # saveRDS(res_1, "dev/res_1_with_confounder.rds")
#
# input_2 =
# 	res_1 %>%
# 	ungroup %>%
# 	select(symbol, `sample wise data`) %>%
# 	unnest(cols = `sample wise data`) %>%
# 	select(symbol, sample, mean_2, .lower_2, .upper_2 , `generated quantities`) %>%
#
# 	# Filter first draw
# 	mutate(`generated quantities` = map(`generated quantities`, ~ .x %>% filter(`.draw` == 1))) %>%
#
# 	# unpack
# 	unnest(cols =  `generated quantities`) %>%
# 	select(-`.chain`, -`.iteration`, -`.draw`) %>%
# 	rename(value = `.value`) %>%
# 	left_join( ppcseq::counts %>% select(-value) %>% distinct()	) %>%
#
# 	# Add negative controls
# 	bind_rows(
# 		ppcseq::counts %>%
# 			inner_join(
# 				(.) %>%
# 					arrange(PValue) %>%
# 					distinct(symbol) %>%
# 					tail(n=2000)
# 			)
# 	) %>%
# 	mutate(value = value %>% as.integer) %>%
# 	mutate(is_significant = FDR < FDR_threshold)
#
# saveRDS(input_2, "dev/input_2_with_confounder.rds")

input_2 = readRDS("dev/input_2_with_confounder.rds")

es =
	expand.grid(
		fp = c(seq(0.2, 0.9, 0.1), seq(1, 10, 1)),
		run = 1:3
	) %>%
	as_tibble() %>%
	mutate(`data source` = list(input_2)) %>%
	mutate(
		`false positive predicted` =
			map2(fp, `data source`, ~
					 	.y %>%
					 	identify_outliers(
					 		formula = ~ Label + W,
					 		sample, symbol, value,
					 		.significance = PValue,
					 		.do_check = is_significant,
					 		percent_false_positive_genes = .x,
					 		approximate_posterior_analysis = T
					 	) %>%
					 	filter( `tot deleterious outliers`>0) %>%
					 	nrow %>%
					 	divide_by( input_2 %>% filter(is_significant) %>% distinct(symbol) %>% nrow )
			)
	)

#save(list=c("es", "input_2", "res_1"), file="dev/false_positive_study_3_runs_with_confounder.RData")
save(list=c("es", "input_2"), file="dev/false_positive_study_3_runs_with_confounder.RData")


(
	es %>%
		mutate(fp = fp %>% divide_by(100)) %>%
		unnest(`false positive predicted`) %>%
		ggplot(aes(x=fp, y= `false positive predicted`)) +
		geom_abline(linetype = "dotted", color="red") +
		#geom_jitter(width = 0.01) +
		geom_smooth(method = "lm") +
		geom_point() +
		ggpubr::stat_cor(method = "pearson") +
		xlab("false positive") +
		my_theme
) %>%
	ggsave(plot=., "dev/false_positive_study_with_confounder.pdf", device = cairo_pdf, width=89, height = 90, units = "mm" )


