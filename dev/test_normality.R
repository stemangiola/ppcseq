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
		#axis.text.x = element_text(angle = 90, hjust = 1, size = 5),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		aspect.ratio=1
	)

# res =
# 	ppcSeq::ppc_seq(
# 		dplyr::mutate(ppcSeq::counts,  is_significant = FDR < 0.05 ),
# 		formula = ~ Label,
# 		significance_column = PValue,
# 		do_check_column  = is_significant,
# 		value_column = value,
# 		percent_false_positive_genes = "1%",
# 		approximate_posterior_inference = F,
# 		approximate_posterior_analysis = F,
# 		how_many_negative_controls = 500,
# 		cores=30,
# 		pass_fit = T, additional_parameters_to_save = c("intercept", "sigma_raw")
# 	)
#
# res %>% saveRDS("dev/fat_ppc_norality_test.rds")

res = readRDS("dev/fat_ppc_norality_test.rds")

res_tidy =
	res %>% attr("fit 2") %>%
	tidybayes::gather_draws(intercept[G], exposure_rate[S], sigma_raw[G])

p =
	res_tidy %>%
	ungroup %>%
	nest(data = -c(G, .variable, S)) %>%
	mutate(
		`P-value` = future_map(data, ~ shapiro.test(sample(.x$.value, size=5000) %>% scale) %>% unlist %>% `[` (2) %>% as.numeric )
	) %>%
	unnest(cols = c(`P-value`)) %>%
	ggplot(aes(`P-value`)) +
	geom_density() +
	geom_vline(xintercept = 0.05) +
	scale_x_log10() +
	ylab("Density") +
	my_theme

p %>%
	ggsave(
		"dev/fat_normality_test.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 89
	)
