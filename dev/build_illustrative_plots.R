# Build illustrative plots

library(tidyverse)
library(ppcSeq)

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



# res =
# 	ppcSeq::counts %>%
# 	mutate(is_significant = FDR < "0.2") %>%
# 	ppc_seq(
# 		formula = ~ Label,
# 		significance_column = PValue,
# 		do_check_column = is_significant,
# 		value_column = value,
# 		save_generated_quantities = T,
# 		percent_false_positive_genes = "5%",
# 		approximate_posterior_inference = F,
# 		approximate_posterior_analysis = F,
# 		pass_fit = T,
# 		additional_parameters_to_save = c("sigma_raw", "lambda_mu", "exposure_rate", "sigma_slope", "sigma_intercept", "intercept", "lambda_log_param")
# 	)
#
# res %>% saveRDS("dev/res_for_build-illustrative_plots.rds", compress = "gzip")

res %>% saveRDS("dev/res_for_build-illustrative_plots.rds")

# MU-SIGMA plot
mu_sigma =
	res %>%
	attr("fit 1") %>%
	tidybayes::spread_draws(sigma_raw[G], intercept[G])

(
	mu_sigma %>%
		ggplot(aes(x = intercept, y = sigma_raw, group=G)) +
		stat_ellipse(geom = "polygon", fill="grey", alpha = 0.075, color = NA) +
		# stat_ellipse(geom = "polygon", fill="grey", alpha = 0.075, color = NA, level = 0.8) +
		# stat_ellipse(geom = "polygon", fill="grey", alpha = 0.075, color = NA, level = 0.6) +
		stat_ellipse(geom = "polygon", fill="#397FB9", alpha = 0.075, color = NA, level = 0.4) +
		geom_point(data = mu_sigma %>% tidybayes::median_qi(), aes(x = intercept, y = sigma_raw), shape=".", color = "red") +
		xlab("Log. of abundance mean") +
		ylab("- Log. of abundance overdispersion") +
		my_theme
) %>%
	ggsave(
		"dev/mu_sigma.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183/2
	)

# EXPOSURE plot
(
	res %>%
	attr("fit 1") %>%
	tidybayes::spread_draws(exposure_rate[S], intercept[G]) %>%
	filter(S==1) %>%
	mutate(`Estimated log abundance for a sample` = intercept + exposure_rate * 30) %>%
	filter(.draws< 100) %>%
	ggplot(aes(`Estimated log abundance for a sample`, group=.draw)) +
	geom_density(alpha = 0.01, color="black") +
	ylab("Density") +
	my_theme
)  %>%
	ggsave(
		"dev/exposure_rate.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183/2
	)

# REGRESSION plot
reg =
	res %>%
	attr("fit 1") %>%
	tidybayes::spread_draws(lambda_log_param[S,G])


(
	res %>%
		ungroup() %>%
		left_join(ppcSeq::counts) %>%
		#tidyBulk::scale_abundance(sample, symbol, value) %>%
		filter(symbol == "ART3") %>%
		mutate(Label = Label %>% as.factor %>% as.integer) %>%
		ggplot(aes(
			x = factor(Label),
			y = value + 1,
			fill = factor(Label)
		)) +
		geom_boxplot() +
		geom_line(
			data = reg %>%
				filter(G == 3) %>%
				left_join(
					ppcSeq::counts %>%
						mutate(S = sample %>% as.factor %>% as.integer) %>%
						distinct(S, Label) %>%
						mutate(Label = Label %>% as.factor %>% as.integer)
				),
			aes(x = Label, y = exp(lambda_log_param - 1.5), group = .draw), alpha = 0.1
		) +
		scale_fill_brewer(palette = "Set1") +
		scale_y_log10() +
		my_theme
) %>%
	ggsave(
		"dev/regression.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183/2
	)

# GENERATED QUANTITIES plot
(
	res %>%
	select(`sample wise data`) %>%
	unnest(`sample wise data`) %>%
	slice(1) %>%
	unnest(`generated quantities`) %>%
	ungroup() %>%
	filter(G == 203) %>%
	ggplot(aes(.value + 1)) + geom_density(fill = "grey") + scale_x_log10() +
	xlab("Logarithm of transcript abundance") +
	ylab("Dnesity") +
	my_theme
) %>%
	ggsave(
		"dev/outlier_density.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183/2
	)


