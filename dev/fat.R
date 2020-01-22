library(tidyverse)
library(magrittr)
library(ppcSeq)

# res =
# 	ppcSeq::ppc_seq(
# 		dplyr::mutate(ppcSeq::counts,  is_significant = FDR < 0.05 ),
# 		formula = ~ Label,
# 		significance_column = PValue,
# 		do_check_column  = is_significant,
# 		value_column = value,
# 		percent_false_positive_genes = "0.5%",
# 		approximate_posterior_inference = F,
# 		approximate_posterior_analysis = F,
# 		how_many_negative_controls = 500,
# 		cores=30,
# 		pass_fit = T
# 	)
#
# save(res, file = "dev/fat_ppc.rda")

load("dev/fat_ppc.rda")

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
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

p1 =
	res %>% left_join(ppcSeq::counts %>% distinct(symbol, FDR)) %>%
	select(-`sample wise data`, - plot) %>%
	ungroup() %>%
	mutate(symbol = factor(symbol, unique(symbol))) %>%
	ggplot(aes(x=symbol, y=`ppc samples failed` )) + geom_point() + my_theme

legend_p2 =
	res %>%
	filter(`ppc samples failed` == 1) %>%
	pull(plot) %>%
	`[[` (1) %>%
	cowplot::get_legend()

p2 =
	res %>% filter(`ppc samples failed` == 1) %>%
	pull(plot) %>%
	map(~ .x + theme(
		legend.position="none",
		axis.title.x  = element_blank(),
		axis.title.y  = element_blank(),
		axis.text.x = element_text(angle = 90, hjust = 1, size = 5)
		)
	) %>%

	# eliminate x labels. Anonymous function
	{
		x = (.)
		t = theme(
			axis.title.x = element_blank(),
			axis.text.x = element_blank(),
			axis.ticks.x = element_blank()
		)
		x[[1]] = x[[1]] + t
		x[[2]] = x[[2]] + t
		x[[3]] = x[[3]] + t
		x
	} %>%
	cowplot::plot_grid(plotlist = ., align = "h", ncol = 3, axis="b" )


cowplot::plot_grid(p1, p2, align = "h", ncol = 1, axis="b" , rel_heights = c(1,4)) %>%
	ggsave(
		"dev/fat_ppc.pdf",
		plot = .,
		useDingbats=FALSE,
		units = c("mm"),
		width = 183 ,
		height = 347,
		limitsize = FALSE
	)

res %>% filter(`ppc samples failed` == 1) %>%
	pull(plot) %>% `[[` (1) + theme( legend.position = "bottom")

# How many outliers after first passage
res %>% filter(`tot deleterious outliers` == 0) %>% unnest(`sample wise data`) %>% rowwise() %>% filter(! between(value, .lower_2, .upper_2)) %>% ungroup()

# What is the n_eff I achieve on real data
res %>% attr("fit 2") %>%
	rstan::summary() %$% summary %>%
	as_tibble(rownames="par") %>%
	filter(grepl("counts_rng", par)) %>%
	ggplot(aes(n_eff)) + geom_histogram() +
	my_theme
