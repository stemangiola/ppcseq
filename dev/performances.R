library(tidyverse)

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

read_csv("dev/elapsed_data.csv") %>%
	filter(FDR<1) %>%
	ggplot(aes(x = draws, y=memory)) +
	geom_point(aes(color=fp_rate)) +
	geom_smooth(method = "lm") +
	facet_grid(approximate_posterior_inference + approximate_posterior_analysis ~cores) +
	my_theme

(read_csv("dev/elapsed_data.csv") %>%
	filter(approximate_posterior_analysis==F) %>%
	filter(FDR<1) %>%
	ggplot(aes(x = draws, y=elapsed)) +
	geom_point(aes(color=fp_rate)) +
	geom_smooth(method = "lm") +
	facet_grid(approximate_posterior_inference + approximate_posterior_analysis ~cores) +
	scale_colour_brewer(palette = "Set1") +
	my_theme) %>%
	ggsave(
		"dev/performances.pdf",
		plot = .,
		device = "pdf",
		useDingbats=FALSE,
		width=183,
		height = 120,
		units = "mm"
	)

read_csv("dev/elapsed_data.csv") %>%
	filter(FDR<1) %>%
	nest(data = -c(bayes)) %>%
	mutate(fit = map(.x=data,~lm(memory ~ draws, data = .x))) %>%
	{map(.$fit, summary)}
