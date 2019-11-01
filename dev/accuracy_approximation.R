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
	filter(
		approximate_posterior_inference == F &
			approximate_posterior_analysis == F &
			fp_rate == "1%" &
			genes == 164
	) %>%
	slice(1) %>%
	map(.$)



"/stornext/Bioinf/data/bioinf-data/Papenfuss_lab/projects/ppc-benchmarking/archive/v0.4.1/tmp_kzehrku.rds" %>%
	readRDS() %>%
	attributes %>%
	names

	ggplot(aes(x = draws, y=memory)) +
	geom_point(aes(color=fp_rate)) +
	geom_smooth(method = "lm") +
	facet_grid(approximate_posterior_inference + approximate_posterior_analysis ~cores) +
	my_theme
