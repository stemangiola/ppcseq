library(tidyverse)

read_csv("dev/elapsed_data.csv") %>%
	filter(FDR<1) %>%
	ggplot(aes(x = draws, y=memory)) +
	geom_point(aes(color=fp_rate)) +
	geom_smooth(method = "lm") +
	facet_grid(approximate_posterior_inference~cores)

read_csv("dev/elapsed_data.csv") %>%
	filter(FDR<1) %>%
	ggplot(aes(x = draws, y=elapsed)) +
	geom_point(aes(color=approximate_posterior_analysis)) +
	geom_smooth(method = "lm") +
	facet_grid(approximate_posterior_inference ~cores)

read_csv("dev/elapsed_data.csv") %>%
	filter(FDR<1) %>%
	nest(data = -c(bayes)) %>%
	mutate(fit = map(.x=data,~lm(memory ~ draws, data = .x))) %>%
	{map(.$fit, summary)}
