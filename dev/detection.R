library(tidyverse)
library(magrittr)
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

FDR_threshold = 0.05

res =
	ppcSeq::counts %>%
	mutate(is_significant = FDR < FDR_threshold) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column  = is_significant,
		value_column = value,
		percent_false_positive_genes = "1%",
		tol_rel_obj = 0.01,
		approximate_posterior_inference = T,
		approximate_posterior_analysis = T
	)

(
	es %>%
		mutate(fp = fp %>% divide_by(10)) %>%
		unnest(`false positive predicted`) %>%
		ggplot(aes(x=fp, y= `false positive predicted`)) +
		geom_jitter(width = 0.01) +
		geom_smooth(method = "lm") +
		xlab("false positive") +
		my_theme
) %>%
	ggsave(plot=., "dev/false_positive_study.pdf", device = cairo_pdf, width=89, height = 90, units = "mm" )

save(list=c("es", "input_2", "res_1"), file="dev/false_positive_study_3_runs.RData")


res_2 %>% unnest(`sample wise data`)
