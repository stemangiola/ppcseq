library(tidyverse)
library(magrittr)
library(ppcSeq)
library(furrr)
plan(multicore)

FDR_threshold = 0.2

res_1 =
	ppcSeq::counts %>%
	mutate(is_significant = FDR < FDR_threshold) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column = is_significant,
		value_column = value,
		save_generated_quantities = T,
		percent_false_positive_genes = "5%",
		full_bayes = T
	)

input_2 =
	res_1 %>%
	select(symbol, `sample wise data`) %>%
	unnest %>%
	select(symbol, sample, `generated quantities`) %>%
	unnest %>%
	filter(`.draw` == 1) %>%
	select(-`.chain`, -`.iteration`, -`.draw`) %>%
	rename(value = `.value`) %>%
	left_join( ppcSeq::counts %>% select(-value) %>% distinct()	) %>%

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
		mutate(`data source` = list(input_2)) %>%
		mutate(
			`false positive predicted` =
				future_map2(fp, `data source`, ~
						.y %>%
						ppc_seq(
							formula = ~ Label,
							significance_column = PValue,
							do_check_column = is_significant,
							value_column = value, full_bayes = F, percent_false_positive_genes = sprintf("%s%%", .x)
						) %>%
						filter( `tot deleterious outliers`>0) %>%
						nrow %>%
						divide_by( input_2 %>% filter(is_significant) %>% distinct(symbol) %>% nrow )
				 	)
		)

(es %>% mutate(fp = fp %>% divide_by(10)) %>% unnest(`false positive predicted`) %>% ggplot(aes(x=fp, y= `false positive predicted`)) + geom_jitter() + xlab("false positive") + my_theme) %>% ggsave(plot=., "dev/false_positive_study.pdf", device="pdf")

save(list=c("es", "input_2", "res_1"), file="dev/false_positive_study_3_runs.RData")


res_2 %>% unnest(`sample wise data`)
