library(tidyverse)
library(magrittr)
library(ppcSeq)

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

res_2 =
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
	mutate(is_significant = FDR < FDR_threshold) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column = is_significant,
		value_column = value, full_bayes = F, percent_false_positive_genes = "5%"
	)


res_2 %>% filter( `tot deleterious outliers`>0)
res_2 %>% unnest(`sample wise data`)
