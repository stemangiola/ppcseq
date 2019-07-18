library(tidyverse)
library(magrittr)
library(ppcSeq)

FDR_threshold = 0.1

res_1 =
	ppcSeq::counts %>%
	mutate(is_significant = FDR < FDR_threshold) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column = is_significant,
		value_column = value,
		save_generated_quantities = T
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
		value_column = value
	)



