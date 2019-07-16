library(tidyverse)
library(magrittr)
library(ppcSeq)

fit_in_out =
	ppcSeq::counts %>%
	mutate(is_significant = FDR < 0.05) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column = is_significant,
		value_column = value,
		save_generated_quantities = T,
		full_bayes = T,
		additional_parameters_to_save = c("exposure_rate","lambda_log_param",	"sigma"),
	)

fit = fit_in_out %$% fit
output = fit_in_out %$% output
input = fit_in_out %$% input


counts_rng =
	fit %>%
	tidybayes::gather_draws(counts_rng[S, G])

counts_rng %>% filter(G==1 & S==1) %>% pull(`.value`) %>% MASS::fitdistr("Negative Binomial")

counts_rng %$% fit %>%
	summary(c("lambda_log_param[1,1]", "exposure_rate[1]", "sigma[1]")) %$% summary %>%
	`[` (1:3) %>%
	enframe %>%
	spread(name, value) %>%
	mutate(exp(`1` + `2`))




output_2 =
	counts_rng %>%
	ungroup %>%
	filter(.draw ==1 ) %>%
	left_join(
		input %>% distinct(symbol, G, sample, S, Label, is_significant)
	) %>%
	left_join(
		output %>%
			distinct(symbol, `of which deleterious`, PValue)
	) %>%
	filter(`of which deleterious` == 0) %>%
	select(-`of which deleterious`,  -`.chain`, -`.iteration`,- `.draw`,    - `S` , -   `G`, -`.variable`) %>%
	mutate(`value` = `.value` %>% as.integer) %>%
	select(-`.value`) %>%
	ppc_seq(
		formula = ~ Label,
		significance_column = PValue,
		do_check_column = is_significant,
		value_column = value,
		full_bayes = T
	)


	summary("alpha_sub_1") %$%
	summary %>%
	as_tibble(rownames = ".variable") %>%
	separate(.variable, c(".variable", "G"), sep="[\\[,\\]]", extra="drop") %>%
	mutate(G = G %>% as.integer) %>%
	summarise(mean = mean %>% mean, sd = mean %>% sd)
