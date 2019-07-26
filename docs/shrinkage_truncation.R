
library(tidyverse)
library(foreach)
library(doParallel)
library(rstan)
doParallel::registerDoParallel()



mod <- stan_model(model_code = 'data{int X; int x[X];} parameters {real mu; real sigma;} model {x ~ neg_binomial_2_log(mu, 1/exp(sigma)); mu ~ normal(0,3); sigma ~ normal(0,3); }' )

shrinkage =
	foreach(mu = seq(10, 1000, 10), .combine = bind_rows) %dopar%
	{
		foreach(sigma = seq(1, 100, 1), .combine = bind_rows) %do%	{

			x =
				rnbinom(1000, mu = mu, size = sigma) %>%
				enframe %>%
				filter(value <= qnbinom(0.975, mu=mu, size=sigma)) %>%
				filter(value >= qnbinom(0.025, mu=mu, size=sigma)) %>%
				pull(value)

			X = x %>% length

			fit = optimizing(mod, data=list(X = X, x = x), hessian = TRUE)
			tibble(
				mu = !!mu,
				sigma = !!sigma,
				sigma_pred = 1 / exp(fit$par[2]),
				mu_pred = exp(fit$par[1])
			)

		}
	} %>%
	mutate(mu_shinkage = mu_pred / mu, sigma_shinkage = sigma_pred / sigma)

shrinkage %>% filter(sigma_shinkage < 100) %>% ggplot(aes(mu, sigma)) + geom_raster(aes(fill = sigma_shinkage))
shrinkage %>% filter(sigma_shinkage < 100)  %>% summarise(median(sigma_shinkage))
shrinkage %>% ggplot(aes(mu, sigma)) + geom_raster(aes(fill = mu_shinkage))
shrinkage %>% summarise(median(mu_shinkage))
