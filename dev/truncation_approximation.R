
library(tidyverse)
library(foreach)
library(doParallel)
library(rstan)
doParallel::registerDoParallel()



mo <- stan_model(model_code = 'data{int X; int x[X];} parameters {real mu; real sigma;} model {x ~ neg_binomial_2_log(mu, 1/exp(sigma)); mu ~ normal(0,3); sigma ~ normal(0,3); }' )

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

			tibble(
				mu = !!mu,
				sigma = !!sigma,
				shrink = optimizing(mo, data=list(X = X, x = x), hessian = TRUE) %>% { sigma / (1 / exp((.)$par[2]))	}
			)

		}
	}

library(viridis)
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
(shrinkage %>% ggplot(aes(mu, sigma)) + geom_raster(aes(fill = shrink)) + scale_fill_viridis(	option="magma") + my_theme) %>%
	ggsave(plot = .,
				 "dev/truncatin_approximation.pdf",
				 useDingbats=FALSE,
				 units = c("mm"),
				 width = 183
	)
