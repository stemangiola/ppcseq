functions{

	int[] get_elements_per_shard(int lenth_v, int shards){

	// Returned integer(max_size, last_element_size)
	int tentative_size = lenth_v / shards;
	int tentative_remaining = lenth_v - (tentative_size * shards);
	int elements_per_shard = tentative_remaining > 0 ? tentative_size + 1 : tentative_size;
	int remaining =  (elements_per_shard * shards) - lenth_v;

	int length_obj[shards];

	for(s in 1:shards) {
		length_obj[s] =
			s != shards ?
			elements_per_shard :
			elements_per_shard - remaining;  // Actual elements in it for last object
	}

 	return length_obj;

}


		vector[] get_mu_sigma_vector_MPI(vector mus, vector sigmas, int shards){

		int elements_per_shard[shards] = get_elements_per_shard(rows(mus), shards); // Length of the returned object
		int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
		vector[size_MPI_obj * 2] v_MPI[shards] ; // Set values to -999 for the ones that are not filled

		int i = 0; // Index sweeping the vector

		for(s in 1:shards){

			// If last shard fill in
			if(s == shards) v_MPI[s] = rep_vector(-999.0, size_MPI_obj * 2);

			v_MPI[s, 1:elements_per_shard[s]] = mus[ (i + 1) : i + elements_per_shard[s] ];
			v_MPI[s, (elements_per_shard[s]+1):(elements_per_shard[s]+elements_per_shard[s])] = sigmas[ (i + 1) : i + elements_per_shard[s] ];

			i += elements_per_shard[s];
		}

		return v_MPI;
	}
}

data {

	// Reference matrix inference
  int<lower=0> N;
	int y[N];
	int s[N];
	int g[N];
	int S;
	int G;

}
transformed data{
	shards = 10
}
parameters {

  // Gene-wise properties of the data
  vector[G] intercept;
  real lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma;
  real lambda_skew;

  // Overdispersion
  vector[G] sigma_raw;
  real<upper=0> sigma_slope;
  real sigma_intercept;
  real<lower=0>sigma_sigma;

  // Exposure rate
  vector[S] exposure_rate;

}
model {

  lambda_mu ~ normal(5,2);
  lambda_sigma ~ normal(0,2);
	lambda_skew ~ normal(0,1);

  sigma_intercept ~ normal(0,2);
  sigma_slope ~ normal(0,2);
  sigma_sigma ~ normal(0,2);

  // Gene-wise properties of the data
  to_vector(intercept) ~ skew_normal(lambda_mu ,lambda_sigma, lambda_skew);

  // Exposure prior
  exposure_rate ~ normal(0,2);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

	target += sum(map_rect(
		lp_reduce_simple ,
		[sigma_intercept, sigma_slope]', // global parameters
		get_mu_sigma_vector_MPI(
			intercept[g] + exposure_rate[s],
			sigma_raw[g],
			shards
		),
		real_data,
		get_int_MPI( y, shards)
	));

}
generated quantities{
	vector[how_many_to_check] counts_rng[S];

	for(g in 1:how_many_to_check) for(s in 1:S)
	// Make the overdispersion bigger making sigma smaller. Because inferring on truncated data with naive NB underestimate overdispersion
		counts_rng[s,g] =	neg_binomial_2_log_rng(exposure_rate[s] + lambda_log_param[s,g],	sigma[g] * truncation_compensation);

}
