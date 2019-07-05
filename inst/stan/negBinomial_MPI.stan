functions{

	real exp_gamma_meanSd_lpdf(vector x_log, real m_log, real s){

    // This function is the  probability of the log gamma function
    // in case you have data that is aleady in log form

    real m = exp(m_log);
    real v = m + square(m) * s;
    real a = square(m) / v;
    real b = m / v;

		vector[rows(x_log)] jacob = x_log; //jacobian
		real norm_constant = a * log(b) -lgamma(a);
		real a_minus_1 = a-1;
		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

	}

	real exp_gamma_meanSd_rng(real m_log, real s){
	  // This function takes care of the two prior choice
	  // without complicating too much the model itself

    real m = exp(m_log);
    real v = m + square(m) * s;
    real a = square(m) / v;
    real b = m / v;

	  return log(gamma_rng(a, b));
	}

	vector[] get_reference_parameters_MPI(int n_shards, int M, int[] G_per_shard, int[,] G_ind, vector lambda_log, vector sigma, vector exposure_rate){

		vector[2*M + cols(exposure_rate)] lambda_sigma_exposure_MPI[n_shards];

		for( i in 1:n_shards ) {

			int size_buffer = (M*2) - (G_per_shard[i]*2) ;
		  vector[ size_buffer] buffer = rep_vector(0.0,size_buffer);

			lambda_sigma_exposure_MPI[i] =
	  		append_row(
	  		  append_row(
		  		    append_row(
		  		    	lambda_log[G_ind[i, 1:G_per_shard[i]]],
		      		  sigma[G_ind[i, 1:G_per_shard[i]]]
		      		),
	      		buffer
	      	),
	      	exposure_rate
	      );
		}

		return(lambda_sigma_exposure_MPI);
	}

	vector lp_reduce( vector global_parameters , vector local_parameters , real[] real_data , int[] int_data ) {

		// Integer data unpack
	 	int M = int_data[1];
	 	int N = int_data[2];
	 	int S = int_data[3];
	 	int G_per_shard = int_data[4];
	 	int symbol_end[M+1] = int_data[(4+1):(4+1+M)];
	 	int sample_idx[N] = int_data[(4+1+M+1):(4+1+M+1+N-1)];
	 	int counts[N] = int_data[(4+1+M+1+N-1+1):size(int_data)];

		// Parameters unpack
	 	vector[G_per_shard] lambda_MPI = local_parameters[1:G_per_shard];
	 	vector[G_per_shard] sigma_MPI = local_parameters[(G_per_shard+1):(G_per_shard*2)];
	 	vector[S] exposure_rate = local_parameters[((M*2)+1):rows(local_parameters)];

		// Vectorise lpmf
		vector[symbol_end[G_per_shard+1]] lambda_MPI_c;
		vector[symbol_end[G_per_shard+1]] sigma_MPI_c;
		for(g in 1:G_per_shard){
			int how_many = symbol_end[g+1] - (symbol_end[g]);
			lambda_MPI_c[(symbol_end[g]+1):symbol_end[g+1]] = rep_vector(lambda_MPI[g], how_many);
			sigma_MPI_c [(symbol_end[g]+1):symbol_end[g+1]] = rep_vector(sigma_MPI[g],  how_many);
		}


		// Return

		return [
			neg_binomial_2_log_lpmf(
    		counts[1:symbol_end[G_per_shard+1]] |
    		exposure_rate[sample_idx[1:symbol_end[G_per_shard+1]]] +
    		lambda_MPI_c,
	    	sigma_MPI_c
    	)
    ]';

  }


}
data {

	// Reference matrix inference
  int<lower=0> N;
  int<lower=0> M;
	int<lower=0> G;
	int<lower=0> S;
  int n_shards;
	int<lower=0> counts[n_shards, N];
	int<lower=0> symbol_end[n_shards, M+1];
	int<lower=0> G_ind[n_shards, M];
	int<lower=0> sample_idx[n_shards, N];
	int<lower=0> G_per_shard[n_shards];
	int<lower=0> G_per_shard_idx[n_shards + 1];

	int<lower=0> CP; // Counts package size
  int<lower=0> counts_package[n_shards, CP];

	int<lower=1> C; // Covariates
	matrix[S,C] X; // Design matrix

	real<lower=0> lambda_mu_mu;

}
transformed data {

  vector[0] global_parameters;
  real real_data[n_shards, 0];

}
parameters {


   // Overall properties of the data
  real lambda_mu; // So is compatible with logGamma prior
  real<lower=0> lambda_sigma;
  vector[S] exposure_rate;

  // Gene-wise properties of the data
  vector[G] lambda_log_param;
  vector[G] sigma_raw_param;

  // Signa linear model

  real<upper=0> sigma_slope;
  real sigma_intercept;
  real<lower=0>sigma_sigma;


}
transformed parameters {
  // Sigma
  vector[G] sigma = 1.0 ./ exp(sigma_raw_param) ;
	vector[G] lambda_log = lambda_log_param;

	//matrix[]  =
}

model {

  lambda_mu ~ normal(lambda_mu_mu,2);
  lambda_sigma ~ normal(0,2);

  //sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

  sigma_intercept ~ normal(0,2);
  sigma_slope ~ normal(0,2);
  sigma_sigma ~ normal(0,2);

    // Gene-wise properties of the data
  lambda_log_param ~ exp_gamma_meanSd(lambda_mu,lambda_sigma);
  sigma_raw_param ~ normal(sigma_slope * lambda_log_param + sigma_intercept,sigma_sigma);

  // Exposure prior
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

	// Gene-wise properties of the data
	target += sum(map_rect(
		lp_reduce ,
		global_parameters ,
		get_reference_parameters_MPI(
			n_shards,
			M,
			G_per_shard,
			G_ind,
			lambda_log,
			sigma,
			exposure_rate
		),
		real_data,
		counts_package
	));


}
generated quantities{
	vector[G] counts_rng[S];

	for(g in 1:G) for(s in 1:S)
		counts_rng[s,g] =	neg_binomial_2_log_rng(exposure_rate[s] + lambda_log_param[g],	sigma[g]);

}
