data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=0> N;
  int<lower=1> M; 
  array[N] vector[J] x;
  array[N] vector[K] y;
  
}
parameters {
  vector[K] alpha;
  matrix[K, J] beta;
  
  matrix[K, M] Gamma; //
  vector<lower=0>[K] lambda;
  real<lower=0> lambda_mean;
}
transformed parameters {
  matrix[K, K] lambda_inv = diag_matrix(1.0 ./ lambda);
  cov_matrix[K] Sigma_inv = lambda_inv - lambda_inv * Gamma * inverse_spd(diag_matrix(rep_vector(1, M)) +  Gamma' * lambda_inv * Gamma) * Gamma' * lambda_inv;
}
model {
  array[N] vector[K] mu;
  
  for(k in 1:K) {
    for(m in 1:M) {
      Gamma[k, m] ~ normal(0, 1);
    }
    
    for(j in 1:J) {
       beta[k, j] ~ normal(0, 1);
    }
  }
  
  alpha ~ normal(0, 10);

  lambda_mean ~ cauchy(0, 1);
  lambda ~ lognormal(log(lambda_mean), 1);
  
  for (n in 1:N)
    mu[n] = alpha + beta * x[n];
  y ~ multi_normal_prec(mu, Sigma_inv);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) {
    log_lik[n] = multi_normal_prec_lpdf(y[n] | alpha + beta * x[n], Sigma_inv);
  }
}
