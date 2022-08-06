data {
  int<lower=1> K;
  int<lower=1> J;
  int<lower=0> N;
  vector[J] x[N];
  vector[K] y[N];
  
}
parameters {
  vector[K] alpha;
  matrix[K, J] beta;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] L_sigma;

} transformed parameters{
  matrix[K, K] L_Sigma = diag_pre_multiply(L_sigma, L_Omega);
}
model {
  vector[K] mu[N];
  
  
  for(k in 1:K) {
    for(j in 1:J) {
       beta[k, j] ~ normal(0, 1);
    }
  }
  
  alpha ~ normal(0, 1);
  
  L_Omega ~ lkj_corr_cholesky(4);
  L_sigma ~ normal(0, 1);
  
  for (n in 1:N)
    mu[n] = alpha + beta * x[n];
  y ~ multi_normal_cholesky(mu, L_Sigma);
}

generated quantities {
  vector[N] log_lik;
  for(n in 1:N) {
    log_lik[n] = multi_normal_cholesky_lpdf(y[n] | alpha + beta * x[n], L_Sigma);
  }
}
