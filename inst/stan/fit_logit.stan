data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=1> d;
  matrix[N, d] X;
  row_vector[p] Z[N];
  int y[N];
  real<lower=0> rho;
}

transformed data {

}

parameters {
  vector[d] beta;
  real<lower=0> lambda;
  vector[N] h_tilde;
}

transformed parameters {
  matrix[N, N] cov = cov_exp_quad(Z, lambda, rho);
  matrix[N, N] L_cov = cholesky_decompose(cov);
  vector[N] h_hat = L_cov * h_tilde;
  vector[N] eta = h_hat + X * beta;
}

model {
  lambda ~ gamma(3, 3);
  beta ~ normal(0, 10);
  h_tilde ~ normal(0, 1);
  y ~ bernoulli_logit(eta);
}
