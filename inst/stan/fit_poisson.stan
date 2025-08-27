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
  //matrix[N, N] cov = cov_exp_quad(Z, lambda, rho);
  //matrix[N, N] L_cov = cholesky_decompose(cov);
  //vector[N] h_hat = L_cov * h_tilde;
  //vector[N] eta = h_hat + X * beta;
  //vector<lower=0>[N] mu = exp(eta);
  vector[N] h_hat = cholesky_decompose(cov_exp_quad(Z, lambda, rho) + diag_matrix(rep_vector(1e-6, N)))*h_tilde;

  vector[N] ystar = h_hat + X * beta;
}

model {
  lambda ~ gamma(1, .1);
  beta ~ normal(0, 100);
  h_tilde ~ normal(0, 1);
  y ~ poisson(exp(ystar));
}
