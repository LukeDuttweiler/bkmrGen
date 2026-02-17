data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=1> d;
  matrix[N, d] X;
  array[N] vector[p] Z;
  int y[N];
}

transformed data {

}

parameters {
  vector[d] beta;
  real<lower=0> lambda;
  vector[N] h_tilde;
  real<lower=1e-4> rho;
}

transformed parameters {
  // matrix[N, N] cov = cov_exp_quad(Z, lambda, rho);
  // matrix[N, N] L_cov = cholesky_decompose(cov);
  // vector[N] h_hat = L_cov * h_tilde;
  // vector[N] eta = h_hat + X * beta;
  // ystar = eta, calculated in two lines to avoid having to save all the extra junk
  // vector[N] h_hat = cholesky_decompose(cov_exp_quad(Z, lambda, rho) + diag_matrix(rep_vector(1e-6, N)))*h_tilde;
  real<lower=0> l = inv(sqrt(2)*pow(rho, 1.0/p)); // transform rho (adatpive parameter) to match stan requirements

  vector[N] h_hat = cholesky_decompose(gp_exp_quad_cov(Z, lambda, l) + diag_matrix(rep_vector(1e-6, N)))*h_tilde;

  vector[N] ystar = h_hat + X * beta;
}

model {
  rho ~ gamma(1, 1);
  lambda ~ gamma(1, 0.1);
  beta ~ normal(0, 100);
  h_tilde ~ normal(0, 1);
  y ~ bernoulli_logit(ystar);
}

generated quantities {

}
