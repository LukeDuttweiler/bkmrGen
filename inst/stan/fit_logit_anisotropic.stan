data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=1> d;
  matrix[N, d] X;
  array[N] vector[p] Z;
  int y[N];
  int<lower=0> M;
}

transformed data {

}

parameters {
  vector[d] beta;
  real<lower=0> lambda;
  vector[N] h_tilde;
  vector[p] zeta;
  simplex[p] theta;
}

transformed parameters {
  // matrix[N, N] cov = cov_exp_quad(Z, lambda, rho);
  // matrix[N, N] L_cov = cholesky_decompose(cov);
  // vector[N] h_hat = L_cov * h_tilde;
  // vector[N] eta = h_hat + X * beta;
  // ystar = eta, calculated in two lines to avoid having to save all the extra junk
  array[p] real<lower=0> l; // transform rho (adatpive parameter) to match stan requirements

  for(i in 1:p){
    l[i] = inv(sqrt(2)*exp(theta[i] * zeta[i])); // rho = exp(zeta)^theta
  }

  vector[N] h_hat = cholesky_decompose(gp_exp_quad_cov(Z, lambda, l) + diag_matrix(rep_vector(1e-6, N)))*h_tilde;

  vector[N] ystar = h_hat + X * beta;
}

model {
  theta ~ dirichlet(rep_vector(1.0, p));
  target += gamma_lpdf(exp(zeta) | 1, 1) + sum(zeta);   // Gamma on rho, Jacobian adjustment
  lambda ~ gamma(1, 0.1);
  beta ~ normal(0, 100);
  h_tilde ~ normal(0, 1);
  target += M * bernoulli_logit_lpmf(y | ystar);
}

generated quantities {

}
