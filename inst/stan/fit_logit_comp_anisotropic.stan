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
  real<lower=0,upper=1> delta[p];
  vector[p] zeta;
  simplex[p] theta;
  vector[N] h_tilde;
}

transformed parameters {
  //matrix[N, N] cov = cov_exp_quad_aug(Z, lambda, r) + diag_matrix(rep_vector(1e-6, N));
  //matrix[N, N] L_cov = cholesky_decompose(cov);
  //vector[N] h_hat = L_cov * h_tilde;
  //vector[N] eta = h_hat + X * beta;
  //ystar = eta, calculated in one line to avoid saving all the extra junk
  array[p] real<lower=0> r;
  array[p] real<lower=0> l; // transform r (adatpive parameter) to match stan requirements

  for(i in 1:p){
    r[i] = sqrt(2)*exp(theta[i] * zeta[i]);
    l[i] = inv(r[i]); // rho = exp(zeta)^theta
  }

  vector[N] h_hat = cholesky_decompose(gp_exp_quad_cov(Z, lambda, l) + diag_matrix(rep_vector(1e-6, N))) * h_tilde;

  vector[N] ystar = h_hat + X * beta;
}

model {
  theta ~ dirichlet(rep_vector(5.0, p));
  for(i in 1:p){
    delta[i] ~ beta(2, 6);
  }

  for(i in 1:p){
    target += log_mix(delta[i],
                      gamma_lpdf(exp(zeta[i]) | 1, 1) + sum(zeta), // Jacobian
                      normal_lpdf(zeta[i] | -30, .5));
  }

  lambda ~ gamma(1, 0.1);
  beta ~ normal(0, 100);
  h_tilde ~ normal(0, 1);
  y ~ bernoulli_logit(ystar);
}
