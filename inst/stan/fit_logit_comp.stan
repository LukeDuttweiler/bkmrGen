functions {
  matrix cov_exp_quad_aug(array[] vector x,
                            real alpha,
                            vector r) {
    int N = size(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    vector sqrt_r = sqrt(r);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha
                      * exp(-1 * dot_self((x[i] - x[j]) * sqrt_r));
        K[j, i] = K[i, j];
      }
    }
    K[N, N] = sq_alpha;
    return K;
  }
}

data {
  int<lower=1> N;
  int<lower=1> p;
  int<lower=1> d;
  matrix[N, d] X;
  row_vector[p] Z[N];
  int y[N];
}

transformed data {

}

parameters {
  vector[d] beta;
  real<lower=0> lambda;
  vector[N] h_tilde;
  int<lower=0,upper=1> delta[p];
  vector[p] rho;
  vector[p] r;
}

transformed parameters {
  for(i in 1:p){
    r[i] = (1/rho[i])*delta[i]
  }
  matrix[N, N] cov = cov_exp_quad_aug(Z, lambda, r);
  matrix[N, N] L_cov = cholesky_decompose(cov);
  vector[N] h_hat = L_cov * h_tilde;
  vector[N] eta = h_hat + X * beta;
}

model {
  for(i in 1:p){
    rho[i] ~ uniform(0,100)
  }
  for(i in 1:p){
    delta[i] ~ bernoulli(.25)
  }
  lambda ~ gamma(1, 0.1);
  beta ~ normal(0, 100);
  h_tilde ~ normal(0, 1);
  y ~ bernoulli_logit(eta);
}
