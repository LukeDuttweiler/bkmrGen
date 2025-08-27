functions {
  matrix cov_exp_quad_aug(matrix x,
                            real alpha,
                            row_vector r) {
    int N = rows(x);
    int D = cols(x);
    matrix[N, N] K;
    real sq_alpha = square(alpha);
    row_vector[D] sqrt_r = sqrt(r);
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha;
      for (j in (i + 1):N) {
        real sq_dist = dot_self((x[i] - x[j]) .* sqrt_r);
        K[i, j] = sq_alpha * exp(-sq_dist);
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
  matrix[N, p] Z;
  int y[N];
}

transformed data {

}

parameters {
  vector[d] beta;
  real<lower=0> lambda;
  real<lower=0,upper=1> delta[p];
  row_vector<lower=.01>[p] rho;
  vector[N] h_tilde;
}

transformed parameters {
  row_vector[p] r;
  for(i in 1:p){
    r[i] = 1/rho[i];
  }

  //matrix[N, N] cov = cov_exp_quad_aug(Z, lambda, r) + diag_matrix(rep_vector(1e-6, N));
  //matrix[N, N] L_cov = cholesky_decompose(cov);
  //vector[N] h_hat = L_cov * h_tilde;
  //vector[N] eta = h_hat + X * beta;
  //ystar = eta, calculated in one line to avoid saving all the extra junk
  vector[N] h_hat = cholesky_decompose(cov_exp_quad_aug(Z, lambda, r) + diag_matrix(rep_vector(1e-6, N))) * h_tilde;

  vector[N] ystar = h_hat + X * beta;
}

model {
  for(i in 1:p){
    delta[i] ~ beta(2, 6);
  }

  for(i in 1:p){
    target += log_mix(delta[i],
                      uniform_lpdf(rho[i] | .01, 100),
                      lognormal_lpdf(rho[i] | 10, .5));
  }

  lambda ~ gamma(1, 0.1);
  beta ~ normal(0, 100);
  h_tilde ~ normal(0, 1);
  y ~ bernoulli_logit(ystar);
}
