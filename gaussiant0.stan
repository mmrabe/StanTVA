

real psi(real t, real l, vector args) {
  return l;
}

real log_psi(real t, real l, vector args) {
  return log(l);
}

real Psi(real t, real l, vector args) {
  return t*l;
}

real tva_t_lpdf(real x, real lambda, vector t0_args) {
  real mu = t0_args[1];
  real sigma = t0_args[2];
  return exp_mod_normal_lpdf(x | mu, sigma, lambda);
}

real tva_t_lcdf(real x, real lambda, vector t0_args) {
  real mu = t0_args[1];
  real sigma = t0_args[2];
  return exp_mod_normal_lcdf(x | mu, sigma, lambda);
}

real tva_t_lccdf(real x, real lambda, vector t0_args) {
  real mu = t0_args[1];
  real sigma = t0_args[2];
  return exp_mod_normal_lccdf(x | mu, sigma, lambda);
}

vector tva_t_rng(vector lambda, vector t0_args) {
  real mu = t0_args[1];
  real sigma = t0_args[2];
  return to_vector(exponential_rng(lambda)) + normal_rng(mu, sigma);
}


