
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
  if(x <= 0) return negative_infinity();
  return exponential_lpdf(x | lambda);
}

real tva_t_lcdf(real x, real lambda, vector t0_args) {
  if(x <= 0) return negative_infinity();
  return exponential_lcdf(x | lambda);
}

real tva_t_lccdf(real x, real lambda, vector t0_args) {
  if(x <= 0) return 0.0;
  return exponential_lccdf(x | lambda);
}

vector tva_t_rng(vector lambda, vector t0_args) {
  return to_vector(exponential_rng(lambda));
}
