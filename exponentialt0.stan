
real log_psi(real x, real v, vector args) {
  return log(v);
}

real tva_t_lpdf(real x, real v, vector t0_args) {
  real lambda1 = t0_args[1];
  if(is_nan(lambda1) || lambda1 < 0) reject("t0_lambda is ", lambda1, " but should be non-zero positive finite!");
  else if(is_nan(v) || v < 0) reject("v is ", v, " but should be non-zero positive finite!");
  else if(lambda1 <= machine_precision() || v <= machine_precision()) return negative_infinity();
  else if(is_inf(lambda1) && is_inf(v)) reject("t0_lambda and/or v must be finite but both are infinite!");
  else if(is_inf(lambda1)) return exponential_lpdf(x | v);
  else if(is_inf(v)) return exponential_lpdf(x | lambda1);
  else if(lambda1 >= v + machine_precision()) {
    real log_lambda1 = log(lambda1);
    real log_v = log(v);
    real d2 = log_diff_exp(-x*v,-x*lambda1);
    if(d2 != negative_infinity()) return log_lambda1+log_v-log(lambda1-v)+d2;
    else return gamma_lpdf(x | 2.0, lambda1);
  } else if(lambda1 <= v - machine_precision()) {
    real log_lambda1 = log(lambda1);
    real log_v = log(v);
    real d2 = log_diff_exp(-x*lambda1,-x*v);
    if(d2 != negative_infinity()) return log_lambda1+log_v-log(v-lambda1)+d2;
    else return gamma_lpdf(x | 2.0, lambda1);
  }
  else return gamma_lpdf(x | 2.0, lambda1);
}

real tva_t_lcdf(real x, real v, vector t0_args) {
  real lambda1 = t0_args[1];
  if(is_nan(lambda1) || lambda1 < 0) reject("t0_lambda is ", lambda1, " but should be non-zero positive finite!");
  else if(is_nan(v) || v < 0) reject("v is ", v, " but should be non-zero positive finite!");
  else if(is_inf(lambda1) && is_inf(v)) reject("t0_lambda and/or v must be finite but both are infinite!");
  else if(lambda1 <= machine_precision() || v <= machine_precision()) return negative_infinity();
  else if(is_inf(lambda1)) return exponential_lcdf(x | v);
  else if(is_inf(v)) return exponential_lcdf(x | lambda1);
  else if(v < t0_args[1] || v > t0_args[1]) return log1m_exp(tva_t_lccdf(x | v, t0_args));
  else return gamma_lcdf(x | 2.0, lambda1);
}

real tva_t_lccdf(real x, real v, vector t0_args) {
  real lambda1 = t0_args[1];
  if(is_nan(lambda1) || lambda1 < 0) reject("t0_lambda is ", lambda1, " but should be non-zero positive finite!");
  else if(is_nan(v) || v < 0) reject("v is ", v, " but should be non-zero positive finite!");
  else if(is_inf(lambda1) && is_inf(v)) reject("t0_lambda and/or v must be finite but both are infinite!");
  else if(lambda1 <= machine_precision() || v <= machine_precision()) return negative_infinity();
  else if(is_inf(lambda1)) return exponential_lccdf(x | v);
  else if(is_inf(v)) return exponential_lccdf(x | lambda1);
  else if(lambda1 >= v + machine_precision()) return log_diff_exp(log(lambda1)-v*x, log(v)-lambda1*x)-log(lambda1-v);
  else if(lambda1 <= v - machine_precision()) return log_diff_exp(log(v)-lambda1*x, log(lambda1)-v*x)-log(v-lambda1);
  else return gamma_lccdf(x | 2.0, lambda1);
}


vector tva_t_rng(vector v, vector t0_args) {
  real lambda1 = t0_args[1];
  return to_vector(exponential_rng(v)) + exponential_rng(lambda1);
}


