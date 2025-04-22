
real tva_K_lpmf(int x, vector K_args, data int max_K) {
  real aK = K_args[1];
  real bK = K_args[2];
  if(x < 0 || x > max_K) return negative_infinity();
  else return beta_binomial_lpmf(x | max_K, aK, bK);
}

int tva_K_rng(vector K_args, data int max_K) {
  real aK = K_args[1];
  real bK = K_args[2];
  return beta_binomial_rng(max_K, aK, bK);
}
