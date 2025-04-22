
real tva_K_lpmf(int x, vector K_args, data int max_K) {
  if(x < 0 || x > max_K) return negative_infinity();
  else return binomial_lpmf(x | max_K, K_args[1]);
}

int tva_K_rng(vector K_args, data int max_K) {
  return binomial_rng(max_K, K_args[1]);
}
