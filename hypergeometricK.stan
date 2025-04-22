
real tva_K_lpmf(int x, vector K_args, data int max_K) {
  real gK = K_args[1];
  real bK = K_args[2];
  if(x < 0 || x > max_K) return negative_infinity();
  else if(x == max_K) return hypergeometric_lccdf(max_K-1 | gK, bK);
  else return hypergeometric_lpmf(x | gK, bK);
}

int tva_K_rng(vector K_args, data int max_K) {
  vector[max_K+1] ps;
  for(i in 0:max_K) ps[i+1] = exp(tva_K_lpmf(i, K_args, max_K));
  return categorical_rng(ps)-1;
}
