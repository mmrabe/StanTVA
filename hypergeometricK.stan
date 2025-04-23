
real generalized_hypergeometric_lpmf(int n, int N, real a, real b) {
  return lchoose(a,n)+lchoose(b,N-n)-lchoose(a+b,N);
}

real hypergeometric_lcdf(int n, int N, real a, real b) {
  vector[n+1] s;
  for(i in 0:n) s[i+1] = generalized_hypergeometric_lpmf(i | N, a, b);
  return log_sum_exp(s);
}

real hypergeometric_lccdf(int n, int N, real a, real b) {
  return log1m_exp(hypergeometric_lcdf(n | N, a, b));
}

real tva_K_lpmf(int x, vector K_args, data int max_K) {
  real gK = K_args[1];
  real bK = K_args[2];
  if(x < 0 || x > max_K) return negative_infinity();
  else if(x == 0 && max_K == 0) return 0.0;
  else if(x == max_K) return hypergeometric_lccdf(max_K-1 | max_K, gK, bK);
  else return generalized_hypergeometric_lpmf(x | max_K, gK, bK);
}

int tva_K_rng(vector K_args, data int max_K) {
  vector[max_K+1] ps;
  for(i in 0:max_K) ps[i+1] = exp(tva_K_lpmf(i | K_args, max_K));
  return categorical_rng(ps)-1;
}
