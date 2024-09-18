
real tva_K_lpmf(int x, vector K_args) {
  real K = K_args[1];
  if(x <= K && x >= K - 1.0) return log1p(x-K);
  else if(x >= K && x <= K + 1.0) return log1p(K-x);
  else return negative_infinity();
}

int tva_K_rng(vector K_args) {
  int k = 0;
  real K = K_args[1];
  while(k < K) {
    if(k >= K - 1.0) {
      return k + bernoulli_rng(K-k);
    }
    k += 1;
  }
  reject("Illegal K!");
}
