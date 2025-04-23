
real tva_K_lpmf(int x, vector K_args, data int max_K) {
  real mK = K_args[1];
  real sK = K_args[2];
  if(x < 0 || x > max_K) return negative_infinity();
  else {
    vector[max_K] cutpoints = linspaced_vector(max_K, sK, max_K*sK);
    return ordered_probit_lpmf(x + 1 | mK * sK, cutpoints);
  }
}

int tva_K_rng(vector K_args, data int max_K) {
  real mK = K_args[1];
  real sK = K_args[2];
  vector[max_K] cutpoints = linspaced_vector(max_K, sK, max_K*sK);
  return ordered_probit_rng(mK * sK, cutpoints) - 1;
}
