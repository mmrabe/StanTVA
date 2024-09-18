
real tva_K_lpmf(int x, vector K_args) {
  if(x < 0 || x > size(K_args)) return negative_infinity();
  real p = K_args[x+1];
  if(p < machine_precision()) p = machine_precision();
  if(p > 1.0 - machine_precision()) p = 1-machine_precision();
  return log(p);
}

int tva_K_rng(vector K_args) {
  return categorical_rng(K_args)-1;
}
