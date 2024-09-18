
int num_matches(array[] int x) {
  return sum(x);
}

array[] int get_matches(array[] int x) {
  array[size(x)] int ret;
  int pos = 0;
  for(i in 1:size(x)) if(x[i]) {
    pos += 1;
    ret[pos] = i;
  }
  return ret[:pos];
}

array[] int sep_matches(array[] int x) {
  array[size(x)] int ret;
  int pos1 = 1, pos2 = size(x);
  for(i in 1:size(x)) {
    if(x[i]) {
      ret[pos1] = i;
      pos1 += 1;
    } else {
      ret[pos2] = i;
      pos2 -= 1;
    }
  }
  return ret;
}

// == t(combn(indices, k)) in R
array[,] int combinations(array[] int indices, int k) {
  int n = size(indices);
  if(k > n) reject("k (",k,") must not be greater than n (",n,")!");
  else if(k == 0) {
    array[1,0] int m;
    return m;
  } else if(k == 1) {
    array[n,1] int m;
    m[,1] = indices;
    return m;
  } else if(k == n) {
    array[1,n] int m;
    m[1,] = indices;
    return m;
  } else {
    array[choose(n,k), k] int m;
    int o = 1;
    for(i in 1:(n-k+1)) {
      int cs = choose(n-i,k-1);
      for(j in (o):(o+cs-1)) {
        m[j,1] = indices[i];
      }
      m[(o):(o+cs-1),2:] = combinations(indices[(i+1):], k-1);
      o += cs;
    }
    return m;
  }
}

// whole-report probability for R after total exposure time t
real tvawpdf(data array[] int R, real t, vector t0_args, int K, vector v) {
  if(size(R) != size(v)) reject("R and v must have the same length (or number of columns)!");
  int nS = size(R);
  int nR = num_matches(R);
  array[nS] int indices = sep_matches(R);
  array[nR] int Rs = indices[:nR];
  array[nS-nR] int Us = indices[(nR+1):];
  //(Rs, Us) = sep_matches(R);
  if(K == 0 && nR == 0) {
    return 0.0;
  } else if(t <= 0 && nR > 0) {
    return negative_infinity();
  } else if(t <= 0 && nR == 0) {
    return 0.0;
  } else if(K < nR) {
    return negative_infinity();
  } else if(0 == nR && K > 0) {
    return tva_t_lccdf(t | sum(v), t0_args);
  } else if(0 < nR && nR < K && nR < nS && t > 0) {
    real xp = negative_infinity();
    real xm = negative_infinity();
    for(k in 0:nR) {
      int r = k == 0 ? 1 : choose(nR, k);
      array[r,k] int PR = combinations(Rs, k);
      for(l in 1:r) {
        real vsum = sum(v[PR[l,]]) + sum(v[Us]);
        if(k % 2 == 0) {
          xp = log_sum_exp(xp, tva_t_lpdf(t | vsum, t0_args) - log_psi(t, vsum, t0_args));
        } else {
          xm = log_sum_exp(xm, tva_t_lpdf(t | vsum, t0_args) - log_psi(t, vsum, t0_args));
        }
      }
    }
    if(xm > xp) {
      return negative_infinity();
    }
    return log_diff_exp(xp, xm);
  } else if(0 < nR && (nR == K || nR == nS) && t > 0) {
    array[nR] real ll2;
    for(j in 1:nR) {
      array[nR-1] int Rmi = append_array(Rs[:(j-1)], Rs[(j+1):]);
      real xp = negative_infinity();
      real xm = negative_infinity();
      for(k in 0:(nR-1)) {
        int r = k == 0 ? 1 : choose(nR-1, k);
        array[r,k] int PRmi = combinations(Rmi, k);
        for(l in 1:r) {
          real vsum = v[Rs[j]] + sum(v[PRmi[l,]]) + sum(v[Us]);
          if(k % 2 == 0) {
            xp = log_sum_exp(xp, tva_t_lcdf(t | vsum, t0_args) - log_psi(t, vsum, t0_args));
          } else {
            xm = log_sum_exp(xm, tva_t_lcdf(t | vsum, t0_args) - log_psi(t, vsum, t0_args));
          }
        }
      }
      if(xm > xp) {
        return negative_infinity();
      }
      ll2[j] = log_psi(t, v[Rs[j]], t0_args) + log_diff_exp(xp, xm);
    }
    return log_sum_exp(ll2);
  } else {
    reject("Unspecified scenario (R=",R,",t=",t,",t0=",t0_args,",K=",K,",v=",v,")");
  }
}

real tva_wr_log(data array[] int R, data array[] int S, real t, vector t0_args, vector K_args, vector v) {
  int nR = num_matches(R);
  int nS = num_matches(S);
  int max_nS = size(S);
  array[nS] int Ss = get_matches(S);
  vector[max_nS+1] ll;
  for(K in 0:max_nS) {
    ll[K+1] = tva_K_lpmf(K | K_args);
    if(ll[K+1] > negative_infinity()) {
      ll[K+1] += tvawpdf(R[Ss], t, t0_args, K, v);
    }
  }
  return log_sum_exp(ll);
}

real tva_pr_score_log(int score, data array[] int S, data array[] int D, real t, vector t0_args, vector K_args, vector v) {
  int nS = num_matches(S);
  int nD = num_matches(D);
  int nT = nS - nD;
  array[nD] int Ds = get_matches(D);
  array[size(S)] int T = S;
  if(nD > 0) T[Ds] = rep_array(0, nD);
  array[nT] int Ts = get_matches(T);
  if(score < 0 || score > nT) return negative_infinity();
  int r = score == 0 ? 1 : choose(nT, score);
  array[r,score] int PR = combinations(Ts, score);
  array[r] real ll; 
  for(i in 1:r) {
    array[size(S)] int R = rep_array(0, size(S));
    R[PR[i,]] = rep_array(1, score);
    ll[i] = tva_pr_log(R, S, D, t, t0_args, K_args, v);
  }
  return log_sum_exp(ll);
}

real tva_wr_score_log(int score, data array[] int S, real t, vector t0_args, vector K_args, vector v) {
  int nS = num_matches(S);
  int nT = nS;
  array[nS] int Ts = get_matches(S);
  if(score < 0 || score > nT) return negative_infinity();
  int r = score == 0 ? 1 : choose(nT, score);
  array[r,score] int PR = combinations(Ts, score);
  array[r] real ll; 
  for(i in 1:r) {
    array[size(S)] int R = rep_array(0, size(S));
    R[PR[i,]] = rep_array(1, score);
    ll[i] = tva_wr_log(R, S, t, t0_args, K_args, v);
  }
  return log_sum_exp(ll);
}

real tva_pr_score_predict(data array[] int S, data array[] int D, real t, vector t0_args, vector K_args, vector v) {
  real ret = 0.0;
  for(score in 0:size(S)) {
    ret += score * exp(tva_pr_score_log(score, S, D, t, t0_args, K_args, v));
  }
  return ret;
}

real tva_wr_score_predict(data array[] int S, real t, vector t0_args, vector K_args, vector v) {
  return tva_pr_score_predict(S, rep_array(0, size(S)), t, t0_args, K_args, v);
}


// partial-report probability for R after total exposure time t
real tvappdf(data array[] int R, data array[] int D, real t, vector t0_args, int K, vector v) {
  int nD = num_matches(D);
  int nR = num_matches(R);
  if(nD == 0) {
    return tvawpdf(R, t, t0_args, K, v);
  } else if(nR <= K) {
    array[size(R)] int M;
    array[nD] int Ds;
    Ds = get_matches(D);
    int q = min(K-nR, nD);
    array[q+1] real ls1;
    for(k in 0:q) {
      int r = k == 0 ? 1 : choose(nD, k);
      array[r,k] int RDs = combinations(Ds, k);
      array[r] real ls2;
      for(j in 1:r) {
        M[] = R[];
        for(m in 1:k) {
          M[RDs[j,m]] = 1;
        }
        ls2[j] = tvawpdf(M, t, t0_args, K, v);
      }
      ls1[k+1] = log_sum_exp(ls2);
    }
    return log_sum_exp(ls1);
  } else {
    return negative_infinity();
  }
}

real tva_pr_log(data array[] int R, data array[] int S, data array[] int D, real t, vector t0_args, vector K_args, vector v) {
  int nR = num_matches(R);
  int nS = num_matches(S);
  int max_nS = size(S);
  array[nS] int Ss = get_matches(S);
  vector[max_nS+1] ll;
  for(K in 0:max_nS) {
    ll[K+1] = tva_K_lpmf(K | K_args);
    if(ll[K+1] > negative_infinity()) {
      ll[K+1] += tvappdf(R[Ss], D[Ss], t, t0_args, K, v);
    }
  }
  return log_sum_exp(ll);
}

array[] int tva_wr_rng(array[] int S, real t, vector t0_args, vector K_args, vector v) {
  int nS = num_matches(S);
  if(size(v) != nS) reject("v must have as many rates as there are items!");
  array[nS] int Ss = get_matches(S);
  array[size(S)] int R = rep_array(0, size(S));
  int K = tva_K_rng(K_args);
  vector[nS] processing_times = tva_t_rng(v, t0_args);
  array[nS] int processing_order = sort_indices_asc(processing_times);
  for(j in 1:min(nS,K)) {
    if(processing_times[processing_order[j]] > t) break;
    R[Ss[processing_order[j]]] = 1;
  }
  return R;
}

int tva_wr_score_rng(array[] int S, real t, vector t0_args, vector K_args, vector v) {
  return sum(tva_wr_rng(S, t, t0_args, K_args, v));
}

array[] int tva_pr_rng(array[] int S, array[] int D, real t, vector t0_args, vector K_args, vector v) {
  int nS = size(S);
  if(size(D) != size(S)) reject("D and S must have the same length!");
  array[nS] int R = tva_wr_rng(S, t, t0_args, K_args, v);
  for(i in 1:nS) if(D[i]) R[i] = 0;
  return R;
}

int tva_pr_score_rng(array[] int S, array[] int D, real t, vector t0_args, vector K_args, vector v) {
  return sum(tva_pr_rng(S, D, t, t0_args, K_args, v));
}

