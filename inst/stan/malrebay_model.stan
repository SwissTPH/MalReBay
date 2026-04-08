// =============================================================================
// inst/stan/malrebay_model.stan
// MalReBay -- length-polymorphic Bayesian recrudescence model
// =============================================================================

functions {
  real dist_log_prob_logic(data real dist, vector log_dvect, data int method,
                           real l_q_cf, real l1m_q_cf, data real threshold) {
    if (method == 1) {
      int d = to_int(round(dist));
      if (d < 0 || d >= num_elements(log_dvect)) return -23.0;
      return log_dvect[d + 1];
    } else if (method == 2) {
      return dist <= threshold ? l1m_q_cf : l_q_cf;
    } else {
      return dist < 0.5 ? l1m_q_cf : l_q_cf;
    }
  }
}
data {
  int<lower=1> N; int<lower=1> J; int<lower=1> maxMOI; int<lower=1> max_K; int<lower=1> max_dist;
  array[J] int<lower=1> K;
  array[N, J*maxMOI] int<lower=0> recoded0; array[N, J*maxMOI] int<lower=0> recodedf;
  array[N, J*maxMOI] int<lower=0,upper=1> hidden0; array[N, J*maxMOI] int<lower=0,upper=1> hiddenf;
  array[N] int<lower=1> MOI0; array[N] int<lower=1> MOIf;
  array[J, max_K, max_K] real<lower=0> dist_array;
  array[J] int<lower=1,upper=3> method_int;
  vector[J] threshold;
  array[N, J] int<lower=0,upper=1> comparable;
  array[J, max_K] int<lower=0> additional_counts;
}
parameters {
  real<lower=0,upper=1> qq;
  real<lower=0,upper=1> qq_crossfamily;
  real<lower=0,upper=1> d_param;
  // BUG 1 FIX: theta_recrud removed -- implicit 50/50 prior per patient,
  // matching the old R MCMC. A shared theta caused a positive feedback loop
  // pulling all p_recrud toward 1.
  array[J] simplex[max_K] freq;
}
transformed parameters {
  // --- Distance distribution ---
  vector[max_dist] log_dvect;
  {
    vector[max_dist] dvect;
    for (d in 1:max_dist) dvect[d] = (1.0 - d_param)^(d - 1) * d_param;
    dvect     = dvect / (sum(dvect) + 1e-10);
    log_dvect = log(dvect + 1e-10);
  }

  real l_q_cf   = log(qq_crossfamily + 1e-10);
  real l1m_q_cf = log1m(qq_crossfamily);
  // Note: l_q (log qq) was removed -- BUG 3 FIX showed qq cancels in all LR cases.

  // --- Per-locus precomputed matrices ---
  array[J] vector[max_K]        log_freq;
  array[J] matrix[max_K, max_K] log_dist_mat;
  array[J] vector[max_K]        log_dist_row_sums;

  // Fix #3: locus-level constants -- computed once per leapfrog step, not per patient.
  // log_both_hidden[j] : LR numerator when both alleles are hidden (h0=1, hf=1).
  // log_col_sum[j][kf] : LR numerator when only h0 is hidden (h0=1, hf=0), indexed by kf.
  array[J] real           log_both_hidden;
  array[J] vector[max_K]  log_col_sum;

  for (j in 1:J) {
    log_dist_mat[j]      = rep_matrix(-23.0, max_K, max_K);
    log_dist_row_sums[j] = rep_vector(-23.0, max_K);
    log_freq[j]          = log(freq[j] + 1e-10);
    log_col_sum[j]       = rep_vector(-23.0, max_K);

    for (k1 in 1:K[j]) {
      for (k2 in 1:K[j]) {
        log_dist_mat[j][k1, k2] = dist_log_prob_logic(
          dist_array[j, k1, k2], log_dvect, method_int[j],
          l_q_cf, l1m_q_cf, threshold[j]
        );
      }
      log_dist_row_sums[j][k1] = log_sum_exp(log_dist_mat[j][k1, 1:K[j]]);
    }

    // Precompute per-locus constants (fix #3)
    log_both_hidden[j] = log_sum_exp(log_freq[j, 1:K[j]] + log_dist_row_sums[j, 1:K[j]]);
    for (kf in 1:K[j])
      log_col_sum[j][kf] = log_sum_exp(log_freq[j, 1:K[j]] + col(log_dist_mat[j], kf)[1:K[j]]);
  }

  // Fix #2: move per-patient SLR into transformed parameters so the model block
  // is trivial and generated quantities requires no re-computation.
  vector[N] slr_vec;
  for (i in 1:N) {
    real slr = 0.0;
    for (j in 1:J) {
      if (comparable[i, j] == 1) {
        int off = (j - 1) * maxMOI;
        vector[MOI0[i] * MOIf[i]] lpr;
        int p_idx = 1;
        for (a0 in 1:MOI0[i]) {
          for (af in 1:MOIf[i]) {
            int k0 = recoded0[i, off + a0]; int kf = recodedf[i, off + af];
            int h0 = hidden0[i, off + a0];  int hf = hiddenf[i, off + af];
            if (h0 == 0 && hf == 0) {
              lpr[p_idx] = log_dist_mat[j, k0, kf] - log_freq[j, kf];
            } else if (h0 == 1 && hf == 0) {
              // BUG 3 FIX: l_q cancelled; use precomputed log_col_sum[j][kf]
              lpr[p_idx] = log_col_sum[j][kf] - log_freq[j, kf];
            } else if (h0 == 0 && hf == 1) {
              // BUG 3 FIX: l_q cancelled; denominator sums to 1
              lpr[p_idx] = log_dist_row_sums[j, k0];
            } else {
              // BUG 3 FIX: both hidden; use precomputed log_both_hidden[j]
              lpr[p_idx] = log_both_hidden[j];
            }
            p_idx += 1;
          }
        }
        slr += log_sum_exp(lpr) - log(MOI0[i] * MOIf[i]);
      }
    }
    slr_vec[i] = slr;
  }
}
model {
  qq             ~ beta(1, 1);
  qq_crossfamily ~ beta(1, 1000);
  d_param        ~ beta(2, 2);
  // BUG 1 FIX: theta_recrud prior removed
  for (j in 1:J) {
    vector[max_K] alpha = rep_vector(0.1, max_K);
    for (k in 1:K[j]) alpha[k] = 1.0;
    freq[j]              ~ dirichlet(alpha);
    additional_counts[j] ~ multinomial(freq[j]);
  }
  // BUG 1 FIX: log_sum_exp(slr, 0.0) = log(exp(slr) + 1), equivalent to theta=0.5
  for (i in 1:N)
    target += log_sum_exp(slr_vec[i], 0.0);
}
generated quantities {
  // Fix #2: slr_vec already computed in transformed parameters -- trivial here.
  // Fix #4 (consistency): inv_logit(slr) = LR/(1+LR), equivalent to theta=0.5
  vector[N] p_recrud = inv_logit(slr_vec);
}
