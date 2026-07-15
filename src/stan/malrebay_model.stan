// =============================================================================
// malrebay_model.stan
// Bayesian recrudescence model for genotyping markers
//
// Goal: classify each recurrence as recrudescence (same infection) or
// reinfection (new infection) using Day 0 vs recurrence genotypes.
//
// For each patient, we compute a summary log-likelihood ratio (SLR):
//   SLR > 0 → supports recrudescence
//   SLR < 0 → supports reinfection
//
// Posterior probability = inv_logit(SLR)
// Prior P(recrudescence) = 0.5 (data-driven classification)
// =============================================================================

functions {
  // Allele similarity:
  // 1 = microsatellite (distance-based decay)
  // 2 = MSP/GLURP (same vs different family)
  // 3 = exact match only
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
  // N patients, J loci, up to maxMOI clones per locus
  // Alleles encoded as integers (1..K[j])
  // hidden0/hiddenf: unobserved allele slots
  // comparable[i,j]: locus observed at both timepoints
  // additional_counts: background data for allele frequencies
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
  // qq: within-family mismatch (genotyping error)
  // qq_crossfamily: cross-family similarity (rare)
  // d_param: decay rate for microsatellite distance
  // freq[j]: allele frequencies per locus
  real<lower=0,upper=1> qq;
  real<lower=0,upper=1> qq_crossfamily;
  real<lower=0,upper=1> d_param;
  array[J] simplex[max_K] freq;
}
transformed parameters {
  // 1. Build similarity distribution from d_param
  // 2. Precompute locus-level log-probability matrices
  // 3. Compute SLR per patient:
  //    - compare all allele pairs
  //    - handle hidden alleles via marginalisation
  //    - divide by freq → likelihood ratio
  //    - sum across loci (assumes independence)
  vector[max_dist] log_dvect;
  {
    vector[max_dist] dvect;
    for (d in 1:max_dist) dvect[d] = (1.0 - d_param)^(d - 1) * d_param;
    dvect     = dvect / (sum(dvect) + 1e-10);
    log_dvect = log(dvect + 1e-10);
  }

  real l_q_cf   = log(qq_crossfamily + 1e-10);
  real l1m_q_cf = log1m(qq_crossfamily);

  // Per-locus precomputed matrices
  array[J] vector[max_K]        log_freq;
  array[J] matrix[max_K, max_K] log_dist_mat;
  array[J] vector[max_K]        log_dist_row_sums;

  // Per-locus constants computed once per leapfrog step rather than per patient. 
  // log_both_hidden[j] : LR contribution when both Day-0 and recurrence alleles are hidden. 
  // log_col_sum[j][kf] : LR contribution when only the Day-0 allele is hidden, indexed by kf.
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

    // Precompute per-locus constants
    log_both_hidden[j] = log_sum_exp(log_freq[j, 1:K[j]] + log_dist_row_sums[j, 1:K[j]]);
    for (kf in 1:K[j])
      log_col_sum[j][kf] = log_sum_exp(log_freq[j, 1:K[j]] + col(log_dist_mat[j], kf)[1:K[j]]);
  }

  // Per-patient summary log-likelihood ratios (SLR) computed here so the model 
  // block stays concise and generated quantities avoids redundant computation.
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
              lpr[p_idx] = log_col_sum[j][kf] - log_freq[j, kf];
            } else if (h0 == 0 && hf == 1) {
              lpr[p_idx] = log_dist_row_sums[j, k0];
            } else {
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
  // Priors:
  // qq ~ Beta(1,1)
  // qq_crossfamily ~ Beta(1,1000)
  // d_param ~ Beta(2,2)
  // freq ~ Dirichlet (with additional data)

  // Likelihood:
  // 50/50 mixture of recrudescence vs reinfection
  // log_sum_exp(slr, 0) implements this
  qq             ~ beta(1, 1);
  qq_crossfamily ~ beta(1, 1000);
  d_param        ~ beta(2, 2);
  for (j in 1:J) {
    vector[max_K] alpha = rep_vector(0.1, max_K);
    for (k in 1:K[j]) alpha[k] = 1.0;
    freq[j]              ~ dirichlet(alpha);
    additional_counts[j] ~ multinomial(freq[j]);
  }
  for (i in 1:N)
    target += log_sum_exp(slr_vec[i], 0.0);
}


generated quantities {
  // Convert each patient's SLR back to a probability on [0, 1].
  // inv_logit(slr) is equivalent to LR / (1 + LR), which is the posterior
  // probability of recrudescence under the equal 0.5 / 0.5 prior.
  vector[N] p_recrud = inv_logit(slr_vec);
}
