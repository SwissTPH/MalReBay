// =============================================================================
// malrebay_ampseq_model.stan - Clean "Best-Match" Logic
// Matches original R Gibbs math precisely
// =============================================================================

data {
  int<lower=1> N;
  int<lower=1> J;
  int<lower=1> maxMOI;
  int<lower=1> max_K;
  array[J] int<lower=1> K;
  array[N, J] int<lower=0> MOI0;
  array[N, J] int<lower=0> MOIf;
  array[N, J*maxMOI] int<lower=0> recoded0;
  array[N, J*maxMOI] int<lower=0> recodedf;
  array[N, J*maxMOI] int<lower=0, upper=1> hidden0;
  array[N, J*maxMOI] int<lower=0, upper=1> hiddenf;
  array[N, J] int<lower=0, upper=1> comparable;
  array[J, max_K] int<lower=0> additional_counts;
}
transformed data {
  // All integer quantities here depend only on data arrays (no parameters).
  // Moved from model block so they are computed ONCE at startup, not at every
  // leapfrog step.

  array[N, J] int td_n_lost      = rep_array(0, N, J);
  array[N, J] int td_m_i_0       = rep_array(0, N, J);
  array[N, J] int td_m_i         = rep_array(0, N, J);
  array[N, J] int td_n_alleles_f = rep_array(0, N, J);
  // td_has_match[i, off+af] = 1 if recurrence allele af at locus j has an
  // observed Day-0 match for patient i; 0 otherwise (or if af is hidden/zero).
  array[N, J*maxMOI] int td_has_match = rep_array(0, N, J * maxMOI);

  for (i in 1:N) {
    for (j in 1:J) {
      if (comparable[i, j] == 0 || MOI0[i, j] == 0 || MOIf[i, j] == 0)
        continue;

      int off = (j - 1) * maxMOI;

      // n_lost: observed Day-0 alleles absent from all observed recurrence alleles
      for (a0 in 1:MOI0[i, j]) {
        int k0 = recoded0[i, off + a0];
        if (k0 == 0 || hidden0[i, off + a0] == 1) continue;
        int found = 0;
        for (af in 1:MOIf[i, j]) {
          if (recodedf[i, off + af] == k0 && hiddenf[i, off + af] == 0) {
            found = 1; break;
          }
        }
        if (!found) td_n_lost[i, j] += 1;
      }

      // BUG 1 FIX: split hidden counts into m_i_0 (Day-0 only) and m_i (total)
      for (a0 in 1:MOI0[i, j]) {
        if (hidden0[i, off + a0] == 1) { td_m_i_0[i, j] += 1; td_m_i[i, j] += 1; }
      }
      for (af in 1:MOIf[i, j]) {
        if (hiddenf[i, off + af] == 1) td_m_i[i, j] += 1;
      }

      // n_alleles_f and has_match per recurrence allele slot
      for (af in 1:MOIf[i, j]) {
        int kf = recodedf[i, off + af];
        if (kf == 0 && hiddenf[i, off + af] == 0) continue;
        td_n_alleles_f[i, j] += 1;

        if (hiddenf[i, off + af] == 0) {
          int match = 0;
          for (a0 in 1:MOI0[i, j]) {
            if (recoded0[i, off + a0] == kf && hidden0[i, off + a0] == 0) {
              match = 1; break;
            }
          }
          td_has_match[i, off + af] = match;
        }
        // hidden recurrence allele: td_has_match stays 0 (not used in that branch)
      }
    }
  }
}
parameters {
  real<lower=0, upper=1> q_mismatch;
  real<lower=0, upper=1> q_loss;
  real<lower=0, upper=1> q_dropout;
  array[J] simplex[max_K] freq;
}
transformed parameters {
  array[J] vector[max_K] log_freq;
  for (j in 1:J) log_freq[j] = log(freq[j] + 1e-10);

  real l_match    = log1m(q_mismatch);
  real l_mismatch = log(q_mismatch + 1e-10);
  real l_dropout  = log(q_dropout  + 1e-10);
  real l_loss     = log(q_loss     + 1e-10);  // Fix #5: precompute alongside the others

  // Fix #2: compute per-patient SLR once here; model and generated quantities
  // blocks become trivial single-line loops.
  vector[N] slr_vec;
  for (i in 1:N) {
    real slr = 0.0;

    for (j in 1:J) {
      if (comparable[i, j] == 0 || MOI0[i, j] == 0 || MOIf[i, j] == 0) continue;

      int off = (j - 1) * maxMOI;
      real locus_sum_lr = 0;

      for (af in 1:MOIf[i, j]) {
        int kf = recodedf[i, off + af];
        if (kf == 0 && hiddenf[i, off + af] == 0) continue;

        if (hiddenf[i, off + af] == 1) {
          // Hidden recurrence allele: marginalises to LR = 1
          locus_sum_lr += 1.0;
        } else {
          if (td_has_match[i, off + af]) {
            // Observed Day-0 match: LR = (1 - q_mismatch) / freq[kf]
            locus_sum_lr += exp(l_match - log_freq[j][kf]);
          } else {
            // No observed match: marginalise over hidden Day-0 alleles.
            // BUG 1 FIX: use td_m_i_0 (Day-0 hidden count) not total hidden count.
            real phm = 1.0 - pow(1.0 - freq[j][kf], td_m_i_0[i, j]);
            locus_sum_lr += phm * exp(l_match - log_freq[j][kf])
                          + (1.0 - phm) * exp(l_mismatch);
          }
        }
      }

      if (td_n_alleles_f[i, j] > 0) {
        slr += td_n_lost[i, j]  * l_loss
             + td_m_i[i, j]     * l_dropout
             + log(locus_sum_lr / td_n_alleles_f[i, j] + 1e-10);
      }
    }
    slr_vec[i] = slr;
  }
}
model {
  // Priors match implicit priors of the old R Gibbs sampler posterior updates:
  //   q_mismatch: rbeta(1, 1+n_mismatch, 50+n_match)   → Beta(1, 50)
  //   q_loss:     rbeta(1, 1+n_lost,      1+n_retained) → Beta(1,  1)  [uniform]
  //   q_dropout:  rbeta(1, 1+n_hid,      50+n_obs)      → Beta(1, 50)
  q_mismatch ~ beta(1, 50);
  q_loss     ~ beta(1,  1);
  q_dropout  ~ beta(1, 50);

  for (j in 1:J) {
    vector[max_K] alpha = rep_vector(0.1, max_K);
    for (k in 1:K[j]) alpha[k] = 1.0;
    freq[j]              ~ dirichlet(alpha);
    additional_counts[j] ~ multinomial(freq[j]);
  }

  // theta_recrud = 0.5 fixed -- log_sum_exp(slr, 0) = log(exp(slr)+1), equivalent to theta=0.5
  for (i in 1:N)
    target += log_sum_exp(slr_vec[i], 0.0);
}
generated quantities {
  // Fix #2: slr_vec already in transformed parameters -- trivial here.
  // Fix #4: use inv_logit(slr) -- numerically stable, avoids redundant exp() calls,
  //         consistent with malrebay_model.stan.
  vector[N] p_recrud = inv_logit(slr_vec);
}
