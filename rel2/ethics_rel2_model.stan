// this stan is used when relationship has two values, -1 and 1 

data {
  int<lower=1> N;
  int<lower=1> I;
  int<lower=2> K;
  int<lower=1> C;
  int<lower=2> R;

  array[N] int<lower=1, upper=I> subj;
  array[N] int<lower=1, upper=C> ctx;
  array[N] int<lower=1, upper=R> y;

  // relationship coding:
  // -1 = kin/friends in kinVSnon and genetic kin in socialVSgenetic
  //  1 = non-kin/non-friends in kinVSnon and social kin in socialVSgenetic
  array[N] int<lower=-1, upper=1> rel;

  array[C, K, R] real<lower=0, upper=1> emit;

  // 0 = posterior fit
  // 1 = prior-only fit
  int<lower=0, upper=1> prior_only;
}

transformed data {
  // Prior target: uniform across all 10 ethical principles
  // Since K = 10, this centers the prior around 0.1 for each principle

  vector[K] base_logits;

  // Equal logits imply equal prior probabilities after centering + softmax
  base_logits = rep_vector(0, K);
}

parameters {
  matrix[I, K] alpha_raw;
  vector[K] mu_alpha;
  vector<lower=0>[K] sigma_alpha;
  vector[K] beta;
  real<lower=0> sigma_beta;
}

transformed parameters {
  matrix[I, K] alpha;

  for (i in 1:I) {
    for (k in 1:K) {
      alpha[i, k] = mu_alpha[k] + sigma_alpha[k] * alpha_raw[i, k];
    }
  }
}

model {
  mu_alpha ~ normal(base_logits, 0.8);
  sigma_alpha ~ normal(0, 1);
  to_vector(alpha_raw) ~ normal(0, 1);

  beta ~ normal(0, sigma_beta);
  sigma_beta ~ normal(0, 1);

  if (prior_only == 0) {
    for (n in 1:N) {
      vector[K] eta;
      vector[K] eta_centered;
      vector[K] log_pi;
      vector[K] log_mix_terms;

      for (k in 1:K) {
        eta[k] = alpha[subj[n], k] + beta[k] * rel[n];
      }

      eta_centered = eta - mean(eta);
      log_pi = log_softmax(eta_centered);

      for (k in 1:K) {
        log_mix_terms[k] = log_pi[k] + log(emit[ctx[n], k, y[n]]);
      }

      target += log_sum_exp(log_mix_terms);
    }
  }
}

generated quantities {
  vector[N] log_lik;
  array[N] int y_rep;

  // relationship = -1  -> kin/friends or genetic kin
  vector[K] mu_alpha_rel_minus1_centered;
  vector[K] pop_prob_rel_minus1;

  // relationship = 1   -> non-kin/non-friends or social kin
  vector[K] mu_alpha_rel_plus1_centered;
  vector[K] pop_prob_rel_plus1;

  mu_alpha_rel_minus1_centered = (mu_alpha - beta) - mean(mu_alpha - beta);
  pop_prob_rel_minus1 = softmax(mu_alpha_rel_minus1_centered);

  mu_alpha_rel_plus1_centered = (mu_alpha + beta) - mean(mu_alpha + beta);
  pop_prob_rel_plus1 = softmax(mu_alpha_rel_plus1_centered);

  for (n in 1:N) {
    vector[K] eta;
    vector[K] eta_centered;
    vector[K] log_pi;
    vector[K] log_mix_terms;
    vector[R] p_y;

    for (k in 1:K) {
      eta[k] = alpha[subj[n], k] + beta[k] * rel[n];
    }

    eta_centered = eta - mean(eta);
    log_pi = log_softmax(eta_centered);

    if (prior_only == 0) {
      for (k in 1:K) {
        log_mix_terms[k] = log_pi[k] + log(emit[ctx[n], k, y[n]]);
      }
      log_lik[n] = log_sum_exp(log_mix_terms);

      for (r in 1:R) {
        real total_prob;
        total_prob = 0;
        for (k in 1:K) {
          total_prob += exp(log_pi[k]) * emit[ctx[n], k, r];
        }
        p_y[r] = total_prob;
      }

      p_y = p_y / sum(p_y);
      y_rep[n] = categorical_rng(p_y);
    } else {
      log_lik[n] = 0;
      y_rep[n] = 1;
    }
  }
}
