## Prior sensitivity check for the no-pooling model's beta ~ normal(0, scale)
## prior. That scale (2) was calibrated by analogy for a different contrast
## (food-vs-cash crop type, in a smaller subsample) - see
## ethics_rel2_model.stan's comment. This script refits the same model at a
## few alternative scales to see how much the headline non-kin/non-friends
## tradition/reward-resource/generosity numbers actually depend on it.
##
## Does not modify run_stan_kinVSnon.R's canonical scale = 2 run; writes each
## scale's output to its own subfolder instead of a separate top-level
## project folder, since these runs are disposable diagnostics, not a
## standing model variant.

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(cmdstanr)
library(posterior)
library(here)
library(ggplot2)

# ------------------------------------------------------------
# 0. Compile Stan model (same model file as run_stan_kinVSnon.R;
# beta_prior_scale is data, not hardcoded, so one compile serves all scales)
# ------------------------------------------------------------
mod <- cmdstan_model("ethics_rel2_model.stan",
                     exe_file = file.path(getwd(), "model_no_random.exe"),
                     force_recompile = TRUE)

# ------------------------------------------------------------
# 1. Load and prep data (identical to run_stan_kinVSnon.R)
# ------------------------------------------------------------

dat_raw <- read_excel(here("ethics_rel_deid.xlsx"))

emit_raw <- read.csv(
  here("emission_probs_softened_moderate.csv"),
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

dat <- dat_raw %>%
  filter(Name != "Note") %>%
  transmute(
    subject = as.character(Name),
    relationship = as.integer(relationship),
    context = as.character(Contribution_summary),
    response = as.character(best_division_whoMore)
  )

dat <- dat %>%
  mutate(
    context = str_trim(context),
    response = str_trim(response),
    context = if_else(context == "ALLEQ", "ALLEQL", context)
  )

valid_responses <- c("CFP", "HD", "SBJ")
dat <- dat %>% filter(response %in% valid_responses)

if (!all(dat$relationship %in% c(-1, 1))) {
  stop("relationship must be coded as -1 or 1")
}

emit_contexts <- emit_raw$Context
data_contexts <- unique(dat$context)
missing_contexts <- setdiff(data_contexts, emit_contexts)
if (length(missing_contexts) > 0) {
  stop("Some contexts in the data are not in the emission table.")
}

principles <- colnames(emit_raw)[colnames(emit_raw) != "Context"]
# see run_stan_kinVSnon.R: emit_raw's columns are authored in (SBJ, HD, CFP)
# order, not (CFP, HD, SBJ)
responses <- c("SBJ", "HD", "CFP")

C <- nrow(emit_raw)
K <- length(principles)
R <- length(responses)

emit_array <- array(
  NA_real_,
  dim = c(C, K, R),
  dimnames = list(emit_raw$Context, principles, responses)
)

for (c_ix in seq_len(C)) {
  for (k_ix in seq_len(K)) {
    cell <- emit_raw[[principles[k_ix]]][c_ix]
    probs <- as.numeric(str_split(cell, ",", simplify = TRUE))
    emit_array[c_ix, k_ix, ] <- probs
  }
}

subject_levels <- sort(unique(dat$subject))
context_levels <- emit_raw$Context
response_levels <- responses

dat <- dat %>%
  mutate(
    subject_id = match(subject, subject_levels),
    context_id = match(context, context_levels),
    response_id = match(response, response_levels)
  )

stan_data_base <- list(
  N = nrow(dat),
  I = length(subject_levels),
  K = K,
  C = C,
  R = R,
  subj = dat$subject_id,
  ctx = dat$context_id,
  y = dat$response_id,
  rel = dat$relationship,
  emit = emit_array
)

# ------------------------------------------------------------
# 2. Labels (identical to run_stan_kinVSnon.R)
# ------------------------------------------------------------

principle_names <- c(
  "selfish", "generous", "tradition", "rewardLBonly", "rewardLDonly",
  "rewardRSonly", "rewardLB_LD", "rewardLB_RS", "rewardLD_RS",
  "rewardALL"
)

pretty_labels <- c(
  "selfish" = "selfishness",
  "generous" = "generosity",
  "tradition" = "tradition\n compliance",
  "rewardLBonly" = "reward\n labor",
  "rewardLDonly" = "reward\n land",
  "rewardRSonly" = "reward\n resource",
  "rewardLB_LD" = "reward labor\n and land",
  "rewardLB_RS" = "reward labor\n and resource",
  "rewardLD_RS" = "reward land\n and resource",
  "rewardALL" = "reward\n all three"
)

desired_order <- c(
  "generosity",
  "tradition\n compliance",
  "reward\n labor",
  "reward\n land",
  "reward\n resource",
  "reward labor\n and land",
  "reward labor\n and resource",
  "reward land\n and resource",
  "reward\n all three",
  "selfishness"
)

# ------------------------------------------------------------
# 3. Loop over alternative beta prior scales
# ------------------------------------------------------------

scale_values <- c(1, 2, 3)
headline_principles <- c("generosity", "tradition\n compliance", "reward\n resource")

for (scale_value in scale_values) {

  out_dir <- file.path("output", paste0("scale_", scale_value))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  stan_data_post <- c(stan_data_base, list(prior_only = 0, beta_prior_scale = scale_value))
  stan_data_prior <- stan_data_post
  stan_data_prior$prior_only <- 1

  fit_prior <- mod$sample(
    data = stan_data_prior,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 0
  )

  fit <- mod$sample(
    data = stan_data_post,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    refresh = 0
  )

  param_summary <- fit$summary(c("mu_alpha", "sigma_alpha", "beta"))
  write.csv(param_summary, file.path(out_dir, "param_summary.csv"), row.names = FALSE)

  prior_draws_rel_minus1 <- as_draws_matrix(fit_prior$draws("pop_prob_rel_minus1"))
  prior_draws_rel_plus1  <- as_draws_matrix(fit_prior$draws("pop_prob_rel_plus1"))
  colnames(prior_draws_rel_minus1) <- principle_names
  colnames(prior_draws_rel_plus1)  <- principle_names

  post_draws_rel_minus1 <- as_draws_matrix(fit$draws("pop_prob_rel_minus1"))
  post_draws_rel_plus1  <- as_draws_matrix(fit$draws("pop_prob_rel_plus1"))
  colnames(post_draws_rel_minus1) <- principle_names
  colnames(post_draws_rel_plus1)  <- principle_names

  summary_post <- bind_rows(
    data.frame(
      ethical_principle = principle_names,
      mean = colMeans(post_draws_rel_minus1),
      ci_lower = apply(post_draws_rel_minus1, 2, quantile, 0.05),
      ci_upper = apply(post_draws_rel_minus1, 2, quantile, 0.95),
      relationship = "kin/friends"
    ),
    data.frame(
      ethical_principle = principle_names,
      mean = colMeans(post_draws_rel_plus1),
      ci_lower = apply(post_draws_rel_plus1, 2, quantile, 0.05),
      ci_upper = apply(post_draws_rel_plus1, 2, quantile, 0.95),
      relationship = "non-kin/non-friends"
    )
  ) %>%
    mutate(
      ethical_principle = recode(ethical_principle, !!!pretty_labels),
      ethical_principle = factor(ethical_principle, levels = desired_order),
      relationship = factor(relationship, levels = c("kin/friends", "non-kin/non-friends"))
    )

  write.csv(summary_post, file.path(out_dir, "summary_post.csv"), row.names = FALSE)

  diff_draws <- post_draws_rel_plus1 - post_draws_rel_minus1
  summary_diff <- data.frame(
    ethical_principle = principle_names,
    mean = colMeans(diff_draws),
    ci_lower = apply(diff_draws, 2, quantile, 0.05),
    ci_upper = apply(diff_draws, 2, quantile, 0.95),
    prob_positive = apply(diff_draws, 2, function(x) mean(x > 0)),
    prob_negative = apply(diff_draws, 2, function(x) mean(x < 0))
  ) %>%
    mutate(
      ethical_principle = recode(ethical_principle, !!!pretty_labels),
      ethical_principle = factor(ethical_principle, levels = desired_order)
    )
  write.csv(summary_diff, file.path(out_dir, "summary_diff.csv"), row.names = FALSE)

  prior_diff_draws <- prior_draws_rel_plus1 - prior_draws_rel_minus1
  summary_prior_diff <- data.frame(
    ethical_principle = principle_names,
    mean = colMeans(prior_diff_draws),
    ci_lower = apply(prior_diff_draws, 2, quantile, 0.05),
    ci_upper = apply(prior_diff_draws, 2, quantile, 0.95),
    prob_positive = apply(prior_diff_draws, 2, function(x) mean(x > 0)),
    prob_negative = apply(prior_diff_draws, 2, function(x) mean(x < 0))
  ) %>%
    mutate(
      ethical_principle = recode(ethical_principle, !!!pretty_labels),
      ethical_principle = factor(ethical_principle, levels = desired_order)
    )
  write.csv(summary_prior_diff, file.path(out_dir, "summary_prior_diff.csv"), row.names = FALSE)

  p_post <- ggplot(
    summary_post,
    aes(x = ethical_principle, y = mean, color = relationship)
  ) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbar(
      aes(ymin = ci_lower, ymax = ci_upper),
      width = 0.15,
      position = position_dodge(width = 0.5)
    ) +
    labs(
      title = paste0("Posterior Distribution by Relationship (beta scale = ", scale_value, ")"),
      y = "Posterior probability",
      x = NULL,
      color = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(hjust = 0.5),
      legend.position = c(0.8, 0.8),
      text = element_text(size = 14)
    )
  ggsave(file.path(out_dir, "plot_post.png"), p_post, width = 10.5, height = 4.5, dpi = 150)

  cat("\n=== scale =", scale_value, "- non-kin/non-friends headline principles ===\n")
  print(
    summary_post %>%
      filter(relationship == "non-kin/non-friends", ethical_principle %in% headline_principles) %>%
      select(ethical_principle, mean, ci_lower, ci_upper)
  )
}

cat("\nDone. Compare output/scale_1, output/scale_2, output/scale_3.\n")
