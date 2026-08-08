library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(cmdstanr)
library(posterior)
library(loo)

root_dir <- "C:/Users/Minhua Yan/ASU Dropbox/Minhua Yan/DDRIG/CFpartnerChoice/Latent variable analysis for principles guiding would-suggest divisions"
nopool_dir <- file.path(root_dir, "noPoolingSDRelOnPrincpl")

# ------------------------------------------------------------
# Load data once (files confirmed byte-identical between the two folders)
# ------------------------------------------------------------
dat_raw <- read_excel(file.path(root_dir, "ethics_rel_deid.xlsx"))
emit_raw <- read.csv(file.path(root_dir, "emission_probs_softened_moderate.csv"), sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

dat <- dat_raw %>%
  filter(Name != "Note") %>%
  transmute(
    subject = as.character(Name),
    relationship = as.integer(relationship),
    context = as.character(Contribution_summary),
    response = as.character(best_division_whoMore)
  ) %>%
  mutate(
    context = str_trim(context),
    response = str_trim(response),
    context = if_else(context == "ALLEQ", "ALLEQL", context)
  ) %>%
  filter(response %in% c("CFP", "HD", "SBJ"))

principles <- colnames(emit_raw)[colnames(emit_raw) != "Context"]
responses <- c("SBJ", "HD", "CFP")
C <- nrow(emit_raw); K <- length(principles); R <- length(responses)

emit_array <- array(NA_real_, dim = c(C, K, R), dimnames = list(emit_raw$Context, principles, responses))
for (c_ix in seq_len(C)) for (k_ix in seq_len(K)) {
  probs <- as.numeric(str_split(emit_raw[[principles[k_ix]]][c_ix], ",", simplify = TRUE))
  emit_array[c_ix, k_ix, ] <- probs
}

subject_levels <- sort(unique(dat$subject))
dat <- dat %>%
  mutate(
    subject_id = match(subject, subject_levels),
    context_id = match(context, emit_raw$Context),
    response_id = match(response, responses)
  )

stan_data_base <- list(
  N = nrow(dat), I = length(subject_levels), K = K, C = C, R = R,
  subj = dat$subject_id, ctx = dat$context_id, y = dat$response_id,
  rel = dat$relationship, emit = emit_array, prior_only = 0
)

# ------------------------------------------------------------
# Fit pooled model (root ethics_rel2_model.stan, beta[k] ~ normal(0, sigma_beta))
# ------------------------------------------------------------
mod_pooled <- cmdstan_model(
  file.path(root_dir, "ethics_rel2_model.stan"),
  exe_file = file.path(root_dir, "model_no_random.exe"),
  force_recompile = TRUE
)

fit_pooled <- mod_pooled$sample(
  data = stan_data_base,
  seed = 1234, chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 2000, refresh = 0
)

# ------------------------------------------------------------
# Fit no-pooling model, beta_prior_scale = 2 (the scale used in
# run_stan_kinVSnon.R, and the one the sensitivity check supported)
# ------------------------------------------------------------
mod_nopool <- cmdstan_model(
  file.path(nopool_dir, "ethics_rel2_model.stan"),
  exe_file = file.path(nopool_dir, "model_no_random.exe"),
  force_recompile = TRUE
)

stan_data_nopool <- c(stan_data_base, list(beta_prior_scale = 2))

fit_nopool <- mod_nopool$sample(
  data = stan_data_nopool,
  seed = 1234, chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 2000, refresh = 0
)

# ------------------------------------------------------------
# LOO comparison
# ------------------------------------------------------------
loo_pooled <- fit_pooled$loo()
loo_nopool <- fit_nopool$loo()

comp <- loo_compare(list(pooled = loo_pooled, no_pooling_scale2 = loo_nopool))

out_dir <- file.path(root_dir, "output")
dir.create(out_dir, showWarnings = FALSE)

capture.output({
  cat("=== Pooled model (sigma_beta shared across principles) LOO ===\n")
  print(loo_pooled)
  cat("\n=== No-pooling model (beta_prior_scale = 2) LOO ===\n")
  print(loo_nopool)
  cat("\n=== loo_compare (positive elpd_diff favors the top model) ===\n")
  print(comp)
}, file = file.path(out_dir, "loo_comparison.txt"))

cat("\n\n===== SUMMARY =====\n")
print(loo_pooled)
cat("\n")
print(loo_nopool)
cat("\n")
print(comp)
cat("\nSaved output/loo_comparison.txt\n")
