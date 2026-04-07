## how ethical principles governing best division responses shift based on relationship (social kin/friends or genetic kin)

library(readxl)
library(dplyr)
library(stringr)
library(cmdstanr)
library(posterior)
library(here)
library(ggplot2)

# ------------------------------------------------------------
# 0. Compile Stan model
# ------------------------------------------------------------

mod <- cmdstan_model("ethics_rel2_model.stan", force_recompile = TRUE)

# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

dat_raw <- read_excel(here("ethics_rel_deid.xlsx"))

emit_raw <- read.csv(
  here("emission_probs_softened_moderate.csv"),
  sep = ";",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# 2. Keep variables used in the model
# ------------------------------------------------------------

dat <- dat_raw %>%
  filter(Name != "Note") %>%
  transmute(
    subject = as.character(Name),
    relationship_kin = relationship_kin3,
    context = as.character(Contribution_summary),
    response = as.character(best_division_whoMore)
  )

# ------------------------------------------------------------
# 3. Clean strings and drop rows with missing relationship_kin
# genetic kin = -1
# social kin = 1
# ------------------------------------------------------------

dat <- dat %>%
  mutate(
    context = str_trim(context),
    response = str_trim(response),
    context = if_else(context == "ALLEQ", "ALLEQL", context),
    relationship_kin = as.numeric(relationship_kin)
  ) %>%
  filter(!is.na(relationship_kin))

valid_responses <- c("CFP", "HD", "SBJ")

dat <- dat %>%
  filter(response %in% valid_responses)

if (!all(dat$relationship_kin %in% c(-1, 1))) {
  stop("relationship_kin must be coded as -1 or 1 after removing NAs")
}

# ------------------------------------------------------------
# 4. Check contexts match emission table
# ------------------------------------------------------------

emit_contexts <- emit_raw$Context
data_contexts <- unique(dat$context)

missing_contexts <- setdiff(data_contexts, emit_contexts)

if (length(missing_contexts) > 0) {
  stop("Some contexts in the data are not in the emission table.")
}

# ------------------------------------------------------------
# 5. Parse emission probability table
# ------------------------------------------------------------

principles <- colnames(emit_raw)[colnames(emit_raw) != "Context"]
responses <- c("CFP", "HD", "SBJ")

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

# ------------------------------------------------------------
# 6. Convert variables to integer indices for Stan
# ------------------------------------------------------------

subject_levels <- sort(unique(dat$subject))
context_levels <- emit_raw$Context
response_levels <- responses

dat <- dat %>%
  mutate(
    subject_id = match(subject, subject_levels),
    context_id = match(context, context_levels),
    response_id = match(response, response_levels)
  )

# ------------------------------------------------------------
# 7. Build Stan data lists
# ------------------------------------------------------------

stan_data_post <- list(
  N = nrow(dat),
  I = length(subject_levels),
  K = K,
  C = C,
  R = R,
  subj = dat$subject_id,
  ctx = dat$context_id,
  y = dat$response_id,
  rel = dat$relationship_kin,
  emit = emit_array,
  prior_only = 0
)

stan_data_prior <- stan_data_post
stan_data_prior$prior_only <- 1

# ------------------------------------------------------------
# 8. Fit PRIOR
# ------------------------------------------------------------

fit_prior <- mod$sample(
  data = stan_data_prior,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 1000,
  refresh = 200
)

# ------------------------------------------------------------
# 9. Fit POSTERIOR
# ------------------------------------------------------------

fit <- mod$sample(
  data = stan_data_post,
  seed = 1234,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1500,
  iter_sampling = 1500,
  refresh = 200
)

# ------------------------------------------------------------
# 10. Parameter summaries
# ------------------------------------------------------------

print(fit$summary(c("mu_alpha", "sigma_alpha", "beta", "sigma_beta")))

# ------------------------------------------------------------
# 11. Labels for plotting
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
# 12. PRIOR plot
# genetic kin = -1
# social kin = 1
# ------------------------------------------------------------

prior_draws_rel_minus1 <- as_draws_matrix(fit_prior$draws("pop_prob_rel_minus1"))
prior_draws_rel_plus1  <- as_draws_matrix(fit_prior$draws("pop_prob_rel_plus1"))

colnames(prior_draws_rel_minus1) <- principle_names
colnames(prior_draws_rel_plus1)  <- principle_names

summary_prior_rel_minus1 <- data.frame(
  ethical_principle = principle_names,
  mean = colMeans(prior_draws_rel_minus1),
  ci_lower = apply(prior_draws_rel_minus1, 2, quantile, 0.025),
  ci_upper = apply(prior_draws_rel_minus1, 2, quantile, 0.975),
  relationship = "genetic kin"
)

summary_prior_rel_plus1 <- data.frame(
  ethical_principle = principle_names,
  mean = colMeans(prior_draws_rel_plus1),
  ci_lower = apply(prior_draws_rel_plus1, 2, quantile, 0.025),
  ci_upper = apply(prior_draws_rel_plus1, 2, quantile, 0.975),
  relationship = "social kin"
)

summary_prior <- bind_rows(summary_prior_rel_minus1, summary_prior_rel_plus1) %>%
  mutate(
    ethical_principle = recode(ethical_principle, !!!pretty_labels),
    ethical_principle = factor(ethical_principle, levels = desired_order),
    relationship = factor(relationship, levels = c("genetic kin", "social kin"))
  )

p_prior <- ggplot(
  summary_prior,
  aes(x = ethical_principle, y = mean, color = relationship)
) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    width = 0.15,
    position = position_dodge(width = 0.5)
  ) +
  labs(
    title = "Prior Distribution of Ethical Principles by Kin Type",
    y = "Prior probability",
    x = NULL,
    color = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

print(p_prior)

# ------------------------------------------------------------
# 13. POSTERIOR plot
# genetic kin = -1
# social kin = 1
# ------------------------------------------------------------

post_draws_rel_minus1 <- as_draws_matrix(fit$draws("pop_prob_rel_minus1"))
post_draws_rel_plus1  <- as_draws_matrix(fit$draws("pop_prob_rel_plus1"))

colnames(post_draws_rel_minus1) <- principle_names
colnames(post_draws_rel_plus1)  <- principle_names

summary_post_rel_minus1 <- data.frame(
  ethical_principle = principle_names,
  mean = colMeans(post_draws_rel_minus1),
  ci_lower = apply(post_draws_rel_minus1, 2, quantile, 0.05),
  ci_upper = apply(post_draws_rel_minus1, 2, quantile, 0.95),
  relationship = "genetic kin"
)

summary_post_rel_plus1 <- data.frame(
  ethical_principle = principle_names,
  mean = colMeans(post_draws_rel_plus1),
  ci_lower = apply(post_draws_rel_plus1, 2, quantile, 0.05),
  ci_upper = apply(post_draws_rel_plus1, 2, quantile, 0.95),
  relationship = "social kin"
)

summary_post <- bind_rows(summary_post_rel_minus1, summary_post_rel_plus1) %>%
  mutate(
    ethical_principle = recode(ethical_principle, !!!pretty_labels),
    ethical_principle = factor(ethical_principle, levels = desired_order),
    relationship = factor(relationship, levels = c("genetic kin", "social kin"))
  )

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
    title = "Posterior Distribution of Ethical Principles by Kin Type",
    y = "Posterior probability",
    x = NULL,
    color = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

print(p_post)

# ------------------------------------------------------------
# 14. Posterior comparison plot for all ethical principles
# difference = social kin - genetic kin
# ------------------------------------------------------------

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

p_diff <- ggplot(summary_diff, aes(x = ethical_principle, y = mean)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.15) +
  labs(
    title = "Posterior Difference in Ethical Principle Probabilities",
    subtitle = "social kin minus genetic kin",
    x = NULL,
    y = "Difference in posterior probability"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

print(p_diff)
print(summary_diff)

# ------------------------------------------------------------
# 15. Relationship coverage check
# ------------------------------------------------------------

coverage <- dat %>%
  distinct(subject, relationship_kin) %>%
  count(subject)

table(coverage$n)

# ------------------------------------------------------------
# 16. Count observations by kin type
# ------------------------------------------------------------

table(dat$relationship_kin, useNA = "always")

