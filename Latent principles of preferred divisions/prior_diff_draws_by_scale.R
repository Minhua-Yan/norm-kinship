library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(cmdstanr)
library(posterior)
library(here)
library(ggplot2)

setwd("C:/Users/Minhua Yan/ASU Dropbox/Minhua Yan/DDRIG/CFpartnerChoice/Latent variable analysis for principles guiding would-suggest divisions/noPoolingSDRelOnPrincpl")

mod <- cmdstan_model("ethics_rel2_model.stan",
                     exe_file = file.path(getwd(), "model_no_random.exe"),
                     force_recompile = TRUE)

dat_raw <- read_excel(here("ethics_rel_deid.xlsx"))
emit_raw <- read.csv(here("emission_probs_softened_moderate.csv"), sep = ";", stringsAsFactors = FALSE, check.names = FALSE)

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
  rel = dat$relationship, emit = emit_array, prior_only = 1
)

principle_names <- c(
  "selfish", "generous", "tradition", "rewardLBonly", "rewardLDonly",
  "rewardRSonly", "rewardLB_LD", "rewardLB_RS", "rewardLD_RS", "rewardALL"
)
pretty_labels <- c(
  "selfish" = "selfishness", "generous" = "generosity",
  "tradition" = "tradition\n compliance", "rewardLBonly" = "reward\n labor",
  "rewardLDonly" = "reward\n land", "rewardRSonly" = "reward\n resource",
  "rewardLB_LD" = "reward labor\n and land", "rewardLB_RS" = "reward labor\n and resource",
  "rewardLD_RS" = "reward land\n and resource", "rewardALL" = "reward\n all three"
)
desired_order <- c(
  "generosity", "tradition\n compliance", "reward\n labor", "reward\n land",
  "reward\n resource", "reward labor\n and land", "reward labor\n and resource",
  "reward land\n and resource", "reward\n all three", "selfishness"
)

scale_values <- c(1, 2, 3)

all_draws <- bind_rows(lapply(scale_values, function(s) {
  stan_data_prior <- c(stan_data_base, list(beta_prior_scale = s))
  fit_prior <- mod$sample(
    data = stan_data_prior, seed = 1234, chains = 4, parallel_chains = 4,
    iter_warmup = 1000, iter_sampling = 1000, refresh = 0
  )
  d_minus1 <- as_draws_matrix(fit_prior$draws("pop_prob_rel_minus1"))
  d_plus1  <- as_draws_matrix(fit_prior$draws("pop_prob_rel_plus1"))
  colnames(d_minus1) <- principle_names
  colnames(d_plus1)  <- principle_names
  diff_mat <- d_plus1 - d_minus1
  df <- as.data.frame(diff_mat)
  df$scale <- s
  df
}))

long_draws <- all_draws %>%
  pivot_longer(cols = all_of(principle_names), names_to = "ethical_principle", values_to = "diff") %>%
  mutate(
    ethical_principle = recode(ethical_principle, !!!pretty_labels),
    ethical_principle = factor(ethical_principle, levels = desired_order),
    scale = factor(paste0("scale = ", scale), levels = paste0("scale = ", scale_values))
  )

p <- ggplot(long_draws, aes(x = ethical_principle, y = diff, color = scale)) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray30") +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6),
    size = 0.5, alpha = 0.04
  ) +
  labs(
    title = "Prior Difference Draws in Ethical Principle Probabilities, by beta_prior_scale",
    subtitle = "non-kin/non-friends minus kin/friends; each point is one prior draw",
    x = NULL,
    y = "Difference in prior probability",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    text = element_text(size = 13)
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))

ggsave("output/plot_prior_diff_draws_by_scale.png", p, width = 13, height = 5.5, dpi = 150)
write.csv(long_draws, "output/prior_diff_draws_by_scale.csv", row.names = FALSE)
cat("Saved output/plot_prior_diff_draws_by_scale.png\n")
