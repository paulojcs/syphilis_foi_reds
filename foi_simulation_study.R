# =============================================================================
# SIMULATION STUDY: VALIDATING RANDOM WALK FOI MODEL (MULTI-YEAR)
# =============================================================================
# Purpose: Test if the model can recover known FOI parameters from simulated
#          multi-year serosurvey data (matching the real data structure)
#
# Workflow:
#   1. Define TRUE FOI trajectory
#   2. Simulate seroprevalence data for multiple survey years x ages
#   3. Fit the random walk model to simulated data
#   4. Compare estimated FOI vs true FOI
# =============================================================================

rm(list = ls())
set.seed(42)

# =============================================================================
# SECTION 1: SETUP
# =============================================================================

required_packages <- c("tidyverse", "rstan", "bayesplot", "patchwork")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# Stan configuration (adjust cores for your machine)
options(mc.cores = parallel::detectCores() - 1)
rstan_options(auto_write = TRUE)

cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("FOI MODEL SIMULATION STUDY (MULTI-YEAR)\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

# =============================================================================
# SECTION 2: DEFINE SIMULATION PARAMETERS
# =============================================================================

# --- Multi-year design matching real data ---
survey_years <- 2008:2023
min_age <- 18
max_age <- 65
min_birth_year <- min(survey_years) - max_age  # 1943
n_years <- max(survey_years) - 1 - min_birth_year + 1  # 80
years <- seq(min_birth_year, max(survey_years) - 1)     # 1943:2022

cat(sprintf("Survey years: %d to %d\n", min(survey_years), max(survey_years)))
cat(sprintf("FOI years: %d to %d (%d years)\n", min(years), max(years), n_years))
cat(sprintf("Age range: %d to %d\n", min_age, max_age))

# --- Sample sizes (per year-age cell, smaller for local testing) ---
n_per_age <- 200  # Adjust based on your machine (200-500)

# =============================================================================
# SECTION 3: DEFINE TRUE FOI TRAJECTORY
# =============================================================================
#
# We'll test with a realistic scenario:
# - Low historical FOI (1943-2000)
# - Gradual increase (2001-2015)
# - Higher recent FOI (2016-2022)
#
# This mimics the syphilis resurgence pattern observed in Brazil
# =============================================================================

cat("\n--- Defining TRUE FOI trajectory ---\n")

# -----------------------------------------------------------------------
# NOTE: All FOI values calibrated so that average seroprevalence across
# ages 18-65 and survey years 2008-2023 is approximately 3%, matching
# the observed prevalence in the REDS blood donor data.
# -----------------------------------------------------------------------

# Option 1: Piecewise pattern (easier to verify)
true_foi_piecewise <- case_when(
  years <= 2000 ~ 0.0003, # 0.3 per 1000
  years <= 2010 ~ 0.0006, # 0.6 per 1000
  years <= 2015 ~ 0.0010, # 1.0 per 1000
  TRUE ~ 0.0015           # 1.5 per 1000
)

# Option 2: Smooth increasing pattern
# Logistic growth from ~0.3 to ~1.5 per 1000
true_foi_smooth <- 0.0003 + 0.0012 * plogis((years - 2005) / 5)

# Option 3: Constant (simplest test - should recover well)
true_foi_constant <- rep(0.0007, n_years) # 0.7 per 1000

# Option 4: Epidemic waves (realistic syphilis pattern)
# - Low baseline in the 1940s-50s (~0.3 per 1000)
# - Big increase peaking in the mid-1970s (~1.8 per 1000)
# - Slow decline through 1980s (~1.2-1.7), gradual fall in 1990s
# - Trough late 1990s-early 2000s
# - Resurgence from ~2005 onwards (~1.3 per 1000 by 2022)
# Asymmetric Gaussian: SD=7 on the rise, SD=12 on the decline
wave1_sd <- ifelse(years <= 1975, 7, 12)
true_foi_epidemic <- 0.0003 +
  0.0015 * exp(-0.5 * ((years - 1975) / wave1_sd)^2) +  # asymmetric 1970s peak
  0.0010 * plogis((years - 2012) / 3)                     # 2005+ logistic resurgence

# --- SELECT WHICH SCENARIO TO USE ---
# Change this to test different scenarios
scenario <- "epidemic" # Options: "piecewise", "smooth", "constant", "epidemic"

true_foi <- switch(
  scenario,
  "piecewise" = true_foi_piecewise,
  "smooth" = true_foi_smooth,
  "constant" = true_foi_constant,
  "epidemic" = true_foi_epidemic
)

cat(sprintf("Scenario: %s\n", scenario))
cat(sprintf(
  "FOI range: %.1f to %.1f per 1000\n",
  min(true_foi) * 1000,
  max(true_foi) * 1000
))

# Visualize true FOI
true_foi_df <- data.frame(
  year = years,
  foi = true_foi,
  foi_per_1000 = true_foi * 1000
)

p_true <- ggplot(true_foi_df, aes(x = year, y = foi_per_1000)) +
  geom_line(color = "darkred", linewidth = 1.2) +
  geom_point(color = "darkred", size = 0.8) +
  labs(
    title = sprintf("TRUE Force of Infection (%s scenario)", scenario),
    x = "Calendar Year",
    y = "FOI (per 1,000 person-years)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

print(p_true)

# =============================================================================
# SECTION 4: SIMULATE MULTI-YEAR SEROPREVALENCE DATA
# =============================================================================

cat("\n--- Simulating multi-year seroprevalence data ---\n")

# Function to calculate P(seropositive) given age, survey year, and FOI history
calc_seroprev <- function(age, survey_year, foi_vector, min_year) {
  birth_year <- survey_year - age
  prob_neg <- 1.0

  for (y in 1:age) {
    cal_year <- birth_year + y - 1
    year_idx <- cal_year - min_year + 1

    if (year_idx >= 1 && year_idx <= length(foi_vector)) {
      prob_neg <- prob_neg * (1 - foi_vector[year_idx])
    }
  }

  return(1 - prob_neg)
}

# Generate multi-year data: expand.grid over survey_years x ages
sim_data <- expand.grid(
  survey_year = survey_years,
  age = min_age:max_age
) %>%
  as.data.frame() %>%
  rowwise() %>%
  mutate(
    # True probability of seropositivity at this age in this survey year
    true_prev = calc_seroprev(age, survey_year, true_foi, min_birth_year),

    # Sample size per cell
    n_sample = n_per_age,

    # Simulate observed positives (binomial)
    n_pos = rbinom(1, n_sample, true_prev),

    # Observed prevalence
    obs_prev = n_pos / n_sample
  ) %>%
  ungroup()

cat(sprintf("Simulated %d observations (year x age)\n", nrow(sim_data)))
cat(sprintf("Survey years: %d\n", length(unique(sim_data$survey_year))))
cat(sprintf("Age groups: %d\n", length(unique(sim_data$age))))
cat(sprintf(
  "Total sample size: %s\n",
  format(sum(sim_data$n_sample), big.mark = ",")
))
cat(sprintf(
  "Total positives: %s\n",
  format(sum(sim_data$n_pos), big.mark = ",")
))
cat(sprintf(
  "Overall prevalence: %.2f%%\n",
  sum(sim_data$n_pos) / sum(sim_data$n_sample) * 100
))

# Plot simulated data for a subset of years
show_years <- c(2008, 2012, 2016, 2020, 2023)

p_data <- ggplot(sim_data %>% filter(survey_year %in% show_years), aes(x = age)) +
  geom_line(
    aes(y = true_prev * 100),
    color = "darkred",
    linewidth = 0.8,
    linetype = "dashed",
    alpha = 0.7
  ) +
  geom_point(
    aes(y = obs_prev * 100),
    color = "steelblue",
    size = 1.5,
    alpha = 0.5
  ) +
  facet_wrap(~survey_year, ncol = 3) +
  labs(
    title = "Simulated Age-Seroprevalence Data",
    subtitle = "Points: observed | Dashed line: true expected",
    x = "Age (years)",
    y = "Seroprevalence (%)"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

print(p_data)

# =============================================================================
# SECTION 5: PREPARE STAN DATA
# =============================================================================

cat("\n--- Preparing Stan data ---\n")

stan_data <- list(
  N = nrow(sim_data),
  n_sample = sim_data$n_sample,
  n_pos = sim_data$n_pos,
  age = sim_data$age,
  survey_year = sim_data$survey_year,  # VECTOR: per observation
  n_years = n_years,
  min_year = min_birth_year,

  # Model settings
  include_seroreversion = 0L,

  # Priors (calibrated for ~3% seroprevalence; exp(-7) ≈ 0.0009)
  log_foi_mean_prior = -7.0,
  log_foi_sd_prior = 1.5,
  sigma_rw_upper = 0.5
)

cat("Stan data prepared:\n")
cat(sprintf("  N (observations): %d\n", stan_data$N))
cat(sprintf("  n_years: %d\n", stan_data$n_years))
cat(sprintf("  min_year: %d\n", stan_data$min_year))
cat(sprintf("  survey_year: vector of length %d\n", length(stan_data$survey_year)))

# =============================================================================
# SECTION 6: LOAD STAN MODEL FROM FILE
# =============================================================================
# Using the actual .stan file ensures we test exactly what runs in production

# --- Auto-detect HPC vs local path ---
if (dir.exists("/scratch/7631403/reds_data/")) {
  stan_file <- "/scratch/7631403/reds_data/stan_models/foi_random_walk_v2.stan"
} else {
  stan_file <- "C:/Users/paulo/OneDrive/Documentos/AnalistaZen/REDS/transfer/new/stan_models/foi_random_walk_v2.stan"
}

if (!file.exists(stan_file)) {
  stop("Stan file not found: ", stan_file)
}

# Compile model from file
cat("\n--- Compiling Stan model ---\n")
cat(sprintf("Using: %s\n", stan_file))
model <- stan_model(file = stan_file)
cat("Model compiled\n")

# =============================================================================
# SECTION 7: FIT MODEL TO SIMULATED DATA
# =============================================================================

cat("\n--- Fitting model ---\n")

fit <- sampling(
  model,
  data = stan_data,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  seed = 42,
  control = list(adapt_delta = 0.97, max_treedepth = 13),
  refresh = 500 # Print progress every 500 iterations
)

cat("\nModel fitted\n")

# =============================================================================
# SECTION 8: DIAGNOSTICS
# =============================================================================

cat("\n--- MCMC Diagnostics ---\n")

# Check Rhat and ESS
fit_summary <- summary(fit, pars = c("sigma_rw", "log_foi_mean"))$summary
print(fit_summary)

# Check for divergences
sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
cat(sprintf("\nDivergent transitions: %d\n", n_divergent))

# Max Rhat
all_summary <- summary(fit)$summary
max_rhat <- max(all_summary[, "Rhat"], na.rm = TRUE)
cat(sprintf("Max Rhat: %.4f (should be < 1.01)\n", max_rhat))

# =============================================================================
# SECTION 9: COMPARE ESTIMATED VS TRUE FOI
# =============================================================================

cat("\n--- Comparing Estimated vs True FOI ---\n")

# Extract posterior
posterior <- rstan::extract(fit)

# Create comparison dataframe
foi_comparison <- data.frame(
  year = years,
  true_foi = true_foi * 1000,
  est_mean = apply(posterior$foi, 2, mean) * 1000,
  est_median = apply(posterior$foi, 2, median) * 1000,
  est_lower_95 = apply(posterior$foi, 2, quantile, 0.025) * 1000,
  est_upper_95 = apply(posterior$foi, 2, quantile, 0.975) * 1000,
  est_lower_50 = apply(posterior$foi, 2, quantile, 0.25) * 1000,
  est_upper_50 = apply(posterior$foi, 2, quantile, 0.75) * 1000
)

# Check coverage: does 95% CrI contain true value?
foi_comparison <- foi_comparison %>%
  mutate(
    covered_95 = true_foi >= est_lower_95 & true_foi <= est_upper_95,
    covered_50 = true_foi >= est_lower_50 & true_foi <= est_upper_50,
    bias = est_mean - true_foi,
    rel_bias = (est_mean - true_foi) / true_foi * 100
  )

# Summary statistics
cat("\n=== RECOVERY STATISTICS ===\n")
cat(sprintf(
  "95%% CrI coverage: %.1f%% (expected ~95%%)\n",
  mean(foi_comparison$covered_95) * 100
))
cat(sprintf(
  "50%% CrI coverage: %.1f%% (expected ~50%%)\n",
  mean(foi_comparison$covered_50) * 100
))
cat(sprintf("Mean bias: %.3f per 1000\n", mean(foi_comparison$bias)))
cat(sprintf(
  "Mean absolute bias: %.3f per 1000\n",
  mean(abs(foi_comparison$bias))
))
cat(sprintf("RMSE: %.3f per 1000\n", sqrt(mean(foi_comparison$bias^2))))

# Correlation between true and estimated
cor_foi <- cor(foi_comparison$true_foi, foi_comparison$est_mean)
cat(sprintf("Correlation (true vs estimated): %.4f\n", cor_foi))

# =============================================================================
# SECTION 10: VISUALIZATION
# =============================================================================

cat("\n--- Generating comparison plots ---\n")

# --- Plot 1: FOI trajectory comparison ---
p1 <- ggplot(foi_comparison, aes(x = year)) +
  # 95% CrI
  geom_ribbon(
    aes(ymin = est_lower_95, ymax = est_upper_95),
    fill = "steelblue",
    alpha = 0.2
  ) +
  # 50% CrI
  geom_ribbon(
    aes(ymin = est_lower_50, ymax = est_upper_50),
    fill = "steelblue",
    alpha = 0.3
  ) +
  # Estimated mean
  geom_line(aes(y = est_mean, color = "Estimated"), linewidth = 1) +
  # True FOI
  geom_line(
    aes(y = true_foi, color = "True"),
    linewidth = 1.2,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("Estimated" = "steelblue", "True" = "darkred"),
    name = ""
  ) +
  labs(
    title = "FOI Recovery: Estimated vs True",
    subtitle = sprintf(
      "%s scenario | 95%% coverage: %.0f%%",
      scenario,
      mean(foi_comparison$covered_95) * 100
    ),
    x = "Calendar Year",
    y = "FOI (per 1,000 person-years)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

# --- Plot 2: Scatter plot of true vs estimated ---
p2 <- ggplot(foi_comparison, aes(x = true_foi, y = est_mean)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
  geom_errorbar(
    aes(ymin = est_lower_95, ymax = est_upper_95),
    width = 0,
    alpha = 0.3,
    color = "steelblue"
  ) +
  geom_point(aes(color = year), size = 2) +
  scale_color_viridis_c(name = "Year") +
  labs(
    title = "True vs Estimated FOI",
    subtitle = sprintf("Correlation: %.3f", cor_foi),
    x = "True FOI (per 1,000)",
    y = "Estimated FOI (per 1,000)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# --- Plot 3: Bias over time ---
p3 <- ggplot(foi_comparison, aes(x = year, y = bias)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(aes(color = covered_95), size = 2) +
  scale_color_manual(
    values = c("TRUE" = "steelblue", "FALSE" = "red"),
    name = "95% CrI\ncovers true"
  ) +
  labs(
    title = "Estimation Bias Over Time",
    subtitle = "Bias = Estimated - True",
    x = "Calendar Year",
    y = "Bias (per 1,000)"
  ) +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold"))

# --- Plot 4: Posterior predictive check (faceted by year) ---
pred_samples <- posterior$n_pos_pred

ppc_sim <- sim_data %>%
  mutate(
    pred_mean = apply(pred_samples, 2, mean) / n_sample * 100,
    pred_lower = apply(pred_samples, 2, quantile, 0.025) / n_sample * 100,
    pred_upper = apply(pred_samples, 2, quantile, 0.975) / n_sample * 100,
    obs_prev_pct = obs_prev * 100,
    true_prev_pct = true_prev * 100
  )

p4 <- ggplot(ppc_sim %>% filter(survey_year %in% show_years), aes(x = age)) +
  geom_ribbon(
    aes(ymin = pred_lower, ymax = pred_upper),
    fill = "steelblue",
    alpha = 0.3
  ) +
  geom_line(aes(y = pred_mean, color = "Predicted"), linewidth = 0.8) +
  geom_line(
    aes(y = true_prev_pct, color = "True"),
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  geom_point(aes(y = obs_prev_pct, color = "Observed"), size = 1, alpha = 0.5) +
  facet_wrap(~survey_year, ncol = 3) +
  scale_color_manual(
    values = c(
      "Predicted" = "steelblue",
      "True" = "darkred",
      "Observed" = "black"
    ),
    name = ""
  ) +
  labs(
    title = "Posterior Predictive Check (Multi-Year)",
    subtitle = "Age-seroprevalence curves by survey year",
    x = "Age (years)",
    y = "Seroprevalence (%)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold")
  )

# Combine plots
combined_plot <- (p1 + p2) /
  (p3 + p4) +
  plot_annotation(
    title = sprintf(
      "Simulation Study Results - %s Scenario (Multi-Year)",
      str_to_title(scenario)
    ),
    subtitle = sprintf(
      "n = %s | %d observations | %d survey years | %d FOI years",
      format(sum(sim_data$n_sample), big.mark = ","),
      nrow(sim_data),
      length(survey_years),
      n_years
    ),
    theme = theme(
      plot.title = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

print(combined_plot)

# Save plot
ggsave(
  "foi_simulation_results_multiyear.png",
  combined_plot,
  width = 16,
  height = 14,
  dpi = 300
)
cat("Plot saved to: foi_simulation_results_multiyear.png\n")

# =============================================================================
# SECTION 11: PERIOD-SPECIFIC RECOVERY (for piecewise scenario)
# =============================================================================

if (scenario %in% c("piecewise", "epidemic")) {
  cat("\n--- Period-Specific Recovery ---\n")

  if (scenario == "piecewise") {
    period_labels <- case_when(
      foi_comparison$year <= 2000 ~ "1943-2000",
      foi_comparison$year <= 2010 ~ "2001-2010",
      foi_comparison$year <= 2015 ~ "2011-2015",
      TRUE ~ "2016-2022"
    )
  } else {
    # epidemic scenario: periods matching the double-wave pattern
    period_labels <- case_when(
      foi_comparison$year <= 1960 ~ "1943-1960 (baseline)",
      foi_comparison$year <= 1980 ~ "1961-1980 (1st wave)",
      foi_comparison$year <= 2000 ~ "1981-2000 (decline)",
      foi_comparison$year <= 2010 ~ "2001-2010 (trough/early rise)",
      TRUE ~ "2011-2022 (resurgence)"
    )
  }

  period_recovery <- foi_comparison %>%
    mutate(period = period_labels) %>%
    group_by(period) %>%
    summarise(
      true_mean = mean(true_foi),
      est_mean = mean(est_mean),
      bias = mean(bias),
      coverage_95 = mean(covered_95) * 100,
      .groups = "drop"
    )

  print(period_recovery)
}

# =============================================================================
# SECTION 12: SAVE RESULTS
# =============================================================================

cat("\n--- Saving results ---\n")

# Save comparison data
write.csv(foi_comparison, "foi_simulation_comparison_multiyear.csv", row.names = FALSE)

# Save summary
sink("foi_simulation_summary_multiyear.txt")
cat("=" |> rep(60) |> paste(collapse = ""), "\n")
cat("FOI SIMULATION STUDY - MULTI-YEAR SUMMARY\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Scenario: %s\n", scenario))
cat(sprintf("Survey years: %d-%d (%d years)\n",
            min(survey_years), max(survey_years), length(survey_years)))
cat(sprintf("Sample size per cell: %d\n", n_per_age))
cat(sprintf("Total observations: %d\n", nrow(sim_data)))
cat(sprintf("Total sample: %d\n", sum(sim_data$n_sample)))
cat(sprintf("Age range: %d-%d years\n", min_age, max_age))
cat(sprintf("FOI years estimated: %d-%d (%d)\n\n", min(years), max(years), n_years))

cat("RECOVERY METRICS:\n")
cat(sprintf(
  "  95%% CrI coverage: %.1f%%\n",
  mean(foi_comparison$covered_95) * 100
))
cat(sprintf(
  "  50%% CrI coverage: %.1f%%\n",
  mean(foi_comparison$covered_50) * 100
))
cat(sprintf("  Correlation: %.4f\n", cor_foi))
cat(sprintf("  Mean bias: %.4f per 1000\n", mean(foi_comparison$bias)))
cat(sprintf("  RMSE: %.4f per 1000\n", sqrt(mean(foi_comparison$bias^2))))

cat("\nMCMC DIAGNOSTICS:\n")
cat(sprintf("  Divergent transitions: %d\n", n_divergent))
cat(sprintf("  Max Rhat: %.4f\n", max_rhat))
sink()

cat("Results saved\n")

# =============================================================================
# SECTION 13: FINAL SUMMARY
# =============================================================================

cat("\n", "=" |> rep(60) |> paste(collapse = ""), "\n")
cat("SIMULATION STUDY COMPLETE\n")
cat("=" |> rep(60) |> paste(collapse = ""), "\n\n")

cat("Key findings:\n")
if (mean(foi_comparison$covered_95) >= 0.90) {
  cat("  Good 95% coverage - model is well-calibrated\n")
} else {
  cat("  Low 95% coverage - model may be overconfident\n")
}

if (cor_foi >= 0.90) {
  cat("  High correlation - model captures FOI trends well\n")
} else if (cor_foi >= 0.70) {
  cat("  Moderate correlation - model captures general trends\n")
} else {
  cat("  Low correlation - model struggles to recover FOI\n")
}

if (abs(mean(foi_comparison$bias)) < 0.5) {
  cat("  Low bias - estimates are accurate on average\n")
} else {
  cat("  Notable bias - systematic over/underestimation\n")
}

cat("\nOutput files:\n")
cat("  - foi_simulation_results_multiyear.png (visual comparison)\n")
cat("  - foi_simulation_comparison_multiyear.csv (detailed data)\n")
cat("  - foi_simulation_summary_multiyear.txt (summary statistics)\n")

cat("\nALL COMPLETE\n")
