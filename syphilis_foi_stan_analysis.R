# =============================================================================
# SYPHILIS FORCE OF INFECTION ANALYSIS - CUSTOM STAN MODELS (MULTI-YEAR)
# =============================================================================
# Multi-year analysis (2008-2023) using 4 custom Stan models:
#   1. Constant FOI (baseline)
#   2. Piecewise-3 periods (REDS eras)
#   3. Piecewise-5 periods (finer resolution)
#   4. Smooth Random Walk (most flexible)
#
# Data: Raw CSV files from REDS-II, REDS-III, REDS-IV
# Aggregation: by (survey_year, age) for multi-year design
# =============================================================================

rm(list = ls())
gc()
set.seed(42)

options(bitmapType = "cairo")

# =============================================================================
# SECTION 1: SETUP AND CONFIGURATION
# =============================================================================

cat("=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SYPHILIS FOI ANALYSIS - CUSTOM STAN MODELS (MULTI-YEAR)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Required packages
required_packages <- c(
  "tidyverse", "data.table", "rstan", "loo", "bayesplot",
  "ggpubr", "patchwork", "scales", "viridis", "posterior"
)

install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cloud.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
invisible(sapply(required_packages, install_if_missing))

# HPC configuration
.libPaths(c("/scratch/7631403/R_packages", .libPaths()))

n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
options(mc.cores = n_cores)
rstan_options(auto_write = TRUE)

cat(sprintf("Using %d cores\n", n_cores))
cat(sprintf("Stan version: %s\n", stan_version()))

# Paths - dual support for HPC and local
if (dir.exists("/scratch/7631403/reds_data/")) {
  base_path <- "/scratch/7631403/reds_data/"
  output_path <- "/scratch/7631403/reds_data/outputs/foi_stan/"
  stan_path <- "/scratch/7631403/reds_data/stan_models/"
} else {
  base_path <- "c:/Users/paulo/OneDrive/Documentos/AnalistaZen/REDS/"
  output_path <- paste0(base_path, "output_stan/")
  stan_path <- paste0(base_path, "transfer/new/stan_models/")
}

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("Base path: %s\n", base_path))
cat(sprintf("Output path: %s\n", output_path))
cat(sprintf("Stan models path: %s\n", stan_path))

# =============================================================================
# SECTION 2: LOAD RAW DATA FROM CSV FILES
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("LOADING RAW DATA FROM CSV FILES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

YEAR_MIN <- 2008
YEAR_MAX <- 2023

datasets_config <- list(
  list(study = "REDS-II",
       donation = paste0(base_path, "reds_II_donation.csv"),
       short = paste0(base_path, "reds_II_short.csv")),
  list(study = "REDS-III",
       donation = paste0(base_path, "reds_III_donation.csv"),
       short = paste0(base_path, "reds_III_short.csv")),
  list(study = "REDS-IV",
       donation = paste0(base_path, "reds_IV_donation_07212025_FINAL.csv"),
       short = paste0(base_path, "reds_IV_short_07212025_FINAL.csv"))
)

cat("Checking input files...\n")
for (config in datasets_config) {
  for (f in c(config$donation, config$short)) {
    if (file.exists(f)) {
      cat(sprintf("  OK %s\n", basename(f)))
    } else {
      cat(sprintf("  MISSING: %s\n", f))
    }
  }
}

# --- Variable mappings ---
center_mapping <- c(
  "31" = "Hemope, Recife", "32" = "Hemominas, Belo Horizonte",
  "33" = "FPS, Sao Paulo", "35" = "Hemoam, Manaus",
  "37" = "Hemorio, Rio de Janeiro"
)

center_state <- c(
  "Hemope, Recife" = "PE", "Hemominas, Belo Horizonte" = "MG",
  "FPS, Sao Paulo" = "SP", "Hemoam, Manaus" = "AM",
  "Hemorio, Rio de Janeiro" = "RJ"
)

sex_mapping <- c("M" = "Male", "F" = "Female", "8" = "Unknown", "9" = "Unknown")

# --- Data loading function ---
load_study_data <- function(config) {
  cat(sprintf("  Loading %s...\n", config$study))

  don <- fread(config$donation, na.strings = c("", "NA", "-", "--", "---"))
  setnames(don, tolower(names(don)))

  donation_cols <- c("centerid", "bloodid", "donorid", "donyr", "donmo", "donda",
                     "birthyr", "sex", "prevscrnhx", "dontype", "syp")
  donation_cols_available <- intersect(donation_cols, names(don))
  don <- don[, ..donation_cols_available]

  cat(sprintf("    Donations: %s records\n", format(nrow(don), big.mark = ",")))

  short <- fread(config$short, na.strings = c("", "NA", "-", "--", "---"))
  setnames(short, tolower(names(short)))

  short_cols <- c("bloodid", "visityr", "visitmo", "visitda", "declared_race", "educatn")
  short_cols_available <- intersect(short_cols, names(short))
  short <- short[, ..short_cols_available]

  cat(sprintf("    Short: %s records\n", format(nrow(short), big.mark = ",")))

  merged <- merge(don, short, by = "bloodid", all.x = TRUE)
  merged[, study := config$study]

  cat(sprintf("    Merged: %s records\n", format(nrow(merged), big.mark = ",")))

  rm(don, short); gc()
  return(merged)
}

# --- Data processing function ---
process_dataset <- function(dt) {
  dt[, donyr := as.integer(donyr)]
  dt[donyr %in% c(9998, 9999), donyr := NA_integer_]

  dt[, visityr := as.integer(visityr)]
  dt[visityr %in% c(9998, 9999), visityr := NA_integer_]

  dt[, sample_year := fifelse(!is.na(visityr), visityr, donyr)]

  dt <- dt[sample_year >= YEAR_MIN & sample_year <= YEAR_MAX]
  cat(sprintf("  After year filter (%d-%d): %s records\n",
              YEAR_MIN, YEAR_MAX, format(nrow(dt), big.mark = ",")))

  dt[, birthyr := as.integer(birthyr)]
  dt[birthyr %in% c(9998, 9999), birthyr := NA_integer_]
  dt[, age := sample_year - birthyr]
  dt[age < 16 | age > 70, age := NA_integer_]

  dt[, age_group := cut(age,
                        breaks = c(18, 25, 35, 45, 55, 66),
                        labels = c("18-24", "25-34", "35-44", "45-54", "55-65"),
                        right = FALSE)]

  dt[, centerid := as.character(centerid)]
  dt[, center_name := center_mapping[centerid]]
  dt[, state := center_state[center_name]]

  dt[, sex := as.character(sex)]
  dt[, sex_mapped := sex_mapping[sex]]
  dt[is.na(sex_mapped), sex_mapped := "Unknown"]

  dt[, syp := toupper(as.character(syp))]
  dt <- dt[syp %in% c("P", "N")]
  cat(sprintf("  After syphilis filter (P/N only): %s records\n",
              format(nrow(dt), big.mark = ",")))

  # Center-specific exclusions
  dt <- dt[!(centerid == "31" & sample_year < 2009)]
  dt <- dt[!(centerid == "37" & sample_year < 2013)]
  cat(sprintf("  After center-specific exclusions: %s records\n",
              format(nrow(dt), big.mark = ",")))

  dt[, syphilis_positive := fifelse(syp == "P", 1L, 0L)]

  # First donation per donor
  dt[, donmo := as.integer(donmo)]
  dt[donmo %in% c(98, 99), donmo := 6L]
  dt[, donda := as.integer(donda)]
  dt[donda %in% c(98, 99), donda := 15L]

  setorder(dt, donorid, sample_year, donmo, donda)
  dt <- dt[, .SD[1], by = donorid]
  cat(sprintf("  First donation per donor: %s unique donors\n",
              format(nrow(dt), big.mark = ",")))

  dt <- dt[!is.na(age) & !is.na(sample_year) & !is.na(center_name)]
  cat(sprintf("  Final clean dataset: %s donors\n",
              format(nrow(dt), big.mark = ",")))

  return(dt)
}

# --- Load and process all data ---
cat("\nLoading all datasets...\n")
all_data <- lapply(datasets_config, load_study_data)
combined_raw <- rbindlist(all_data, fill = TRUE)
rm(all_data); gc()

cat(sprintf("\nCombined raw: %s records\n", format(nrow(combined_raw), big.mark = ",")))

df_clean <- as.data.frame(process_dataset(combined_raw))
rm(combined_raw); gc()

cat(sprintf("\nTotal unique donors: %s\n", format(nrow(df_clean), big.mark = ",")))
cat(sprintf("Syphilis positive: %s (%.2f%%)\n",
            format(sum(df_clean$syphilis_positive), big.mark = ","),
            mean(df_clean$syphilis_positive) * 100))
cat(sprintf("Year range: %d - %d\n", min(df_clean$sample_year), max(df_clean$sample_year)))
cat(sprintf("Age range: %d - %d years\n", min(df_clean$age), max(df_clean$age)))

# =============================================================================
# SECTION 3: PREPARE MULTI-YEAR STAN DATA
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("PREPARING MULTI-YEAR STAN DATA\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Aggregate by (sample_year, age) - multi-year design
serosurvey_multiyear <- df_clean %>%
  filter(age >= 18 & age <= 65) %>%
  group_by(sample_year, age) %>%
  summarise(
    n_sample = n(),
    n_seropositive = sum(syphilis_positive),
    .groups = "drop"
  ) %>%
  filter(n_sample >= 10) %>%
  arrange(sample_year, age)

cat(sprintf("Multi-year aggregated data: %d rows (year x age combinations)\n",
            nrow(serosurvey_multiyear)))
cat(sprintf("Survey years: %d to %d\n",
            min(serosurvey_multiyear$sample_year),
            max(serosurvey_multiyear$sample_year)))
cat(sprintf("Total N: %s\n", format(sum(serosurvey_multiyear$n_sample), big.mark = ",")))
cat(sprintf("Total positive: %s\n", format(sum(serosurvey_multiyear$n_seropositive), big.mark = ",")))

# Compute FOI time frame from data
min_birth_year <- min(serosurvey_multiyear$sample_year - serosurvey_multiyear$age)
max_survey_year <- max(serosurvey_multiyear$sample_year)
n_years <- (max_survey_year - 1) - min_birth_year + 1
years <- seq(min_birth_year, max_survey_year - 1)

cat(sprintf("\nFOI estimation frame:\n"))
cat(sprintf("  Min birth year: %d\n", min_birth_year))
cat(sprintf("  Max survey year: %d\n", max_survey_year))
cat(sprintf("  FOI years: %d to %d (%d years)\n", min(years), max(years), n_years))

# Base Stan data (common to all models)
stan_data_base <- list(
  N = nrow(serosurvey_multiyear),
  n_sample = serosurvey_multiyear$n_sample,
  n_pos = serosurvey_multiyear$n_seropositive,
  age = serosurvey_multiyear$age,
  survey_year = serosurvey_multiyear$sample_year,  # VECTOR: per-observation
  n_years = n_years,
  min_year = min_birth_year
)

cat("\nStan data structure:\n")
cat(sprintf("  N (observations): %d\n", stan_data_base$N))
cat(sprintf("  n_years: %d\n", stan_data_base$n_years))
cat(sprintf("  min_year: %d\n", stan_data_base$min_year))
cat(sprintf("  survey_year: vector of length %d (range %d-%d)\n",
            length(stan_data_base$survey_year),
            min(stan_data_base$survey_year), max(stan_data_base$survey_year)))

# Save aggregated data
write.csv(serosurvey_multiyear, paste0(output_path, "serosurvey_multiyear.csv"),
          row.names = FALSE)

# =============================================================================
# SECTION 4: DEFINE FOI PERIOD INDICES
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("DEFINING FOI PERIOD STRUCTURES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# --- Model 2: Piecewise-3 (REDS program eras) ---
foi_index_3 <- case_when(
  years <= 2011 ~ 1L,
  years <= 2016 ~ 2L,
  TRUE ~ 3L
)

cat("Piecewise-3 periods:\n")
cat(sprintf("  Period 1 (%d-%d): %d years\n", min(years), 2011, sum(foi_index_3 == 1)))
cat(sprintf("  Period 2 (2012-2016): %d years\n", sum(foi_index_3 == 2)))
cat(sprintf("  Period 3 (2017-%d): %d years\n", max(years), sum(foi_index_3 == 3)))

# --- Model 3: Piecewise-5 (finer resolution) ---
foi_index_5 <- case_when(
  years <= 2010 ~ 1L,
  years <= 2013 ~ 2L,
  years <= 2017 ~ 3L,
  years <= 2020 ~ 4L,
  TRUE ~ 5L
)

cat("\nPiecewise-5 periods:\n")
for (k in 1:5) {
  idx <- which(foi_index_5 == k)
  cat(sprintf("  Period %d (%d-%d): %d years\n",
              k, years[min(idx)], years[max(idx)], length(idx)))
}

# Store period labels for later use
period_labels_3 <- c(
  sprintf("%d-2011 (Historical + REDS-II)", min(years)),
  "2012-2016 (REDS-III)",
  sprintf("2017-%d (REDS-IV)", max(years))
)

period_labels_5 <- c(
  sprintf("%d-2010 (Baseline)", min(years)),
  "2011-2013",
  "2014-2017",
  "2018-2020",
  sprintf("2021-%d", max(years))
)

# =============================================================================
# SECTION 5: COMPILE STAN MODELS
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("COMPILING STAN MODELS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

stan_files <- c(
  constant = "foi_constant.stan",
  piecewise = "foi_piecewise.stan",
  random_walk = "foi_random_walk.stan",
  random_walk_v2 = "foi_random_walk_v2.stan"
)

models <- list()

for (model_name in names(stan_files)) {
  stan_file <- paste0(stan_path, stan_files[model_name])

  if (!file.exists(stan_file)) {
    cat(sprintf("Stan file not found: %s\n", stan_file))
    next
  }

  cat(sprintf("Compiling %s model...\n", model_name))
  tryCatch({
    models[[model_name]] <- stan_model(stan_file)
    cat(sprintf("  %s compiled successfully\n", model_name))
  }, error = function(e) {
    cat(sprintf("  ERROR compiling %s: %s\n", model_name, e$message))
  })
}

cat(sprintf("\n%d models compiled successfully\n", length(models)))

# =============================================================================
# SECTION 6: FIT MODEL 1 - CONSTANT FOI
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL 1: CONSTANT FOI\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if ("constant" %in% names(models)) {

  cat("Fitting constant FOI model...\n")
  t1 <- Sys.time()

  fit_constant <- sampling(
    models$constant,
    data = stan_data_base,
    iter = 4000,
    warmup = 1500,
    chains = 4,
    cores = n_cores,
    seed = 42,
    control = list(adapt_delta = 0.95)
  )

  t2 <- Sys.time()
  cat(sprintf("Model 1 fitted in %.1f minutes\n", difftime(t2, t1, units = "mins")))

  # Diagnostics
  cat("\n--- Diagnostics ---\n")
  print(summary(fit_constant, pars = c("foi", "foi_per_1000"))$summary)

  rhat_constant <- summary(fit_constant)$summary[, "Rhat"]
  cat(sprintf("Max Rhat: %.4f\n", max(rhat_constant, na.rm = TRUE)))

} else {
  cat("Constant model not available - skipping\n")
  fit_constant <- NULL
}

# =============================================================================
# SECTION 7: FIT MODEL 2 - PIECEWISE-3
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL 2: PIECEWISE-3 (REDS ERAS)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if ("piecewise" %in% names(models)) {

  stan_data_pw3 <- c(stan_data_base, list(
    K = 3,
    foi_index = foi_index_3
  ))

  cat("Fitting piecewise-3 FOI model...\n")
  t1 <- Sys.time()

  fit_piecewise_3 <- sampling(
    models$piecewise,
    data = stan_data_pw3,
    iter = 5000,
    warmup = 2000,
    chains = 4,
    cores = n_cores,
    seed = 42,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )

  t2 <- Sys.time()
  cat(sprintf("Model 2 fitted in %.1f minutes\n", difftime(t2, t1, units = "mins")))

  cat("\n--- Diagnostics ---\n")
  print(summary(fit_piecewise_3, pars = c("foi", "foi_per_1000"))$summary)

} else {
  cat("Piecewise model not available - skipping\n")
  fit_piecewise_3 <- NULL
}

# =============================================================================
# SECTION 8: FIT MODEL 3 - PIECEWISE-5
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL 3: PIECEWISE-5 (FINE RESOLUTION)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if ("piecewise" %in% names(models)) {

  stan_data_pw5 <- c(stan_data_base, list(
    K = 5,
    foi_index = foi_index_5
  ))

  cat("Fitting piecewise-5 FOI model...\n")
  t1 <- Sys.time()

  fit_piecewise_5 <- sampling(
    models$piecewise,
    data = stan_data_pw5,
    iter = 5000,
    warmup = 2000,
    chains = 4,
    cores = n_cores,
    seed = 42,
    control = list(adapt_delta = 0.95, max_treedepth = 12)
  )

  t2 <- Sys.time()
  cat(sprintf("Model 3 fitted in %.1f minutes\n", difftime(t2, t1, units = "mins")))

  cat("\n--- Diagnostics ---\n")
  print(summary(fit_piecewise_5, pars = c("foi", "foi_per_1000"))$summary)

} else {
  fit_piecewise_5 <- NULL
}

# =============================================================================
# SECTION 9: FIT MODEL 4 - SMOOTH RANDOM WALK
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL 4: SMOOTH RANDOM WALK\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

if ("random_walk_v2" %in% names(models)) {

  stan_data_rw <- c(stan_data_base, list(
    include_seroreversion = 0L,
    log_foi_mean_prior = -7.0,
    log_foi_sd_prior = 1.5,
    sigma_rw_upper = 0.5
  ))

  cat("Fitting random walk FOI model...\n")
  cat("  (This may take longer due to year-by-year estimation)\n")
  t1 <- Sys.time()

  fit_rw <- sampling(
    models$random_walk_v2,
    data = stan_data_rw,
    iter = 5000,
    warmup = 2500,
    chains = 4,
    cores = n_cores,
    seed = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )

  t2 <- Sys.time()
  cat(sprintf("Model 4 fitted in %.1f minutes\n", difftime(t2, t1, units = "mins")))

  cat("\n--- Diagnostics (summary) ---\n")
  rw_summary <- summary(fit_rw, pars = c("sigma_rw", "foi_recent"))$summary
  print(rw_summary)

} else if ("random_walk" %in% names(models)) {

  # Fallback to simpler RW model
  stan_data_rw <- c(stan_data_base, list(
    sigma_rw_prior_mean = 0.02,
    sigma_rw_prior_sd = 0.01
  ))

  cat("Fitting random walk FOI model (v1)...\n")
  t1 <- Sys.time()

  fit_rw <- sampling(
    models$random_walk,
    data = stan_data_rw,
    iter = 5000,
    warmup = 2500,
    chains = 4,
    cores = n_cores,
    seed = 42,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )

  t2 <- Sys.time()
  cat(sprintf("Model 4 fitted in %.1f minutes\n", difftime(t2, t1, units = "mins")))

} else {
  cat("Random walk model not available - skipping\n")
  fit_rw <- NULL
}

# =============================================================================
# SECTION 10: MODEL COMPARISON
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("MODEL COMPARISON (LOO-CV)\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

loo_results <- list()

model_fits <- list(
  "1_Constant" = fit_constant,
  "2_Piecewise3" = fit_piecewise_3,
  "3_Piecewise5" = fit_piecewise_5,
  "4_RandomWalk" = fit_rw
)

for (model_name in names(model_fits)) {
  fit <- model_fits[[model_name]]
  if (!is.null(fit)) {
    cat(sprintf("Computing LOO for %s...\n", model_name))
    tryCatch({
      log_lik <- extract_log_lik(fit, parameter_name = "log_lik")
      loo_results[[model_name]] <- loo(log_lik, cores = n_cores)
      cat(sprintf("  %s: ELPD = %.1f (SE = %.1f)\n",
                  model_name,
                  loo_results[[model_name]]$estimates["elpd_loo", "Estimate"],
                  loo_results[[model_name]]$estimates["elpd_loo", "SE"]))
    }, error = function(e) {
      cat(sprintf("  Error for %s: %s\n", model_name, e$message))
    })
  }
}

if (length(loo_results) >= 2) {
  cat("\n--- LOO Comparison ---\n")
  loo_comp <- loo_compare(loo_results)
  print(loo_comp)

  loo_df <- as.data.frame(loo_comp)
  loo_df$model <- rownames(loo_df)
  write.csv(loo_df, paste0(output_path, "model_comparison_loo.csv"), row.names = FALSE)
  cat("\nLOO comparison saved\n")
}

# =============================================================================
# SECTION 11: EXTRACT FOI ESTIMATES
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("EXTRACTING FOI ESTIMATES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Function to extract FOI posterior summary
extract_foi_summary <- function(fit, years, foi_param = "foi", per_1000 = TRUE) {
  posterior <- rstan::extract(fit)
  foi_samples <- posterior[[foi_param]]
  n_samples <- nrow(foi_samples)
  n_params <- ncol(foi_samples)

  if (n_params == length(years)) {
    result <- data.frame(
      year = years,
      mean = apply(foi_samples, 2, mean),
      median = apply(foi_samples, 2, median),
      lower_95 = apply(foi_samples, 2, quantile, 0.025),
      upper_95 = apply(foi_samples, 2, quantile, 0.975),
      lower_50 = apply(foi_samples, 2, quantile, 0.25),
      upper_50 = apply(foi_samples, 2, quantile, 0.75),
      sd = apply(foi_samples, 2, sd)
    )
  } else {
    result <- data.frame(
      period = 1:n_params,
      mean = apply(foi_samples, 2, mean),
      median = apply(foi_samples, 2, median),
      lower_95 = apply(foi_samples, 2, quantile, 0.025),
      upper_95 = apply(foi_samples, 2, quantile, 0.975),
      lower_50 = apply(foi_samples, 2, quantile, 0.25),
      upper_50 = apply(foi_samples, 2, quantile, 0.75),
      sd = apply(foi_samples, 2, sd)
    )
  }

  if (per_1000) {
    for (col in c("mean", "median", "lower_95", "upper_95", "lower_50", "upper_50", "sd")) {
      result[[col]] <- result[[col]] * 1000
    }
  }

  return(result)
}

foi_estimates <- list()

if (!is.null(fit_constant)) {
  cat("Extracting constant model estimates...\n")
  posterior_const <- rstan::extract(fit_constant)
  foi_estimates$constant <- data.frame(
    model = "Constant",
    foi = mean(posterior_const$foi),
    foi_lower = quantile(posterior_const$foi, 0.025),
    foi_upper = quantile(posterior_const$foi, 0.975),
    foi_per_1000 = mean(posterior_const$foi) * 1000,
    foi_per_1000_lower = quantile(posterior_const$foi, 0.025) * 1000,
    foi_per_1000_upper = quantile(posterior_const$foi, 0.975) * 1000
  )
  cat(sprintf("  Constant FOI: %.2f (%.2f-%.2f) per 1000\n",
      foi_estimates$constant$foi_per_1000,
      foi_estimates$constant$foi_per_1000_lower,
      foi_estimates$constant$foi_per_1000_upper))
}

if (!is.null(fit_piecewise_3)) {
  cat("Extracting piecewise-3 model estimates...\n")
  posterior_pw3 <- rstan::extract(fit_piecewise_3)

  foi_estimates$piecewise_3 <- data.frame(
    model = "Piecewise-3",
    period = 1:3,
    period_label = period_labels_3,
    foi_mean = apply(posterior_pw3$foi, 2, mean),
    foi_median = apply(posterior_pw3$foi, 2, median),
    foi_lower = apply(posterior_pw3$foi, 2, quantile, 0.025),
    foi_upper = apply(posterior_pw3$foi, 2, quantile, 0.975)
  ) %>%
    mutate(
      foi_per_1000 = foi_mean * 1000,
      foi_per_1000_lower = foi_lower * 1000,
      foi_per_1000_upper = foi_upper * 1000
    )

  print(foi_estimates$piecewise_3[, c("period_label", "foi_per_1000", "foi_per_1000_lower", "foi_per_1000_upper")])
}

if (!is.null(fit_piecewise_5)) {
  cat("\nExtracting piecewise-5 model estimates...\n")
  posterior_pw5 <- rstan::extract(fit_piecewise_5)

  foi_estimates$piecewise_5 <- data.frame(
    model = "Piecewise-5",
    period = 1:5,
    period_label = period_labels_5,
    foi_mean = apply(posterior_pw5$foi, 2, mean),
    foi_median = apply(posterior_pw5$foi, 2, median),
    foi_lower = apply(posterior_pw5$foi, 2, quantile, 0.025),
    foi_upper = apply(posterior_pw5$foi, 2, quantile, 0.975)
  ) %>%
    mutate(
      foi_per_1000 = foi_mean * 1000,
      foi_per_1000_lower = foi_lower * 1000,
      foi_per_1000_upper = foi_upper * 1000
    )

  print(foi_estimates$piecewise_5[, c("period_label", "foi_per_1000", "foi_per_1000_lower", "foi_per_1000_upper")])
}

if (!is.null(fit_rw)) {
  cat("\nExtracting random walk model estimates...\n")
  posterior_rw <- rstan::extract(fit_rw)

  foi_estimates$random_walk <- data.frame(
    model = "RandomWalk",
    year = years,
    foi_mean = apply(posterior_rw$foi, 2, mean),
    foi_median = apply(posterior_rw$foi, 2, median),
    foi_lower = apply(posterior_rw$foi, 2, quantile, 0.025),
    foi_upper = apply(posterior_rw$foi, 2, quantile, 0.975),
    foi_lower_50 = apply(posterior_rw$foi, 2, quantile, 0.25),
    foi_upper_50 = apply(posterior_rw$foi, 2, quantile, 0.75)
  ) %>%
    mutate(
      foi_per_1000 = foi_mean * 1000,
      foi_per_1000_lower = foi_lower * 1000,
      foi_per_1000_upper = foi_upper * 1000
    )

  cat("\nRecent FOI (2017-2022):\n")
  print(foi_estimates$random_walk %>%
          filter(year >= 2017) %>%
          select(year, foi_per_1000, foi_per_1000_lower, foi_per_1000_upper))
}

for (nm in names(foi_estimates)) {
  write.csv(foi_estimates[[nm]],
            paste0(output_path, sprintf("foi_estimates_%s.csv", nm)),
            row.names = FALSE)
}
cat("\nFOI estimates saved\n")

# =============================================================================
# SECTION 12: PUBLICATION FIGURES
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("GENERATING PUBLICATION FIGURES\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

# Lancet-style theme
theme_lancet <- theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "gray80")
  )

model_colors <- c(
  "Constant" = "#999999",
  "Piecewise-3" = "#E69F00",
  "Piecewise-5" = "#56B4E9",
  "RandomWalk" = "#009E73"
)

# --- Figure 1: FOI Model Comparison ---
if (!is.null(fit_rw) && !is.null(fit_piecewise_3)) {

  pw3_plot <- data.frame(
    year = c(min(years), 2011, 2012, 2016, 2017, max(years)),
    foi_per_1000 = c(rep(foi_estimates$piecewise_3$foi_per_1000[1], 2),
                     rep(foi_estimates$piecewise_3$foi_per_1000[2], 2),
                     rep(foi_estimates$piecewise_3$foi_per_1000[3], 2)),
    model = "Piecewise-3"
  )

  fig1 <- ggplot() +
    geom_ribbon(data = foi_estimates$random_walk,
                aes(x = year, ymin = foi_per_1000_lower, ymax = foi_per_1000_upper),
                fill = model_colors["RandomWalk"], alpha = 0.2) +
    geom_line(data = foi_estimates$random_walk,
              aes(x = year, y = foi_per_1000, color = "RandomWalk"),
              linewidth = 1) +
    geom_step(data = pw3_plot,
              aes(x = year, y = foi_per_1000, color = "Piecewise-3"),
              linewidth = 1.2) +
    geom_hline(yintercept = foi_estimates$constant$foi_per_1000,
               color = model_colors["Constant"], linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = c(2012, 2017), linetype = "dotted", color = "gray60") +
    scale_color_manual(values = model_colors, name = "Model") +
    scale_x_continuous(breaks = seq(round(min(years) / 10) * 10, max(years), 10)) +
    scale_y_continuous(limits = c(0, NA)) +
    labs(
      title = "Force of Infection Estimates - Syphilis in Brazilian Blood Donors",
      subtitle = sprintf("Comparison of model specifications (%d-%d)", min(years), max(years)),
      x = "Calendar Year",
      y = "Annual FOI (per 1,000 person-years)",
      caption = "Shaded region: 95% credible interval for Random Walk model\nVertical lines: REDS study period boundaries"
    ) +
    theme_lancet

  ggsave(paste0(output_path, "fig1_foi_model_comparison.png"),
         fig1, width = 12, height = 7, dpi = 300)
  ggsave(paste0(output_path, "fig1_foi_model_comparison.pdf"),
         fig1, width = 12, height = 7)
  cat("Figure 1 saved\n")
}

# --- Figure 2: Posterior Predictive Check (faceted by survey year subset) ---
if (!is.null(fit_rw)) {

  posterior_rw_pred <- rstan::extract(fit_rw)
  pred_samples <- posterior_rw_pred$n_pos_pred

  # Add predictions to multi-year data
  ppc_data <- serosurvey_multiyear %>%
    mutate(
      obs_prev = n_seropositive / n_sample,
      pred_mean = apply(pred_samples, 2, mean) / n_sample,
      pred_lower = apply(pred_samples, 2, quantile, 0.025) / n_sample,
      pred_upper = apply(pred_samples, 2, quantile, 0.975) / n_sample
    )

  # Show a subset of years for readability
  show_years <- c(2008, 2011, 2014, 2017, 2020, 2023)
  show_years <- show_years[show_years %in% unique(ppc_data$sample_year)]

  fig2 <- ggplot(ppc_data %>% filter(sample_year %in% show_years)) +
    geom_ribbon(aes(x = age, ymin = pred_lower * 100, ymax = pred_upper * 100),
                fill = "#0072B2", alpha = 0.3) +
    geom_line(aes(x = age, y = pred_mean * 100),
              color = "#0072B2", linewidth = 0.8) +
    geom_point(aes(x = age, y = obs_prev * 100, size = n_sample),
               alpha = 0.5, color = "#D55E00") +
    facet_wrap(~sample_year, ncol = 3) +
    scale_size_continuous(range = c(0.5, 3), name = "Sample Size") +
    labs(
      title = "Posterior Predictive Check - Seroprevalence by Age and Year",
      subtitle = "Observed (orange) vs. Predicted (blue, 95% CrI) - Random Walk model",
      x = "Age (years)",
      y = "Seroprevalence (%)"
    ) +
    theme_lancet +
    theme(strip.text = element_text(face = "bold"))

  ggsave(paste0(output_path, "fig2_posterior_predictive_multiyear.png"),
         fig2, width = 14, height = 10, dpi = 300)
  ggsave(paste0(output_path, "fig2_posterior_predictive_multiyear.pdf"),
         fig2, width = 14, height = 10)
  cat("Figure 2 saved\n")
}

# --- Figure 3: Period-Specific FOI (Piecewise-5) ---
if (!is.null(fit_piecewise_5)) {

  fig3_data <- foi_estimates$piecewise_5 %>%
    filter(period >= 2) %>%
    mutate(period_label = factor(period_label, levels = period_labels_5[-1]))

  fig3 <- ggplot(fig3_data, aes(x = period_label)) +
    geom_errorbar(aes(ymin = foi_per_1000_lower, ymax = foi_per_1000_upper),
                  width = 0.3, linewidth = 0.8, color = "#0072B2") +
    geom_point(aes(y = foi_per_1000), size = 4, color = "#0072B2") +
    labs(
      title = "Period-Specific Force of Infection (Piecewise-5 Model)",
      subtitle = "Brazilian Blood Donors",
      x = "Period",
      y = "Annual FOI (per 1,000 person-years)",
      caption = "Error bars: 95% credible intervals"
    ) +
    theme_lancet +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(paste0(output_path, "fig3_foi_by_period.png"),
         fig3, width = 9, height = 6, dpi = 300)
  ggsave(paste0(output_path, "fig3_foi_by_period.pdf"),
         fig3, width = 9, height = 6)
  cat("Figure 3 saved\n")
}

# --- Figure 4: MCMC Diagnostics ---
if (!is.null(fit_rw)) {

  png(paste0(output_path, "fig4_mcmc_diagnostics.png"),
      width = 12, height = 8, units = "in", res = 300)

  p1 <- mcmc_trace(fit_rw, pars = c("sigma_rw", "log_foi_mean"))
  p2 <- mcmc_dens_overlay(fit_rw, pars = c("sigma_rw", "log_foi_mean"))

  print(p1 / p2)
  dev.off()
  cat("Figure 4 saved\n")
}

# =============================================================================
# SECTION 13: SUPPLEMENTARY - AGE-PREVALENCE SLOPE BY SURVEY YEAR
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SUPPLEMENTARY: AGE-PREVALENCE SLOPE ANALYSIS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

slope_data <- serosurvey_multiyear %>%
  mutate(prevalence = n_seropositive / n_sample)

# Slope by survey year
slope_by_year <- slope_data %>%
  group_by(sample_year) %>%
  summarise(
    n_ages = n(),
    n_total = sum(n_sample),
    slope = coef(lm(prevalence ~ age, weights = n_sample))["age"] * 1000,
    .groups = "drop"
  )

cat("Slope by survey year (per 1000 per year of age):\n")
print(slope_by_year)

write.csv(slope_by_year, paste0(output_path, "age_prevalence_slopes_by_year.csv"),
          row.names = FALSE)

# Slope by birth decade (aggregated across all survey years)
slope_by_decade <- slope_data %>%
  mutate(birth_decade = floor((sample_year - age) / 10) * 10) %>%
  group_by(birth_decade) %>%
  summarise(
    n_ages = n(),
    n_total = sum(n_sample),
    slope = coef(lm(prevalence ~ age, weights = n_sample))["age"] * 1000,
    .groups = "drop"
  )

cat("\nSlope by birth decade:\n")
print(slope_by_decade)

write.csv(slope_by_decade, paste0(output_path, "age_prevalence_slopes_by_decade.csv"),
          row.names = FALSE)

# =============================================================================
# SECTION 14: SAVE WORKSPACE AND SESSION INFO
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("SAVING RESULTS\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

save(
  serosurvey_multiyear, df_clean,
  stan_data_base,
  fit_constant, fit_piecewise_3, fit_piecewise_5, fit_rw,
  foi_estimates, loo_results,
  years, foi_index_3, foi_index_5,
  period_labels_3, period_labels_5,
  min_birth_year, max_survey_year, n_years,
  file = paste0(output_path, "syphilis_foi_stan_analysis.RData")
)

cat("Workspace saved to:", paste0(output_path, "syphilis_foi_stan_analysis.RData\n"))

# Session info
sink(paste0(output_path, "session_info_stan.txt"))
cat("SESSION INFO - STAN FOI ANALYSIS (MULTI-YEAR)\n")
cat(sprintf("Analysis completed: %s\n\n", Sys.time()))
cat(sprintf("Cores used: %d\n\n", n_cores))
cat(sprintf("N observations: %d\n", nrow(serosurvey_multiyear)))
cat(sprintf("N years FOI: %d (%d-%d)\n", n_years, min(years), max(years)))
sessionInfo()
sink()

# =============================================================================
# SECTION 15: FINAL SUMMARY
# =============================================================================

cat("\n", "=" |> rep(70) |> paste(collapse = ""), "\n")
cat("ANALYSIS COMPLETE\n")
cat("=" |> rep(70) |> paste(collapse = ""), "\n\n")

cat("Data structure:\n")
cat(sprintf("  Multi-year observations: %d (year x age)\n", nrow(serosurvey_multiyear)))
cat(sprintf("  Survey years: %d-%d\n", min(serosurvey_multiyear$sample_year),
            max(serosurvey_multiyear$sample_year)))
cat(sprintf("  FOI years estimated: %d-%d (%d years)\n", min(years), max(years), n_years))

cat("\nModels fitted:\n")
cat(sprintf("  1. Constant FOI: %s\n", ifelse(!is.null(fit_constant), "OK", "SKIP")))
cat(sprintf("  2. Piecewise-3:  %s\n", ifelse(!is.null(fit_piecewise_3), "OK", "SKIP")))
cat(sprintf("  3. Piecewise-5:  %s\n", ifelse(!is.null(fit_piecewise_5), "OK", "SKIP")))
cat(sprintf("  4. Random Walk:  %s\n", ifelse(!is.null(fit_rw), "OK", "SKIP")))

cat("\nKey outputs:\n")
cat("  - serosurvey_multiyear.csv: Multi-year aggregated data\n")
cat("  - model_comparison_loo.csv: LOO-CV comparison\n")
cat("  - foi_estimates_*.csv: FOI estimates by model\n")
cat("  - fig1-4_*.png/pdf: Publication figures\n")
cat("  - syphilis_foi_stan_analysis.RData: Full workspace\n")

cat("\nALL COMPLETE\n")
