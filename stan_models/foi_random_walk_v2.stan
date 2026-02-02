// =============================================================================
// SMOOTH RANDOM WALK FOI MODEL - OPTIMIZED PARAMETERIZATION
// =============================================================================
// Uses non-centered parameterization for better MCMC sampling
// Includes optional seroreversion (set rho = 0 to disable)
// =============================================================================

data {
  int<lower=1> N;                      // Number of age groups
  array[N] int<lower=1> n_sample;      // Sample size per age group
  array[N] int<lower=0> n_pos;         // Number seropositive per age group
  array[N] int<lower=1> age;           // Age (years) for each group
  array[N] int<lower=1> survey_year;   // Year of survey for each observation
  int<lower=1> n_years;                // Number of years for FOI estimation
  int<lower=1> min_year;               // First calendar year in FOI series
  
  // Optional seroreversion
  int<lower=0, upper=1> include_seroreversion;
  
  // Prior hyperparameters
  real log_foi_mean_prior;             // Prior mean for average log-FOI (e.g., -5)
  real<lower=0> log_foi_sd_prior;      // Prior SD for average log-FOI (e.g., 1.5)
  real<lower=0> sigma_rw_upper;        // Upper bound for sigma_rw (e.g., 0.5)
}

parameters {
  real log_foi_mean;                   // Mean log-FOI level
  vector[n_years] log_foi_raw;         // Non-centered RW innovations
  real<lower=0, upper=sigma_rw_upper> sigma_rw;  // RW standard deviation
  
  // Seroreversion rate (only used if include_seroreversion = 1)
  real<lower=0, upper=0.1> rho_raw;    
}

transformed parameters {
  vector[n_years] log_foi;
  vector<lower=0>[n_years] foi;
  real<lower=0> rho;
  
  // Non-centered random walk construction
  log_foi[1] = log_foi_mean + sigma_rw * log_foi_raw[1];
  for (t in 2:n_years) {
    log_foi[t] = log_foi[t-1] + sigma_rw * log_foi_raw[t];
  }
  
  // Transform to FOI scale with soft upper bound at 0.5
  for (t in 1:n_years) {
    foi[t] = 0.5 * inv_logit(log_foi[t]);  // Range: 0 to 0.5
  }
  
  // Seroreversion
  rho = include_seroreversion ? rho_raw : 0.0;
}

model {
  // Priors
  log_foi_mean ~ normal(log_foi_mean_prior, log_foi_sd_prior);
  log_foi_raw ~ std_normal();
  sigma_rw ~ exponential(20);  // Favors smooth trajectories
  
  if (include_seroreversion) {
    rho_raw ~ beta(1, 50);  // Prior mean ~0.02
  }
  
  // Likelihood
  for (i in 1:N) {
    int birth_year = survey_year[i] - age[i];
    real prob_pos = 0.0;  // Probability of being seropositive
    
    if (include_seroreversion) {
      // With seroreversion: P(t) = P(t-1)*(1-rho) + (1-P(t-1))*foi(t)
      for (y in 1:age[i]) {
        int cal_year = birth_year + y - 1;
        int year_idx = cal_year - min_year + 1;
        
        if (year_idx >= 1 && year_idx <= n_years) {
          prob_pos = prob_pos * (1.0 - rho) + (1.0 - prob_pos) * foi[year_idx];
        }
      }
    } else {
      // Without seroreversion: P = 1 - prod(1 - foi)
      real prob_neg = 1.0;
      for (y in 1:age[i]) {
        int cal_year = birth_year + y - 1;
        int year_idx = cal_year - min_year + 1;
        
        if (year_idx >= 1 && year_idx <= n_years) {
          prob_neg = prob_neg * (1.0 - foi[year_idx]);
        }
      }
      prob_pos = 1.0 - prob_neg;
    }
    
    // Bound probability to avoid numerical issues
    prob_pos = fmax(1e-10, fmin(1.0 - 1e-10, prob_pos));
    n_pos[i] ~ binomial(n_sample[i], prob_pos);
  }
}

generated quantities {
  array[N] int n_pos_pred;
  vector[N] log_lik;
  vector[N] prob_seropos;
  
  for (i in 1:N) {
    int birth_year = survey_year[i] - age[i];
    real prob_pos = 0.0;
    
    if (include_seroreversion) {
      for (y in 1:age[i]) {
        int cal_year = birth_year + y - 1;
        int year_idx = cal_year - min_year + 1;
        if (year_idx >= 1 && year_idx <= n_years) {
          prob_pos = prob_pos * (1.0 - rho) + (1.0 - prob_pos) * foi[year_idx];
        }
      }
    } else {
      real prob_neg = 1.0;
      for (y in 1:age[i]) {
        int cal_year = birth_year + y - 1;
        int year_idx = cal_year - min_year + 1;
        if (year_idx >= 1 && year_idx <= n_years) {
          prob_neg = prob_neg * (1.0 - foi[year_idx]);
        }
      }
      prob_pos = 1.0 - prob_neg;
    }
    
    prob_pos = fmax(1e-10, fmin(1.0 - 1e-10, prob_pos));
    prob_seropos[i] = prob_pos;
    n_pos_pred[i] = binomial_rng(n_sample[i], prob_pos);
    log_lik[i] = binomial_lpmf(n_pos[i] | n_sample[i], prob_pos);
  }
  
  // Summary statistics
  vector[n_years] foi_per_1000 = foi * 1000;
  
  // Period average (useful for interpretation)
  real foi_recent = mean(foi[(n_years-5):n_years]);  // Last 6 years
}
