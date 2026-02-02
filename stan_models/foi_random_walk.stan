// =============================================================================
// SMOOTH RANDOM WALK FORCE OF INFECTION MODEL
// =============================================================================
// Year-by-year FOI with first-order random walk smoothing
// Most flexible model - captures fine temporal structure
// Smoothness controlled by sigma_rw parameter
// =============================================================================

data {
  int<lower=1> N;                      // Number of age groups
  array[N] int<lower=1> n_sample;      // Sample size per age group
  array[N] int<lower=0> n_pos;         // Number seropositive per age group
  array[N] int<lower=1> age;           // Age (years) for each group
  array[N] int<lower=1> survey_year;   // Year of survey for each observation
  int<lower=1> n_years;                // Number of years for FOI estimation
  int<lower=1> min_year;               // First calendar year in FOI series
  
  // Hyperparameters for smoothness
  real<lower=0> sigma_rw_prior_mean;   // Prior mean for RW sigma (recommend 0.02)
  real<lower=0> sigma_rw_prior_sd;     // Prior SD for RW sigma (recommend 0.01)
}

parameters {
  real log_foi_init;                   // Log FOI for first year
  vector[n_years - 1] log_foi_innovations;  // Random walk innovations
  real<lower=0> sigma_rw;              // Random walk standard deviation
}

transformed parameters {
  vector[n_years] log_foi;
  vector<lower=0, upper=1>[n_years] foi;
  
  // Build log-FOI via random walk
  log_foi[1] = log_foi_init;
  for (t in 2:n_years) {
    log_foi[t] = log_foi[t-1] + sigma_rw * log_foi_innovations[t-1];
  }
  
  // Transform to FOI scale (bounded 0-1)
  for (t in 1:n_years) {
    foi[t] = inv_logit(log_foi[t]);
  }
}

model {
  // Priors
  log_foi_init ~ normal(-5, 2);        // Prior: exp(-5) ≈ 0.007, allows 0.001-0.1
  log_foi_innovations ~ std_normal();  // Standard normal innovations
  sigma_rw ~ normal(sigma_rw_prior_mean, sigma_rw_prior_sd);  // Smoothness control
  
  // Likelihood: Catalytic model
  for (i in 1:N) {
    int birth_year = survey_year[i] - age[i];
    real prob_neg = 1.0;
    
    for (y in 1:age[i]) {
      int cal_year = birth_year + y - 1;
      int year_idx = cal_year - min_year + 1;
      
      if (year_idx >= 1 && year_idx <= n_years) {
        prob_neg = prob_neg * (1.0 - foi[year_idx]);
      }
    }
    
    real prob_pos = 1.0 - prob_neg;
    n_pos[i] ~ binomial(n_sample[i], prob_pos);
  }
}

generated quantities {
  // Posterior predictive
  array[N] int n_pos_pred;
  vector[N] log_lik;
  
  for (i in 1:N) {
    int birth_year = survey_year[i] - age[i];
    real prob_neg = 1.0;
    
    for (y in 1:age[i]) {
      int cal_year = birth_year + y - 1;
      int year_idx = cal_year - min_year + 1;
      
      if (year_idx >= 1 && year_idx <= n_years) {
        prob_neg = prob_neg * (1.0 - foi[year_idx]);
      }
    }
    
    real prob_pos = 1.0 - prob_neg;
    n_pos_pred[i] = binomial_rng(n_sample[i], prob_pos);
    log_lik[i] = binomial_lpmf(n_pos[i] | n_sample[i], prob_pos);
  }
  
  // FOI per 1000 person-years
  vector[n_years] foi_per_1000 = foi * 1000;
  
  // Average FOI for recent period (last 6 years)
  real foi_recent_mean = mean(foi[(n_years-5):n_years]);
  real foi_recent_per_1000 = foi_recent_mean * 1000;
}
