// =============================================================================
// CONSTANT FORCE OF INFECTION MODEL
// =============================================================================
// Assumes FOI is constant across all calendar years
// Simplest model - serves as null/baseline for comparison
// =============================================================================

data {
  int<lower=1> N;                    // Number of age groups
  array[N] int<lower=1> n_sample;    // Sample size per age group
  array[N] int<lower=0> n_pos;       // Number seropositive per age group
  array[N] int<lower=1> age;         // Age (years) for each group
  array[N] int<lower=1> survey_year;  // Year of survey for each observation
}

transformed data {
  int n_years = max(age);            // Number of years for FOI estimation
  // Each age corresponds to years of potential exposure
}

parameters {
  real<lower=0, upper=1> foi;        // Single constant FOI
}

model {
  // Prior on FOI - weakly informative
  foi ~ beta(1, 50);  // Prior mean ~0.02, allows 0-0.1 range
  
  // Likelihood: Catalytic model without seroreversion
  // P(seropositive | age) = 1 - (1 - foi)^age
  for (i in 1:N) {
    real prob_pos = 1.0 - pow(1.0 - foi, age[i]);
    n_pos[i] ~ binomial(n_sample[i], prob_pos);
  }
}

generated quantities {
  // Posterior predictive
  array[N] int n_pos_pred;
  vector[N] log_lik;
  
  for (i in 1:N) {
    real prob_pos = 1.0 - pow(1.0 - foi, age[i]);
    n_pos_pred[i] = binomial_rng(n_sample[i], prob_pos);
    log_lik[i] = binomial_lpmf(n_pos[i] | n_sample[i], prob_pos);
  }
  
  // FOI per 1000 person-years for interpretation
  real foi_per_1000 = foi * 1000;
}
