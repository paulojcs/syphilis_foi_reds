// =============================================================================
// PIECEWISE CONSTANT FORCE OF INFECTION MODEL
// =============================================================================
// FOI is constant within defined periods but can vary between periods
// Flexible: K periods defined by foi_index input
// =============================================================================

data {
  int<lower=1> N;                      // Number of age groups
  array[N] int<lower=1> n_sample;      // Sample size per age group
  array[N] int<lower=0> n_pos;         // Number seropositive per age group
  array[N] int<lower=1> age;           // Age (years) for each group
  array[N] int<lower=1> survey_year;   // Year of survey for each observation
  int<lower=1> n_years;                // Number of years for FOI estimation
  int<lower=1> min_year;               // First calendar year in FOI series
  int<lower=1> K;                      // Number of FOI periods
  array[n_years] int<lower=1, upper=K> foi_index;  // Period index for each year
}

parameters {
  vector<lower=0, upper=1>[K] foi;     // FOI for each period
}

model {
  // Priors - weakly informative
  foi ~ beta(1, 50);  // Prior mean ~0.02
  
  // Likelihood: Catalytic model
  for (i in 1:N) {
    // Calculate probability of seropositivity by integrating over exposure years
    // Person of age a in survey_year was born in (survey_year - a)
    // Could be infected in any year from birth to survey_year - 1
    
    int birth_year = survey_year[i] - age[i];
    real prob_neg = 1.0;  // Probability of remaining seronegative
    
    // Integrate over each year of life
    for (y in 1:age[i]) {
      int cal_year = birth_year + y - 1;  // Calendar year
      // Map calendar year to foi_index
      // foi_index[1] corresponds to oldest year (min_birth_year)
      int year_idx = cal_year - min_year + 1;
      
      if (year_idx >= 1 && year_idx <= n_years) {
        int period = foi_index[year_idx];
        prob_neg = prob_neg * (1.0 - foi[period]);
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
  
  // Year-specific FOI for plotting
  vector[n_years] foi_by_year;
  
  for (i in 1:N) {
    int birth_year = survey_year[i] - age[i];
    real prob_neg = 1.0;
    
    for (y in 1:age[i]) {
      int cal_year = birth_year + y - 1;
      int year_idx = cal_year - min_year + 1;
      
      if (year_idx >= 1 && year_idx <= n_years) {
        int period = foi_index[year_idx];
        prob_neg = prob_neg * (1.0 - foi[period]);
      }
    }
    
    real prob_pos = 1.0 - prob_neg;
    n_pos_pred[i] = binomial_rng(n_sample[i], prob_pos);
    log_lik[i] = binomial_lpmf(n_pos[i] | n_sample[i], prob_pos);
  }
  
  // Map FOI to each year
  for (y in 1:n_years) {
    foi_by_year[y] = foi[foi_index[y]];
  }
  
  // FOI per 1000 person-years
  vector[K] foi_per_1000 = foi * 1000;
}
