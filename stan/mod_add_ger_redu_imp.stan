functions {
  // Function to derive the likelihood of weekly hospitalization counts
  // given: number of new infections on each day (lambda_hosp_day)
  // probability of hospitalizationfor every day of infection (p_hosp_days)
  // distribution of time between infection and hospitalization (tau_hosp)
  vector likel_hosp_week(int n_days, int n_weeks, int hosp_cutoff,
                         vector tau_hosp, vector new_inf_days,
                         vector p_hosp_days) {
    vector[n_days] lambda_hosp_day;
    vector[n_weeks] lambda_hosp_week;
    for (t in 1 : n_days) {
      lambda_hosp_day[t] = sum((new_inf_days[max(1,
                                                 t - hosp_cutoff) : t])
                               .* (p_hosp_days[max(1,
                                                   t - hosp_cutoff) : t])
                               .* reverse(tau_hosp[1 : min(t, hosp_cutoff+1)]));
    }
    for (w in 1 : n_weeks) {
      lambda_hosp_week[w] = sum(lambda_hosp_day[(1 + ((w - 1) * 7)) : (
                                7 + ((w - 1) * 7))]);
    }
    return (lambda_hosp_week);
  }
  // rep_each function
  vector rep_each(vector x, int K) {
    int N = rows(x);
    vector[N * K] y;
    int pos = 1;
    for (n in 1 : N) {
      for (k in 1 : K) {
        y[pos] = x[n];
        pos += 1;
      }
    }
    return y;
  }
}
data {
  // SIR model and associated data
  int<lower=1> n_days; // number of days of covariates
  int<lower=1> n_regions; // number of analysed regions
  int<lower=1> n_cat; // number of behavioral categories
  array[n_regions] matrix[n_days, n_cat] beh_cat; // Designmatrix without intercept for the region specific mobility reports
  int<lower=1> n_mult; // number of mult. effects (3 for temp and alpha and delta variant)
  array[n_regions] matrix[n_days, n_mult] mult_cov;
  array[n_regions] int<lower=1> n_agegroup; // number of age-groups per region used for accessing relevant 
  int n_agegroup_max; // maximum number of age-groups, used to specify padded objects containing -Inf for irrelevant values
  array[n_regions] vector[n_agegroup_max] pop_r_age; // population size per region and age-group at baseline, padded object with -Inf for irrelevant age-groups
  array[n_regions] matrix[n_days, n_agegroup_max] n_vaccinated; // number of vaccinated individuals age age-group and region, padded object with -Inf for irrelevant age-groups
  array[n_regions] vector<upper=1>[n_agegroup_max] p_hosp_age; // hospitalization risk per age-group in region, padded object with -Inf for irrelevant age-groups
  real<lower=0> gamma;
  int<lower=0> hosp_cutoff; // maximum number of days between infection and hosp. in days
  vector<lower=0>[hosp_cutoff+1] tau_hosp; // distribution of time between infection and hospitalizations.
  // Weekly data
  int<lower=1> n_weeks; // Number of weeks
  array[n_regions, n_weeks] int hosp_data; // Weekly hospital data. One row per region
  // Prior params
  //real<lower=0.5> sd_phi;
}
parameters {
  array[n_regions] vector<lower=0, upper=0.5>[n_cat] phi_star;
  vector[n_mult] phi_mult;
  array[n_regions] vector<lower=0, upper=0.5>[n_weeks] u_rt;
  vector<lower=0>[n_regions] sd_u_rt;
}
transformed parameters {
  array[n_regions] matrix[n_days, n_agegroup_max] S; // the susceptibles in each region and age-group (padded)
  array[n_regions] matrix[n_days, n_agegroup_max] I; // infectious in each region and age-group (padded)
  array[n_regions] matrix[n_days, n_agegroup_max] R_n; // removed in each region and age-group (padded)
  array[n_regions] vector[n_days] R_v; // removed by vaccination
  
  array[n_regions] vector[n_days] beta; // beta_{t,r}
  array[n_regions] vector<lower=0>[n_days] new_inf_days; // New infections per day in region
  array[n_regions] vector<lower=0>[n_days] p_hosp_days; // Hospitalization prob. among infected per day in region
  
  array[n_regions] vector<lower=0>[n_weeks] lambda_hosp_week;
  
  // The regression of beta.
  for (r in 1 : n_regions) {
    beta[r] = (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7)) .* 
    exp(mult_cov[r] * phi_mult);
  }

  //initialization (need to be at least one infectious individuals in each region at time t=0)
  //0.1 per 100000 of the population in each region is initially infectious
  for (r in 1 : n_regions) {
    for (j in 1 : n_agegroup[r]) {
      I[r, 1, j] = pop_r_age[r, j] / 100000 * 0.1; // 0.1 per 100000 of population in region and age-groups
      S[r, 1, j] = pop_r_age[r, j] - I[r, 1, j]; // Initial number of susceptibles per age-group
      R_n[r, 1, j] = 0; // no removed 
    }
    R_v[r, 1] = 0;
  }
  
  // The SIR model
  profile("SIR-trans_data") {
  for (r in 1 : n_regions) {
    for (t in 2 : n_days) {
      for (j in 1 : n_agegroup[r]) {
        // Susceptibles at time t: Susceptibles at t-1 minus new infections, minus vaccinations among susceptibles, minus "imported" infections
        S[r, t, j] = S[r, t - 1, j]
                     - ((beta[r, t] ./ sum(pop_r_age[r, 1 : n_agegroup[r]]))
                        .* sum(I[r, t - 1, 1 : n_agegroup[r]])
                        .* S[r, t - 1, j])
                     - (n_vaccinated[r, t - 1, j]
                        .* (S[r, t - 1, j]
                            ./ (S[r, t - 1, j] + R_n[r, t - 1, j])))
                     - (S[r, t - 1, j] / 100000 * 0.1); // "Importation": 0.1 per 100000 susceptibles
        // Infections at time t: Infections at t-1 + new infected, minus recovered, plus "imported" infections  
        I[r, t, j] = I[r, t - 1, j]
                     + ((beta[r, t] ./ sum(pop_r_age[r, 1 : n_agegroup[r]]))
                        .* sum(I[r, t - 1, 1 : n_agegroup[r]])
                        .* S[r, t - 1, j])
                     - (gamma * I[r, t - 1, j])
                     + (S[r, t - 1, j] / 100000 * 0.1); // "Importation": 0.1 per 100000 susceptibles
        // Recovered natural at time t: Recovered natural at t-1, plus newly recovered among infected, minus vaccinated among recovered    
        R_n[r, t, j] = R_n[r, t - 1, j] + (gamma * I[r, t - 1, j])
                       - (n_vaccinated[r, t - 1, j]
                          .* (R_n[r, t - 1, j]
                              ./ (S[r, t - 1, j] + R_n[r, t - 1, j])));
      }
      // Recovered vaccinated at time t: Recovered vaccinated at t-1 plus new number of vaccinated.
      R_v[r, t] = R_v[r, t - 1]
                  + sum(n_vaccinated[r, t - 1, 1 : n_agegroup[r]]);
    }
  }
  }
  
  // Derive quantities from latent compartments
  profile("new_inf_p_hosp-trans_data") {
  for (r in 1 : n_regions) {
    // New infections from one day to the other
    new_inf_days[r, 1] = sum(I[r, 1, 1 : n_agegroup[r]]);
    // Hospitalization probability at day t
    p_hosp_days[r, 1] = sum(pop_r_age[r][1 : n_agegroup[r]]
                            .* p_hosp_age[r][1 : n_agegroup[r]])
                        / sum(pop_r_age[r][1 : n_agegroup[r]]);
    for (t in 2 : n_days) {
      // New infections: Difference in susceptibles yesterday and today (>0), minus vaccinations among susceptibles
      new_inf_days[r, t] = sum(S[r, t - 1, 1 : n_agegroup[r]])
                           - sum(S[r, t, 1 : n_agegroup[r]])
                           - (sum(n_vaccinated[r, t - 1, 1 : n_agegroup[r]]
                                  .* (S[r, t - 1, 1 : n_agegroup[r]]
                                      ./ (S[r, t - 1, 1 : n_agegroup[r]]
                                          + R_n[r, t - 1, 1 : n_agegroup[r]]))));
      // Hospitalisation probability among infected on day t: weighted sum of hospitalization probabilities in age groups by share of age 
      // groups among susceptibles at t-1. Considered variants (alpha and delta) increase hosp. prob by factor 1.9 (not for German model)  
      p_hosp_days[r, t] = ((S[r, t - 1, 1 : n_agegroup[r]]
                            * p_hosp_age[r][1 : n_agegroup[r]])
                           / sum(S[r, t - 1, 1 : n_agegroup[r]]));
    }
  }
  }
  profile("likelihood-trans_data") {
    // Hospitalization likelihood
    for (r in 1 : n_regions) {
      lambda_hosp_week[r] = likel_hosp_week(n_days, n_weeks, hosp_cutoff,
      tau_hosp, new_inf_days[r],
      p_hosp_days[r]);
    }
  }
}
model {
  // Prior distributions
  profile("sampling-prior_u_rt") {
    sd_u_rt ~ normal(0,0.1);
    for (r in 1 : n_regions) {
      u_rt[r] ~ normal(0, sd_u_rt[r]);
    }
  }
  profile("sampling-phi_star") {
    for (r in 1 : n_regions) {
    phi_star[r] ~ normal(0,0.1);
    }
  }
    profile("sampling-phi_mult") {
    phi_mult ~ normal(0, 0.25);
  }
  //Likelihoods hospitalizations
  profile("sampling-likelihood") {
    for (r in 1 : n_regions) {
      hosp_data[r, ] ~ poisson(lambda_hosp_week[r]);
    }
  }
}

generated quantities {
  array[n_regions] vector[n_days] r_t;
  array[n_regions] vector[n_days] share_beta_ex;
  array[n_regions] vector[n_days] share_beta_1_ex;
  array[n_regions] vector[n_days] share_beta_2_ex;
  array[n_regions] vector[n_days] share_beta_3_ex;
  array[n_regions] vector[n_days] share_beta_4_ex;

  array[n_regions] real share_beta_ex_mean;
  array[n_regions] real share_beta_1_ex_mean;
  array[n_regions] real share_beta_2_ex_mean;
  array[n_regions] real share_beta_3_ex_mean;
  array[n_regions] real share_beta_4_ex_mean;

  for (r in 1:n_regions) {
    for (t in 1:n_days) {
      r_t[r][t] = beta[r][t]/gamma*sum(S[r, t, 1 : n_agegroup[r]]) / sum(pop_r_age[r, 1 : n_agegroup[r]]);
    }
  }
  for (r in 1:n_regions) {
    share_beta_ex[r] = (sum(phi_star[r]) + beh_cat[r] * phi_star[r]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7));
    share_beta_ex_mean[r] = mean((sum(phi_star[r]) + beh_cat[r] * phi_star[r]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7)));
    
    share_beta_1_ex[r] = (phi_star[r][1] + beh_cat[r][,1] * phi_star[r][1]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7));
    share_beta_1_ex_mean[r] = mean((phi_star[r][1] + beh_cat[r][,1] * phi_star[r][1]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7)));
    
    share_beta_2_ex[r] = (phi_star[r][2] + beh_cat[r][,2] * phi_star[r][2]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7));
    share_beta_2_ex_mean[r] = mean((phi_star[r][2] + beh_cat[r][,2] * phi_star[r][2]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7)));
    
    share_beta_3_ex[r] = (phi_star[r][3] + beh_cat[r][,3] * phi_star[r][3]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7));
    share_beta_3_ex_mean[r] = mean((phi_star[r][3] + beh_cat[r][,3] * phi_star[r][3]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7)));
    
    share_beta_4_ex[r] = (phi_star[r][4] + beh_cat[r][,4] * phi_star[r][4]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7));
    share_beta_4_ex_mean[r] = mean((phi_star[r][4] + beh_cat[r][,4] * phi_star[r][4]) ./ 
    (sum(phi_star[r]) + beh_cat[r] * phi_star[r] + rep_each(u_rt[r], 7)));
    
  }
}
