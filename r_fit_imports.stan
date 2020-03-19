data {
  int t; // number of time steps
  int <lower = 0> obs_imported[t]; // imported cases
  int <lower = 0> obs_local[t]; // local cases
  int tau; // length of window
  real si_mean; // lognormal mean
  real si_sd; // lognormal std deviation
}

transformed data{
  real infectiousness[t];
  real w[t - 1];
  real ln_location;
  real si_var;
  real ln_scale;
  infectiousness[1] = 0;
  
  
  // Calculate location and scale parameters for lognormal serial interval distribution
  si_var = si_sd^2;
  ln_location = log((si_mean^2) / sqrt(si_var + si_mean^2));
  ln_scale = sqrt(log((si_var / (si_mean^2)) + 1));

  // Discretise serial interval distribution
  for (i in 1:t -1){
    w[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
  }

  // Calculate infectiousness at each timestep
  for (s in 2:t){
    infectiousness[s] = 0;
    for (i in 1:(s - 1)){
      infectiousness[s] += (obs_imported[i] + obs_local[i]) * w[s - i];
    }
  }
}

parameters{
  real <lower = 0> R[t]; // Effective reproduction number over time
  real <lower = 0> phi; // Dispersion of negative binomial distribution
}

model {

  for (s in (tau + 1):t){
    for (i in (s-tau + 1):s){
      target += neg_binomial_2_lpmf(obs_local[i] | R[s] * infectiousness[i], 1 / sqrt(phi));
      // target += poisson_lpmf(obs_local[i] | R[s] * infectiousness[i]);
    }
  }

  // Mean = 2.6, sd = 2.6
  // R ~ lognormal(0.7231051, 0.6817717);
  R ~ gamma(1, 0.2); // Prior used by EpiEstim

  phi ~ normal(0, 1) T[0,];
}

// spare code for checking serial interval distribution
generated quantities {
  real w1[t-1];
  
  for (i in 1:t -1){
    w1[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
  }
}
