data {
  int t; // number of time steps
  int <lower = 0> obs_imported[t]; // imported cases
  int <lower = 0> obs_local[t]; // local cases
  int tau; // length of window
}

parameters{
  real <lower = 0> R[t]; // Effective reproduction number over time
  real <lower = 0> phi; // Dispersion of negative binomial distribution
  real <lower = 0> si_mean;
  real <lower = 0> si_sd;
}

model {
  real infectiousness[t];
  real w[t - 1];
  real ln_location;
  real si_var;
  real ln_scale;
  
  // Calculate location and scale parameters for lognormal serial interval distribution
  si_var = si_sd^2;
  ln_location = log((si_mean^2) / sqrt(si_var + si_mean^2));
  ln_scale = sqrt(log((si_var / (si_mean^2)) + 1));
  infectiousness[1] = 0;
  
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

  for (s in (tau + 1):t){
    for (i in (s-tau + 1):s){
      target += neg_binomial_2_lpmf(obs_local[i] | R[s] * infectiousness[i], 1 / sqrt(phi));
      // target += poisson_lpmf(obs_local[i] | R[s] * infectiousness[i]);
    }
  }

  target +=  lognormal_lpdf(si_mean | 1.525425, 0.2104153) + lognormal_lpdf(si_sd | 1.008535, 0.3351887);
  si_mean ~ lognormal(1.525425, 0.2104153);
  si_sd ~ lognormal(1.008535, 0.3351887);
  

  // Mean = 2.6, sd = 2.6
  // R ~ lognormal(0.7231051, 0.6817717);
  R ~ gamma(1, 5); // Prior used by EpiEstim

  phi ~ normal(0, 1) T[0,];
}

// spare code for checking serial interval distribution
// generated quantities {
//   real ran;
//   ran = lognormal_rng(ln_location, ln_scale);
// }
