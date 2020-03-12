data {
  int t; // number of time steps
  int n_samp; // number of samples from serial interval
  int <lower = 0> obs_imported[t, n_samp]; // imported cases
  int <lower = 0> obs_local[t, n_samp]; // local cases
  int tau; // length of window
  real si_mean[n_samp]; // lognormal mean
  real si_sd[n_samp]; // lognormal std deviation
}

transformed data{
  real infectiousness[t, n_samp];
  real w[t - 1, n_samp];
  real ln_location[n_samp];
  real si_var[n_samp];
  real ln_scale[n_samp];
  
  for(k in 1:n_samp){
    
    infectiousness[1,k] = 0;
  
  
    // Calculate location and scale parameters for lognormal serial interval distribution
    si_var[k] = si_sd[k]^2;
    ln_location[k] = log((si_mean[k]^2) / sqrt(si_var[k] + si_mean[k]^2));
    ln_scale[k] = sqrt(log((si_var[k] / (si_mean[k]^2)) + 1));

    // Discretise serial interval distribution
    for (i in 1:t -1){
      w[i, k] = lognormal_cdf(i + 0.5, ln_location[k], ln_scale[k]) - lognormal_cdf(i  - 0.5, ln_location[k], ln_scale[k]);
    }

    // Calculate infectiousness at each timestep
    for (s in 2:t){
      infectiousness[s, k] = 0;
      for (i in 1:(s - 1)){
        infectiousness[s, k] += (obs_imported[i, k] + obs_local[i, k]) * w[s - i, k];
      }
    }
  }
}

parameters {
  real <lower = 0> R[t, n_samp]; // Effective reproduction number over time
  real <lower = 0> phi; // Dispersion of negative binomial distribution
}

model {

for(k in 1:n_samp){
  R[,k] ~ gamma(1, 5); // Prior used by EpiEstim
  for (s in (tau + 1):t){
    for (i in (s-tau + 1):s){
      target += neg_binomial_2_lpmf(obs_local[i, k] | R[s, k] * infectiousness[i, k], 1 / sqrt(phi));
      // target += poisson_lpmf(obs_local[i] | R[s] * infectiousness[i]);
    }
  }
}


  phi ~ normal(0, 1) T[0,];
}

// spare code for checking serial interval distribution
// generated quantities {
//   real ran;
//   ran = lognormal_rng(ln_location, ln_scale);
// }
