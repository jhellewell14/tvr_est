data {
  int t; // number of time steps
  int <lower = 0> obs_imported[t]; // imported cases
  int <lower = 0> obs_local[t]; // local cases
  int tau; // length of window
  real si_mean; // lognormal mean
  real si_sd; // lognormal std deviation
  real lam_mean; // observed mean of delay distribution
  int N; // number of delay samples
  vector[N] low; // lower bound for delay
  vector[N] up; // upper bound for delay
}

transformed data {
  real w[t - 1]; // discretised serial interval distribution
  real ln_location; // log-normal location parameter
  real si_var; // serial interval distribution variance
  real ln_scale; // log-normal scale parameter
  
  // Calculate location and scale parameters for lognormal serial interval distribution
  si_var = si_sd^2;
  ln_location = log((si_mean^2) / sqrt(si_var + si_mean^2));
  ln_scale = sqrt(log((si_var / (si_mean^2)) + 1));

  // Discretise serial interval distribution
  for (i in 1:t -1){
    w[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
  }
}
parameters{
  real <lower = 0> R[t]; // Effective reproduction number over time
  real <lower = 0> phi; // Dispersion of negative binomial distribution
  real <lower = 0> lambda; // delay distribution exponential parameter
}

model {
  real infectiousness[t]; // infectiousness at each timestep
  real v[t-1]; // convolution of two delay distributions
  real wv[t - 1]; // convolution of w and v

  ////////////////////////
  // DELAY DISTRIBUTION // 
  ////////////////////////
  
  // Estimate delay distribution from onset to confirmation
  lambda ~ uniform(1/(5*lam_mean),1/(0.2*lam_mean));

  for(i in 1:N){
    target += log(exponential_cdf(up[i] , lambda) - exponential_cdf(low[i] , lambda));
  }
  
  // Discretise delay distribution
  for(i in 1:(t-1)){
    v[i] =  exponential_cdf(i + 0.5, lambda) - exponential_cdf(i  - 0.5, lambda);
  }
  
  // Calculate convolution of double delay and serial interval
  wv[1] = 0.000001; // fix this?
  for(s in 2:(t-1)){
    wv[s] = 0;
    for(j in 1:(s-1)){
          wv[s] += w[j]*v[s-j];
    }
  }
  
  /////////////////
  // ESTIMATE R0 //
  /////////////////
  
  // Calculate infectiousness at each timestep
  infectiousness[1] = 0;
  for (s in 2:t){
    infectiousness[s] = 0;
    for (i in 1:(s - 1)){
      infectiousness[s] += (obs_imported[i] + obs_local[i]) * wv[s - i];
    }
  }
  
  for (s in (tau + 1):t){
    for (i in (s-tau + 1):s){
      target += neg_binomial_2_lpmf(obs_local[i] | R[s] * infectiousness[i], 1 / sqrt(phi));
      // target += poisson_lpmf(obs_local[i] | R[s] * infectiousness[i]);
    }
  }

  R ~ gamma(1, 5); // Prior used by EpiEstim
  phi ~ normal(0, 1) T[0,];
}

// spare code for checking serial interval distribution
generated quantities {
  real out[t - 1];
  real out_v[t - 1];
  real out_w[t - 1];
  real out_i[t];
  real out_wv[t - 1];
  real out_i2[t];

  // Discretise "double delay" gamma distribution
  for(i in 1:(t-1)){
    out_v[i] =  exponential_cdf(i + 0.5, lambda) - exponential_cdf(i  - 0.5, lambda);
    out_w[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
  }
  
    // Calculate convolution of double delay and serial interval
  out_wv[1] = 0.000001; // fix this?
  for(s in 2:(t-1)){
    out_wv[s] = 0;
    for(j in 1:(s-1)){
          out_wv[s] += out_w[j]*out_v[s-j];
    }
  }
  
  for(i in 1:(t-1)){
    out_wv[i] = 0 ? 0.0000001 : out_wv[i];
  }
  
  
  out_i[1] = 0;
  for (s in 2:t){
    out_i[s] = 0;
    for (i in 1:(s - 1)){
      out_i[s] += (obs_imported[i] + obs_local[i]) * w[s - i];
    }
  }

out_i2[1] = 0;
  for (s in 2:t){
    out_i2[s] = 0;
    for (i in 1:(s - 1)){
      out_i2[s] += (obs_imported[i] + obs_local[i]) * out_wv[s - i];
    }
  }

}
