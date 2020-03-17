functions {
  real conv_exp_pdf(int x, real lambda) {
    real out;
    out = (lambda / 2) * exp(-lambda * fabs(x));
    return out;
  }
  
  real w_dist(int x, vector w, int t) {
    real out;
    
    if(x > 0){
      if(x < (t-1)){
        out = w[x];
      }
    }else{
      out = 0;
    }
    return out;
  }
  
}

data {
  int t; // number of time steps
  int <lower = 0> obs_imported[t]; // imported cases
  int <lower = 0> obs_local[t]; // local cases
  int tau; // length of window
  int n; // number of samples for serial interval parameters
  real si_mean_samp[n]; // lognormal mean
  real si_sd_samp[n]; // lognormal std deviation
  real lam_mean; // observed mean of delay distribution
  int N; // number of delay samples
  vector[N] low; // lower bound for delay
  vector[N] up; // upper bound for delay
}

parameters{
  real <lower = 0> R[t]; // Effective reproduction number over time
  real <lower = 0> phi; // Dispersion of negative binomial distribution
  real <lower = 0> lambda; // delay distribution exponential parameter
  real <lower = 0> alpha_si_mean;
  real <lower = 0> beta_si_mean;
  real <lower = 0> alpha_si_sd;
  real <lower = 0> beta_si_sd;
  real ln_location; // log-normal location parameter
  real <lower = 0> ln_scale; // log-normal scale parameter
  // real <lower = 0> si_mean;
  // real <lower = 0> si_sd;
}

model {
  real infectiousness[t]; // infectiousness at each timestep
  real v[t-1]; // convolution of two delay distributions
  real wv[t - 1]; // convolution of w and v
  vector[t-1] wvec;
  // real w[t - 1]; // discretised serial interval distribution

  // real si_var; // serial interval distribution variance

  
  /////////////////////
  // SERIAL INTERVAL //
  /////////////////////
  
  alpha_si_mean ~ normal(0, 1) T[0,];
  beta_si_mean ~ normal(0, 1) T[0,];
  alpha_si_sd ~ normal(0, 1) T[0,];
  beta_si_sd ~ normal(0, 1) T[0,];

  for(i in 1:n) {
    target += gamma_lpdf(si_mean_samp[i] | alpha_si_mean, beta_si_mean);
    target += gamma_lpdf(si_sd_samp[i] | alpha_si_sd, beta_si_sd);
  }
  
  ln_location ~ gamma(alpha_si_mean, beta_si_mean);
  ln_scale ~ gamma(alpha_si_sd, beta_si_sd);
  // si_mean ~ gamma(146.7815, 106.6926);
  // si_sd ~ gamma(45.31689, 79.98991);

  
  // Calculate location and scale parameters for lognormal serial interval distribution
  // si_var = si_sd^2;
  // ln_location = log((si_mean^2) / sqrt(si_var + si_mean^2));
  // ln_scale = sqrt(log((si_var / (si_mean^2)) + 1));
  // ln_location = si_mean;
  // ln_scale = si_sd;

  // Discretise serial interval distribution
  // for (i in 1:t -1){
  //   w[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
  // }
  
  ////////////////////////
  // DELAY DISTRIBUTION // 
  ////////////////////////
  
  // Estimate delay distribution from onset to confirmation
  lambda ~ uniform(1/(5*lam_mean),1/(0.2*lam_mean));

  for(i in 1:N){
    target += log(exponential_cdf(up[i] , lambda) - exponential_cdf(low[i] , lambda));
  }
  
  // Calculate convolution of double delay and serial interval
  for(i in 1:(t-1)){
    wvec[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
  }
  
  for(i in 1:(t-1)) {
    wv[i] = 0;
    for(j in 1:50) {
      // for(k in -j:50){
        wv[i] +=  w_dist(j, wvec, t) * conv_exp_pdf(i-j, lambda);
      // }
    }
  }
  
  for(i in 1:(t-1))
    wv[i] = wv[i] / sum(wv);
  
  
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
    }
  }

  R ~ gamma(1, 5); // Prior used by EpiEstim
  // R ~ gamma(1, 0.2);
  phi ~ normal(0, 1) T[0,];
}

// spare code for checking serial interval distribution
generated quantities {
  vector[t-1] out_w_vec;
  real out_wv[t - 1];
  real out_v[t - 1];
  real si_mean;
  real si_sd;
  
  si_mean = exp(ln_location + (ln_scale^2)/2);
  si_sd = sqrt(exp(2*ln_location + ln_scale^2)*(exp(ln_scale^2)-1));

  // Discretise "double delay" gamma distribution
  for(i in 1:(t-1)){
    out_w_vec[i] = lognormal_cdf(i + 0.5, ln_location, ln_scale) - lognormal_cdf(i  - 0.5, ln_location, ln_scale);
    out_v[i] = exponential_cdf(i + 0.5, lambda) - exponential_cdf(i  - 0.5, lambda);
  }
  // Calculate convolution of serial interval and difference between two exponential delays
  for(i in 1:(t-1)) {
    out_wv[i] = 0;
    for(j in 0:50) {
      out_wv[i] +=  w_dist(j, out_w_vec, t) * conv_exp_pdf(i-j, lambda);
    }
  }
  
    
  for(i in 1:(t-1))
    out_wv[i] = out_wv[i] / sum(out_wv);

// out_i2[1] = 0;
//   for (s in 2:t){
//     out_i2[s] = 0;
//     for (i in 1:(s - 1)){
//       out_i2[s] += (obs_imported[i] + obs_local[i]) * out_wv[s - i];
//     }
//   }

}
