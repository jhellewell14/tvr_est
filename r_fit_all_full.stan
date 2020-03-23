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
  
  real wv_pdf(int x, vector w, int t, real lambda) {
    real out;
    
    out = 0;
    for(j in 0:50){
      out += w_dist(j, w, t) * conv_exp_pdf(x - j, lambda);
    }
    
    return out;
  }
  
}

data {
  int t; // number of time steps
  int <lower = 0> obs_imported[t]; // imported cases
  int <lower = 0> obs_local[t]; // local cases
  int tau; // length of window
  real lam_mean; // observed mean of delay distribution
  int N; // number of delay samples
  vector[N] low; // lower bound for delay
  vector[N] up; // upper bound for delay
  real <lower = 0> si_loc;
  real <lower = 0> si_scale;
}

transformed data {
  vector[t-1] wvec;
  // Discretised serial interval distribution
  for(i in 1:(t-1)){
    wvec[i] = lognormal_cdf(i + 0.5, si_loc, si_scale) - lognormal_cdf(i  - 0.5, si_loc, si_scale);
  }
}

parameters{
  real <lower = 0> R[t]; // Effective reproduction number over time
  real <lower = 0> phi; // Dispersion of negative binomial distribution
  real <lower = 0> lambda; // delay distribution exponential parameter
}

model {
  real infectiousness[t]; // infectiousness at each timestep
  real wv[t - 1]; // convolution of w and v
  real neg_wv1[t - 1];
  real zero_wv1;


  ////////////////////////
  // DELAY DISTRIBUTION // 
  ////////////////////////
  
  // Estimate delay distribution from onset to confirmation
  lambda ~ uniform(1/(5*lam_mean),1/(0.2*lam_mean));

  for(i in 1:N){
    target += log(exponential_cdf(up[i] , lambda) - exponential_cdf(low[i] , lambda));
  }
  
  // Convolution of serial interval and two delay distributions
  // For positive time values
  for(i in 1:(t-1)) {
    wv[i] = 0;
    for(j in 0:50) {
        wv[i] +=  w_dist(j, wvec, t) * conv_exp_pdf(i - j, lambda);
    }
  }
  // For negative time values
  for(i in 1:(t-1)) {
    neg_wv1[i] = 0;
    for(j in 0:50) {
      neg_wv1[i] += w_dist(j, wvec, t) * conv_exp_pdf(-i - j, lambda);
    }
  }
  // For time = 0
  zero_wv1 = 0;
  for(i in 0:50){
    zero_wv1 += w_dist(i, wvec, t) * conv_exp_pdf(i, lambda);
  }
  
  /////////////////
  // ESTIMATE R0 //
  /////////////////
  
  // Calculate infectiousness at each timestep
  for (s in 1:t){
    infectiousness[s] = (obs_imported[s] + obs_local[s]) * zero_wv1;
    for (i in 1:(s - 1)){
      infectiousness[s] += (obs_imported[i] + obs_local[i]) * wv[s - i];
    }
    for(k in (s+1):t) {
      infectiousness[s] += (obs_imported[k] + obs_local[k]) * neg_wv1[k - s];
    }
  }
  
  for (s in (tau + 1):t){
    for (i in (s-tau + 1):s){
      target += neg_binomial_2_lpmf(obs_local[i] | R[s] * infectiousness[i], 1 / sqrt(phi));
    }
  }

  R ~ gamma(1, 0.2); // Prior used by EpiEstim
  phi ~ normal(0, 1) T[0,];
}

// spare code for checking serial interval distribution
generated quantities {
  vector[t-1] out_w_vec;
  real out_wv[t - 1];
  real out_i[t];
  real out_i1[t];
  real neg_wv[t - 1]; // negative values for convolution of w and v
  real zero_wv;

  // Discretise "double delay" gamma distribution
  for(i in 1:(t-1)){
    out_w_vec[i] = lognormal_cdf(i + 0.5, si_loc, si_scale) - lognormal_cdf(i  - 0.5, si_loc, si_scale);
    // out_v[i] = exponential_cdf(i + 0.5, lambda) - exponential_cdf(i  - 0.5, lambda);
  }
  // // Calculate convolution of serial interval and difference between two exponential delays
  for(i in 1:(t-1)) {
    out_wv[i] = 0;
    for(j in 0:50) {
      out_wv[i] +=  w_dist(j, out_w_vec, t) * conv_exp_pdf(i-j, lambda);
    }
  }
  // x - j = i
  for(i in 1:(t-1)) { // CI = -i
    neg_wv[i] = 0;
    for(j in 0:50) { // SI = j, 
      neg_wv[i] += w_dist(j, out_w_vec, t) * conv_exp_pdf(-i - j, lambda);
    }
  }

  out_i[1] = 0;
  for (s in 2:t){
    out_i[s] = 0;
    for (i in 1:(s - 1)){
      out_i[s] += (obs_imported[i] + obs_local[i]) * out_wv[s - i];
    }
  }
  
  zero_wv = 0;
  for(i in 0:50){
    zero_wv += w_dist(i, out_w_vec, t) * conv_exp_pdf(i, lambda);
  }
  
  // out_i1[1] = 0;
  out_i1[1] = (obs_imported[1] + obs_local[1]) * zero_wv;
  
  for (s in 2:t){
    // out_i1[s] = 0;
    out_i1[s] = (obs_imported[s] + obs_local[s]) * zero_wv;
    
    // future infectiousness
    for(k in (s+1):t) {
      out_i1[s] += (obs_imported[k] + obs_local[k]) * neg_wv[k - s];
    }
    
    // past infectiousness 
    for (i in 1:(s - 1)){
      out_i1[s] += (obs_imported[i] + obs_local[i]) * out_wv[s - i];
    }

  }

}
