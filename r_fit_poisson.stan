functions {
  real discr_si(int k, real mu, real sigma) {
    real a;
    real b;
    real k1;
    real k2;
    real res;
    
    a = ((mu) / sigma)^2;
    b = sigma^2 / (mu);
    
    k1 = (k - 1) >= 0 ? k - 1 : 0;
    k2 = (k - 2) >= 0 ? k - 2 : 0;
    
    res = k * gamma_cdf(k, a, 1 / b) + (k - 2) * gamma_cdf(k2, a, 1 / b) - 2 * (k - 1) * gamma_cdf(k1, a, 1 / b);
    res += a * b * (2 * gamma_cdf(k1, a + 1, 1 / b) - 
                          gamma_cdf(k2, a + 1, 1 / b) - gamma_cdf(k, a + 1, 1 / b));
                          
    return(max([res, 0]));                     
  }
}

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
  real a;
  real b;

  
  // calculate alpha and beta for gamma distribution
  a = (si_mean / si_sd)^2;
  b = (si_sd^2) / si_mean;

  // Discretise serial interval distribution
  for (i in 1:(t - 1)){
    // w[i] = discr_si(i - 1, si_mean, si_sd);
    w[i] = gamma_cdf(i + 0.5, a, 1 / b) - gamma_cdf(i - 0.5, a, 1 / b);
  }

  // Calculate infectiousness at each timestep
  infectiousness[1] = 0;
  for (s in 2:t){
    infectiousness[s] = 0;
    for (i in 1:(s - 1)){
      infectiousness[s] += (obs_imported[i] + obs_local[i]) * w[s - i];
    }
  }
}

parameters{
  real <lower = 0> R[t]; // Effective reproduction number over time
}

model {

  for (s in (tau + 1):t){
    for (i in (s-tau + 1):s){
      target += poisson_lpmf(obs_local[i] | R[s] * infectiousness[i]);
    }
  }

  R ~ gamma(1, 0.2);

}

// spare code for checking serial interval distribution
generated quantities {

}
