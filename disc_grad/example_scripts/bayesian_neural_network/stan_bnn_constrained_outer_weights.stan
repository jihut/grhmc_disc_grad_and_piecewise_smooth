data {
  int <lower = 1> n;
  int <lower = 1> p;
  matrix[n, p] x;
  vector[n] y;
}

parameters {
  
  real gamma; // gamma = log(sigma ^ 2) 
  real alpha;
  vector[2] w_star; // w_star = log(w) in order to get w >=0
  vector[2] delta;
  vector[p] beta1; // weights in first node
  vector[p] beta2; // weights in second node
  
}

transformed parameters {
  
  real <lower = 0> sigma = exp(0.5 * gamma);
  
}

model {
  
  // Prior
  
  target += (exponential_lpdf(exp(0.5 * gamma) | 1) + 0.5 * gamma); // want sigma to have exponential prior with rate 1
  alpha ~ normal(0, 1);
  w_star ~ normal(0, 1);
  delta ~ normal(0, 1);
  beta1 ~ normal(0, 1);
  beta2 ~ normal(0, 1);
  
  // Neural network with two neurons and one hidden layer
  
  y ~ normal(alpha + exp(w_star[1]) * fmax(delta[1] + x * beta1, 0) + exp(w_star[2]) * fmax(delta[2] + x * beta2, 0), sigma);
  
}
