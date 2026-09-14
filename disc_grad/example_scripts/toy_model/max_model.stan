data {
  int<lower=2> d;
  real c;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real q1;
  real q2;
  vector[d - 2] q;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  q1 ~ normal(0, 1);
  q2 ~ normal(fmax(c * q1, 0), 1);
  if (d > 2) {
    q ~ normal(fmax(c * q1, 0), 1);
  }
}

