data { 
  int<lower=2> d; 
  real<lower=0> eta; 
  real<lower=0> delta; 
}
parameters { 
  real q1;
  real q2;
  vector[d - 2] x; 
}
model {
  x ~ std_normal();
  for (i in 1:(d - 2)) {
    target += eta * fmax(0, 1 - abs(x[i]) / delta);
  }
  q1 ~ std_normal();
  target += eta * fmax(0, 1 - abs(q1) / delta); 
  q2 ~ std_normal();
  target += eta * fmax(0, 1 - abs(q2) / delta); 
}
