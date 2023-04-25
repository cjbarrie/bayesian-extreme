data{
int<lower = 1> n; // total number of observations
int<lower = 1> p; // number of covariates in design matrix
int<lower = 0> y[n]; // vector of labels
matrix[n, p] X; // design matrix
int<lower = 1> small_area_id[n];
int<lower = 1> N_small_area;
int<lower = 1> N_small_area_edges;
int<lower=1, upper=N_small_area> node1_small_area[N_small_area_edges];
int<lower=1, upper=N_small_area> node2_small_area[N_small_area_edges];
real scaling_factor;
int<lower = 1> large_area_id[n];
int<lower = 1> N_large_area; 
vector[N_large_area] log_offset;
matrix[2,N_large_area] theta;
}//end_data

parameters{//begin_param
vector[p] aux_a;
vector<lower = 0>[p] aux_b;


vector[N_small_area] phi;
vector[N_small_area] psi;
real<lower = 0,upper = 1> lambda;
real<lower = 0> sigma_gamma;
vector[N_large_area] eta;
 real<lower = 0> sigma_eta;
}//end_param

transformed parameters{//begin_transform
vector[p] beta = aux_a ./ sqrt(aux_b);

vector[n] mu;
vector[N_small_area] gamma = (sqrt(1-lambda) * phi + sqrt(lambda / scaling_factor) * psi)*sigma_gamma;

mu = log_offset[large_area_id] + eta[large_area_id]*sigma_eta + gamma[small_area_id] + X * beta; 
}//end_transform

model{//begin_model

aux_a ~ normal(0,1);
aux_b[1] ~ gamma(0.5,100*0.5);
aux_b ~ gamma(0.5,0.5);
target += -0.5 * dot_self(psi[node1_small_area] - psi[node2_small_area]);
phi ~ normal(0,1); 
sum(psi) ~ normal(0, 0.01 * N_small_area);
lambda ~ beta(0.5,0.5);
sigma_gamma ~ normal(0,1);
eta ~ normal(0,1);
sigma_eta ~ normal(0,1);
// likelihood
for (i in 1:n) { 
  target += log_mix(1-inv_logit(mu[i]),
                    bernoulli_lpmf(y[i] | theta[1,large_area_id[i]]),
                    bernoulli_lpmf(y[i] | theta[2,large_area_id[i]]));
                    }
}//end_model

generated quantities{//begin_gen
vector[n] log_lik;
vector[n] y_gen;
int r_gen[n];
  
for (i in 1:n) {
log_lik[i] = log_mix(1-inv_logit(mu[i]),bernoulli_lpmf(y[i] | theta[1,large_area_id[i]]),bernoulli_lpmf(y[i] | theta[2,large_area_id[i]]) );
r_gen[i] = bernoulli_rng(inv_logit(mu[i]));
y_gen[i] = bernoulli_rng(theta[r_gen[i]+1,large_area_id[i]]);
} }//end_gen