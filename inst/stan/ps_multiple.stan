// Copyright 2024 Oliver Eales
// Copyright 2021 Oliver Eales
// Copyright 2017 Milad Kharratzadeh
//
// Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


functions {
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order);
  vector build_b_spline(array[] real t, array[] real ext_knots, int ind, int order) {
    // INPUTS:
      //    t:          the points at which the b_spline is calculated
    //    ext_knots:  the set of extended knots
    //    ind:        the index of the b_spline
    //    order:      the order of the b-spline
    vector[size(t)] b_spline;
    vector[size(t)] w1 = rep_vector(0, size(t));
    vector[size(t)] w2 = rep_vector(0, size(t));
    if (order==1)
      for (i in 1:size(t)) // B-splines of order 1 are piece-wise constant
    b_spline[i] = (ext_knots[ind] <= t[i]) && (t[i] < ext_knots[ind+1]);
    else {
      if (ext_knots[ind] != ext_knots[ind+order-1])
        w1 = (to_vector(t) - rep_vector(ext_knots[ind], size(t))) /
          (ext_knots[ind+order-1] - ext_knots[ind]);
      if (ext_knots[ind+1] != ext_knots[ind+order])
        w2 = 1 - (to_vector(t) - rep_vector(ext_knots[ind+1], size(t))) /
          (ext_knots[ind+order] - ext_knots[ind+1]);
      // Calculating the B-spline recursively as linear interpolation of two lower-order splines
      b_spline = w1 .* build_b_spline(t, ext_knots, ind, order-1) +
        w2 .* build_b_spline(t, ext_knots, ind+1, order-1);
    }
    return b_spline;
  }
}



data {
  int num_data;             // number of data points
  int num_knots;            // num of knots
  int num_path;
  vector[num_knots] knots;  // the sequence of knots
  int spline_degree;        // the degree of spline (is equal to order - 1)
  int Y[num_data];
  int P[num_path, num_data];
  real X[num_data];
  int week_effect;          // Number of days in day of week effect? 1=none, 2=weekends?, 7=all days
  int DOW[num_data];        // integer of day of the week
  int<lower = 0, upper = 2> cov_structure; //0 is tau[1], 1 is tau[num_path], 2 is Sigma[num_path, num_path]
  int<lower = 0, upper = 1> noise_structure; //0 is only includes observation noise (same between pathogens), 1 includes noise in individual pathogens as well

  int phi_priors_provided;      // 1=priors not provided, 2=priors provided
  real<lower=0> phi_mean;
  real<lower=0> phi_sd;
  
  int tau_priors_provided;      // 1=priors not provided, 2=priors provided
  real<lower=0> tau_mean[cov_structure==0 ? 1: cov_structure==1? num_path: 0 ];
  real<lower=0> tau_sd[cov_structure==0 ? 1: cov_structure==1? num_path: 0 ];
}

transformed data {
  int num_basis = num_knots + spline_degree - 1; // total number of B-splines
  matrix[num_basis, num_data] B;  // matrix of B-splines
  vector[spline_degree + num_knots] ext_knots_temp;
  vector[2*spline_degree + num_knots] ext_knots; // set of extended knots
  ext_knots_temp = append_row(rep_vector(knots[1], spline_degree), knots);
  ext_knots = append_row(ext_knots_temp, rep_vector(knots[num_knots], spline_degree));
  for (ind in 1:num_basis)
    B[ind,:] = to_row_vector(build_b_spline(X, to_array_1d(ext_knots), ind, spline_degree + 1));

  int cols;
  int rows;

  if(noise_structure==1){
    cols = num_path;
    rows = num_data;

  }
  else{
    cols = 0;
    rows = 0;
  }
}

parameters {
  matrix[num_path, num_basis] a;

  matrix<lower=0>[cols,rows] c_new;
  real<lower=0> eta[noise_structure==1 ? 1: 0 ];

  real<lower=0> phi;

  real<lower=0> tau[cov_structure==0 ? 1: cov_structure==1? num_path: 0 ];

  cov_matrix[cov_structure==2? num_path: 0] Sigma;

  simplex[week_effect] day_of_week_simplex;
}


transformed parameters {
  matrix[num_path, num_data] a_new;

  for (i in 1:num_path)
    a_new[i,] =  a[i,]*B;

}

model {
  //// Priors
  // Second-order random walk prior on b-spline coefficients
  if(cov_structure==2){
    for(i in 3:num_basis)
      a[,i] ~ multi_normal(2*a[,(i-1)] - a[,(i-2)], Sigma);
  }

  else if(cov_structure==1){
    for(i in 1:num_path)
      a[i,3:num_basis] ~ normal(2*a[i,2:(num_basis-1)] - a[i,1:(num_basis-2)], tau[i]);
  }

  else{
    for(i in 1:num_path)
      a[i,3:num_basis] ~ normal(2*a[i,2:(num_basis-1)] - a[i,1:(num_basis-2)], tau[1]);
  }


  // Prior on phi
  if(phi_priors_provided ==2){
    phi ~ normal(phi_mean, phi_sd);
  }
  
  // Prior on tau
  if(tau_priors_provided ==2){
    tau ~ normal(tau_mean, tau_sd);
  }
  
  //// Likelihood
  // This assumes there is some noise in the number of symptomatic cases for each pathogen individually (with shared parameter eta)

  real y_model[num_data];

  if(noise_structure==1){
    for(i in 1:num_path){
      c_new[i,] ~ gamma(exp(a_new[i,])*eta[1], eta[1]);
    }

    for(i in 1:num_data){
      y_model[i] = sum(c_new[,i]);
      P[,i] ~ multinomial(c_new[,i]/y_model[i]);
    }
  }
  else{
    for(i in 1:num_data){
      y_model[i] = sum(exp(a_new[,i]));
      P[,i] ~ multinomial(exp(a_new[,i])/y_model[i]);
    }
  }


  // Total number of cases (Y[i]) is negative-binomially distributed
  // Proportion of each pathogen is multinomially distributed
  if(week_effect==1){
    for(i in 1:num_data){
      Y[i] ~ neg_binomial(y_model[i]*phi, phi);
      }
  }
  else{
    for(i in 1:num_data){
      Y[i] ~ neg_binomial(y_model[i]*phi*week_effect*day_of_week_simplex[DOW[i]], phi);
      }
  }

}

