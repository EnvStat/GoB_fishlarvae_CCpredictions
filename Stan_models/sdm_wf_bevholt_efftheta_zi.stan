data {
    int<lower=1> N; //sampling sites (NO tot. nbr.grid cells)
    int<lower=1> Dx1; //nbr. covariates theta model
    int<lower=1> Dx2; //nbr. covariates eta model
    int<lower=1> Dx3; //nbr. covariates a model
    int<lower=1> Dx1q; //nbr. quad. covariates theta model
    int<lower=1> Dx2q; //nbr. quad. covariates eta model
    int<lower=1> Dx3q; //nbr. quad. covariates a model
    matrix[N,Dx1] x1; //variables to consider for theta model
    matrix[N,Dx2] x2; //variables to consider for eta model
    matrix[N,Dx3] x3; //variables to consider for a model
    matrix[N,Dx1q] x1q; // quadratic variables to consider for theta model
    matrix[N,Dx2q] x2q; //quadratic variables to consider for eta model
    matrix[N,Dx3q] x3q; //quadratic variables to consider for a model
    int<lower=1> Ds; //2 dim
    int<lower=0> y[N];
    vector[N] V;
    int<lower=1> Nj; //nbr regions (18)
    matrix[Nj,3] E; //effort
    real<lower=0> Ai; //grid-cell area = 0.3*0.3 km2
    matrix[Nj,3] R; //log fisheries catch  (3 years) 
    int<lower=1> N_new; // tot. nbr.grid cells (for predictions)
    matrix[Nj*3,N] W; //ICES region nbr for each grid cell      
    matrix[Nj,N_new] W_new; //ICES region nbr for each grid cell    
    matrix[N_new,Dx1] x_new1;  //covariates at each grid cell, for theta_new
    matrix[N_new,Dx1q] x_new1q;  //covariates at each grid cell, for theta_new
    matrix[N_new,Dx2] x_new2;  //covariates at each grid cell, for eta_new
    matrix[N_new,Dx2q] x_new2q;  //covariates at each grid cell, for eta_new
    vector[Nj] weight; // mean weight (kg)
    matrix[Nj,N] Weggs; //ICES region nbr for each sampling site 
    real<lower =0, upper =1> f; //fraction of female whitefish
    matrix[Nj,2] s; //region coordinates (km)
    int nonzero1[Nj-4]; // only fin -> -4 , full GoB -> -17
    int nonzero2[Nj-4];
    int nonzero3[Nj-5];
    }
      transformed data {
    matrix[Nj, Nj] Dist_spatial;

    // off-diagonal elements
    for (i in 1:(Nj-1)) {
      for (j in (i+1):Nj) {
        Dist_spatial[i, j] = pow(dot_self(s[i] - s[j]),0.5)  ;    

        // Fill in the other half
        Dist_spatial[j, i] = Dist_spatial[i, j];
      }
    }
    // diagonal elements
    for (k in 1:Nj){
      Dist_spatial[k, k] = 1e-6; // add also some jitter
    }
  }
  
  parameters {
    //real alpha;
    real alpha_bar;
    real a_bar;
    //linear coefficients
    vector[Dx1] beta_bar;
    vector<upper=0>[Dx1q] beta_bar_q; 
    vector[Dx2] beta;
    vector<upper=0>[Dx2q] beta_q;
    vector[Dx3] b_bar;
    vector<upper=0>[Dx3q] b_bar_q;
    // GIGG covariance
    vector<lower =0>[3] tau;  vector<lower =0>[3] sigma_tau; 
    vector<lower=0>[Dx1] gamma1; vector<lower=0>[Dx2] gamma2; vector<lower=0>[Dx3] gamma3; 
    vector<lower=0>[Dx1] lambda11; vector<lower=0>[Dx2] lambda21; vector<lower=0>[Dx3] lambda31; 
    vector<lower=0>[Dx1] lambda12; vector<lower=0>[Dx2] lambda22; vector<lower=0>[Dx3] lambda32; 
    
    
    real<lower=0> b; //larvae mortality rate 
    
    //spatail random effects
    real<lower=0> sigma_delta;
    real<lower=0> l;
    matrix[Nj,3] z;

    //NB overdispersion
    real<lower=0> r;
    
    //Fisheries
    real<lower=0> lambda;
    real<lower=0> sigma_epsilon;
    
  }
  transformed parameters {
    matrix[Nj,3] delta;
    vector<lower=0, upper=1>[N] theta; //prob. of suitable spawning site
    vector<lower= 0>[N] eta; //spawners density  # machine_precision()
    vector<lower=0>[N] eta_tilde;  //larvae density
    vector[N] log_a;  //proliferation rate
    vector<lower=0, upper=1>[N_new] theta_new; //prob. of suitable spawning cell
    vector<lower=0>[N_new] eta_new; //spawners density in each cell
    matrix[Nj, Nj] Sigma;
    matrix[Nj, Nj] L;
    
    Sigma = sigma_delta^2*exp(-Dist_spatial*inv(l)) + 4;
    L = cholesky_decompose(Sigma);
    delta[,1] = L*z[,1];
    delta[,2] = L*z[,2]; 
    delta[,3] = L*z[,3]; 

    theta = inv_logit(alpha_bar + x1*beta_bar + x1q*beta_bar_q );
    eta = exp( 2.475 + x2*beta + x2q*beta_q+ W' *to_vector(delta)) ;

    log_a = a_bar + x3*b_bar+ x3q*b_bar_q ; 
    eta_tilde = exp(log_a).* ( Weggs' *to_vector(weight) ).* (eta)./(1+b* eta .*  (Weggs' *to_vector(weight ))); //Beverton-Holt model    
    theta_new = inv_logit(alpha_bar + x_new1*beta_bar + x_new1q*beta_bar_q);
    eta_new = exp(2.475 + x_new2*beta + x_new2q*beta_q);
  }
  model {
    vector[Nj] Sum;
    vector[Nj] Sum_theta;

    // iid Gaussian spatial random effects
    log(l) ~ normal(4.2,1); // length scale: IC(0.95)=(10,500) km
    target += -log(l);
    sigma_delta ~ student_t(4,0,sqrt(1)); 
    to_vector(z) ~ normal(0, 1);
       
  
  // A weakly informative prior for linear coeff
    //alpha ~ normal(0, sqrt(10) );
    alpha_bar ~ normal(2.475, 0.25); 
    a_bar ~ normal(9.145,1);  //exp(a_bar)= avg.nbr.eggs.survived, for a 1kg spawner=2000+-500 
    beta ~ normal(0, tau[1]*sqrt(gamma1 .*lambda11 )); 
    beta_bar ~ normal(0, tau[2]*sqrt(gamma2 .*lambda21 ));
    b_bar ~ normal(0, tau[3]*sqrt(gamma3 .*lambda31));
    beta_q ~ normal(0, tau[1]*sqrt(gamma1 .*lambda12 )); 
    beta_bar_q ~ normal(0, tau[2]*sqrt(gamma2 .*lambda22 ));
    b_bar_q ~ normal(0, tau[3]*sqrt(gamma3 .*lambda32));
    // GIGG parameters prior
    gamma1 ~ gamma(0.5,1); lambda11 ~ inv_gamma(0.5,1); lambda12 ~ inv_gamma(0.5,1); 
    gamma2 ~ gamma(0.5,1); lambda21 ~ inv_gamma(0.5,1); lambda22 ~ inv_gamma(0.5,1);
    gamma3 ~ gamma(0.5,1); lambda31 ~ inv_gamma(0.5,1); lambda32 ~ inv_gamma(0.5,1);
    tau ~ cauchy(0,sigma_tau); 
    sigma_tau ~ student_t(4,0,sqrt(1));
    // A weakly informative prior for b
    b ~ student_t(4,0,20); //spawners mortality: P(b>10)=0.01
    
    // A weakly informative prior for catch param
    lambda ~ gamma(2,10 );  //capture rate: 98%CI = (3x10-5,0.006)
    sigma_epsilon ~ student_t(4,0,sqrt(1)); 
  
    // A weakly informative prior for r
    r ~ gamma(1,0.1 ); // overdispersion 
    
    //Larvae abundances zero-inf NB
   for (n in 1:N) {

          if (y[n] == 0)
          target += log_sum_exp(bernoulli_lpmf(0 | theta[n]),
                              bernoulli_lpmf(1 | theta[n])
                                + neg_binomial_2_lpmf( y[n]| (V[n]*(eta_tilde[n])), r)); 
          else
          target += bernoulli_lpmf(1 | theta[n])
                    + neg_binomial_2_lpmf( y[n]| (V[n]*(eta_tilde[n])), r); //*10^6    
        
    }
    // Catches    
     Sum = W_new*(theta_new);
    
    R[nonzero1,1] ~ student_t(4,log(Ai/f*(weight[nonzero1])/1000 .*( W_new[nonzero1,] * ((exp((W_new[nonzero1,]' *(E[nonzero1,1]./Sum[nonzero1])) .* (theta_new)*lambda) -1) .* (theta_new .* eta_new)) .* exp(delta[nonzero1,1]))), sigma_epsilon ); //2008
    R[nonzero2,2] ~ student_t(4,log(Ai/f*(weight[nonzero2])/1000 .*( W_new[nonzero2,] * ((exp((W_new[nonzero2,]' *(E[nonzero2,2]./Sum[nonzero2])) .* (theta_new)*lambda) -1) .* (theta_new .* eta_new)) .* exp(delta[nonzero2,2]))), sigma_epsilon ); //2009
    R[nonzero3,3] ~ student_t(4,log(Ai/f*(weight[nonzero3])/1000 .*( W_new[nonzero3,] * ((exp((W_new[nonzero3,]' *(E[nonzero3,3]./Sum[nonzero3])) .* (theta_new)*lambda) -1) .* (theta_new .* eta_new)) .* exp(delta[nonzero3,3]))), sigma_epsilon ); //2010    

  }
