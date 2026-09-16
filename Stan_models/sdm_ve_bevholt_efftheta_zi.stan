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
    real<lower=0> Ai; //grid-cell area = 0.3*0.3 km2
    matrix[Nj,3] Eg; //effort gill (3 years: 2008,2009,2010)
    matrix[Nj,3] Rg; //log fisheries catch gil (log tons)  
    matrix[Nj,3] Ef; //effort fyke
    matrix[Nj,3] Rf; //log fisheries catch fyke 
    matrix[Nj,3] Et; //effort trawl
    matrix[Nj,3] Rt; //log fisheries catch trawl 
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
    // lengths need to be modified according to the case:
    int nonzero1g[9]; int nonzero1f[5]; int nonzero1t[2];
    int nonzero2g[11]; int nonzero2f[6]; int nonzero2t[3]; 
    int nonzero3g[11]; int nonzero3f[5]; int nonzero3t[3];
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
    vector<upper=0>[Dx1q] beta_bar_q; //no bottom type
    vector[Dx2] beta;
    vector<upper=0>[Dx2q] beta_q;
    vector[Dx3] b_bar;
    vector<upper=0>[Dx3q] b_bar_q;
    // GIGG covariance
    vector<lower =0>[3] tau;  vector<lower =0>[3] sigma_tau; 
    vector<lower=0>[Dx1] gamma1; vector<lower=0>[Dx2] gamma2; vector<lower=0>[Dx3] gamma3; //they are gamma^2
    vector<lower=0>[Dx1] lambda11; vector<lower=0>[Dx2] lambda21; vector<lower=0>[Dx3] lambda31; //they are lambda^2
    vector<lower=0>[Dx1] lambda12; vector<lower=0>[Dx2] lambda22; vector<lower=0>[Dx3] lambda32; //they are lambda^2
    
    
    real<lower=0> b; //larvae mortality rate 
    
    //spatail random effects
    real<lower=0> sigma_delta;
    real<lower=0> l;
    matrix[Nj,3] z;

    //NB overdispersion
    real<lower=0> r;
    
    //Fisheries
    //vector[Nj] log_w; 
    real<lower=0> lambda_g;
    real<lower=0> lambda_f;
    real<lower=0> lambda_t;
    real<lower=0> sigma_g;
    real<lower=0> sigma_f;
    real<lower=0> sigma_t;
    
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
    vector[Nj] sum_lambda_eff1;
    vector[Nj] sum_lambda_eff2;
    vector[Nj] sum_lambda_eff3;
    
    
     sum_lambda_eff1 = (Eg[,1]*lambda_g + Ef[,1]*lambda_f + Et[,1]*lambda_t  );
     sum_lambda_eff2 = (Eg[,2]*lambda_g + Ef[,2]*lambda_f + Et[,2]*lambda_t  );
     sum_lambda_eff3 = (Eg[,3]*lambda_g + Ef[,3]*lambda_f + Et[,3]*lambda_t  );

    
    Sigma = sigma_delta^2*exp(-Dist_spatial*inv(l))+4; // + 10
    L = cholesky_decompose(Sigma);
    delta[,1] = L*z[,1];
    delta[,2] = L*z[,2]; 
    delta[,3] = L*z[,3]; 

    theta = inv_logit(alpha_bar + x1*beta_bar + x1q*beta_bar_q );
    eta = exp( 2.475 + x2*beta + x2q*beta_q+ W' *to_vector(delta)) ;

    log_a = a_bar + x3*b_bar+ x3q*b_bar_q ; //prolif. rate for one vendace (weight=0.03 kg)
    eta_tilde = exp(log_a).* (eta)./ (1+ b* eta );     
    theta_new = inv_logit(alpha_bar + x_new1*beta_bar + x_new1q*beta_bar_q);
    eta_new = exp(2.475 + x_new2*beta + x_new2q*beta_q);
  }
  model {
    vector[Nj] Sum;
    vector[Nj] Sum_reg;

    // iid Gaussian spatial random effects
    log(l) ~ normal(4.2,1); // length scale: IC(0.95)=(10,500) km
    target += -log(l);
    sigma_delta ~ student_t(4,0,sqrt(1)); 
    to_vector(z) ~ normal(0, 1);
       
  
// A weakly informative prior for linear coeff
    alpha_bar ~ normal(2.475, 0.25); // alpha not included in GIGG
    a_bar ~ normal(5.975,1);  //exp(a_bar)= avg.nbr.eggs.survived, for a 1kg spawner=2000+-500
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
    b ~ student_t(4,0,20); //spawners mortality: small
    

    // A weakly informative prior for catch param
    lambda_g ~ gamma(2,10 );  
    lambda_f ~ gamma(2,10 ); 
    lambda_t ~ gamma(2,10 ); 
    sigma_g ~ student_t(4,0,sqrt(1)); 
    sigma_f ~ student_t(4,0,sqrt(1)); 
    sigma_t ~ student_t(4,0,sqrt(1)); 
  
    // A weakly informative prior for r
    r ~ gamma(1,.1 ); // overdispersion
    
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
    
    // Catches    (weight in tons)
     Sum = W_new*(theta_new );

     // Gill
    Rg[nonzero1g,1] ~ student_t(4,log( Ai/f*(weight[nonzero1g])/1000 .* Eg[nonzero1g,1]*lambda_g ./sum_lambda_eff1[nonzero1g] .*( W_new[nonzero1g,] * ((exp((W_new[nonzero1g,]' *(sum_lambda_eff1[nonzero1g] ./Sum[nonzero1g])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero1g,1])), sigma_g ); //2008
    Rg[nonzero2g,2] ~ student_t(4,log( Ai/f*(weight[nonzero2g])/1000 .* Eg[nonzero2g,2]*lambda_g ./sum_lambda_eff2[nonzero2g] .*( W_new[nonzero2g,] * ((exp((W_new[nonzero2g,]' *(sum_lambda_eff2[nonzero2g] ./Sum[nonzero2g])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero2g,2])), sigma_g ); //2009
    Rg[nonzero3g,3] ~ student_t(4,log( Ai/f*(weight[nonzero3g])/1000 .* Eg[nonzero3g,3]*lambda_g ./sum_lambda_eff3[nonzero3g] .*( W_new[nonzero3g,] * ((exp((W_new[nonzero3g,]' *(sum_lambda_eff3[nonzero3g] ./Sum[nonzero3g])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero3g,3])), sigma_g ); //2010
     // Fyke
    Rf[nonzero1f,1] ~ student_t(4,log( Ai/f*(weight[nonzero1f])/1000 .* Ef[nonzero1f,1]*lambda_f ./sum_lambda_eff1[nonzero1f] .*( W_new[nonzero1f,] * ((exp((W_new[nonzero1f,]' *(sum_lambda_eff1[nonzero1f] ./Sum[nonzero1f])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero1f,1])), sigma_f ); //2008
    Rf[nonzero2f,2] ~ student_t(4,log( Ai/f*(weight[nonzero2f])/1000 .* Ef[nonzero2f,2]*lambda_f ./sum_lambda_eff2[nonzero2f] .*( W_new[nonzero2f,] * ((exp((W_new[nonzero2f,]' *(sum_lambda_eff2[nonzero2f] ./Sum[nonzero2f])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero2f,2])), sigma_f ); //2009
    Rf[nonzero3f,3] ~ student_t(4,log( Ai/f*(weight[nonzero3f])/1000 .* Ef[nonzero3f,3]*lambda_f ./sum_lambda_eff3[nonzero3f] .*( W_new[nonzero3f,] * ((exp((W_new[nonzero3f,]' *(sum_lambda_eff3[nonzero3f] ./Sum[nonzero3f])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero3f,3])), sigma_f ); //2010
     // Trawl
    Rt[nonzero1t,1] ~ student_t(4,log( Ai/f*(weight[nonzero1t])/1000 .* Et[nonzero1t,1]*lambda_t ./sum_lambda_eff1[nonzero1t] .*( W_new[nonzero1t,] * ((exp((W_new[nonzero1t,]' *(sum_lambda_eff1[nonzero1t] ./Sum[nonzero1t])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero1t,1])), sigma_t ); //2008
    Rt[nonzero2t,2] ~ student_t(4,log( Ai/f*(weight[nonzero2t])/1000 .* Et[nonzero2t,2]*lambda_t ./sum_lambda_eff2[nonzero2t] .*( W_new[nonzero2t,] * ((exp((W_new[nonzero2t,]' *(sum_lambda_eff2[nonzero2t] ./Sum[nonzero2t])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero2t,2])), sigma_t ); //2009
    Rt[nonzero3t,3] ~ student_t(4,log( Ai/f*(weight[nonzero3t])/1000 .* Et[nonzero3t,3]*lambda_t ./sum_lambda_eff3[nonzero3t] .*( W_new[nonzero3t,] * ((exp((W_new[nonzero3t,]' *(sum_lambda_eff3[nonzero3t] ./Sum[nonzero3t])) .* (theta_new)) -1) .* (theta_new .* eta_new))) .* exp(delta[nonzero3t,3])), sigma_t ); //2010

  }
  
  
