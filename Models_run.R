#######################
# Run Stan models results in paper: 
# "Effect of climate change on reproduction of sea-spawning coregonids in the Baltic Sea"
#######################

# Stan model fits ##

# Load data ####
library("rstan")
library(matrixStats)
options(mc.cores = parallel::detectCores())


setwd("~Data")
load("FishGoBData")  
cov1 = cov2 = cov3 = cov1q = cov2q = cov3q = 2:8

# Whitefish model run ####
#data
dat.time.q = list(Dx1 = length(cov1), 
                  Dx1q = length(cov1q), 
                  x1= x.fin[, cov1],
                  x1q= x.fin[, cov1q]^2, 
                  Dx2 = length(cov2), 
                  Dx2q = length(cov2q), 
                  x2= x.fin[, cov2],
                  x2q= x.fin[, cov2q]^2, 
                  Dx3 = length(cov3), 
                  Dx3q = length(cov3q), 
                  x3= x.fin[, cov3],
                  x3q= x.fin[, cov3q]^2, 
                  N = nrow(x.fin),
                  Ds = ncol(s.fin),
                  s = s.reg,
                  y = y[],
                  V = V[],
                  Nj=18,
                  E = E3,  #net-days
                  Ai = 300*300, #m2  
                  R = log(R3/1000), # tons 
                  region_indexes = reg.index, #reg.index
                  N_new =  nrow(xpred.fin),
                  W = W[,],
                  W_new = W_new[,],  # instead of region index we use the summation matrix
                  x_new1 = xpred.fin[,cov1],
                  x_new1q= xpred.fin[, cov1q]^2,
                  x_new2 = xpred.fin[,cov2],
                  x_new2q= xpred.fin[, cov2q]^2,
                  weight = weight.j, # kg 
                  Weggs= Weggs[,],
                  f=f.ratio, 
                  nonzero1 = which(R3[,1] !=0),
                  nonzero2 = which(R3[,2] !=0),
                  nonzero3 = which(R3[,3] !=0)
)
# Remove all unnecessary objects
rm(list = setdiff(ls(), "dat.time.q")) ; gc()

init1 = list("sigma_epsilon"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3), "lambda"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))
init2 = list("sigma_epsilon"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3), "lambda"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))
init3 = list("sigma_epsilon"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3), "lambda"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))
init4 = list("sigma_epsilon"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3), "lambda"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))

setwd("~Stan_models") 
file.name ="bevholt_efftheta_zi"
fit.time.quad = stan(file =  paste0("sdm_wf_",file.name, ".stan"), data = dat.time.q, warmup=200, iter = 700, thin = 1, chains = 4, init = list(init1, init2, init3, init4)) #warmup=3000, iter=5000 changed on 11/12/24
setwd("~Stan_models") # fit posteriors, save location
saveRDS(fit.time.quad, file = paste0("post_wf_",file.name, ".Rda"))


# Vendace model run ####
#data
dat.time.q = list(Dx1 = length(cov1), 
                  Dx1q = length(cov1q), 
                  x1 = x.fin[, cov1],
                  x1q = x.fin[, cov1q]^2, 
                  Dx2 = length(cov2), 
                  Dx2q = length(cov2q), 
                  x2 = x.fin[, cov2],
                  x2q = x.fin[, cov2q]^2, 
                  Dx3 = length(cov3), 
                  Dx3q = length(cov3q), 
                  x3 = x.fin[, cov3],
                  x3q = x.fin[, cov3q]^2, 
                  N = nrow(x.fin),
                  Ds = ncol(s.fin),
                  s = s.reg,
                  y = y.v[finnish.ind.sample],
                  V = V[finnish.ind.sample],
                  Nj = 18,
                  Eg = E3.ven.gil, Ef = E3.ven.fyke,Et = E3.ven.trawl,  #net-days
                  Ai = 300*300, #m2 
                  Rg = log(R3.ven.gil/1000),Rf = log(R3.ven.fyke/1000),Rt = log(R3.ven.trawl/1000), # tons 
                  region_indexes = reg.index, #reg.index
                  N_new = nrow(xpred.fin),
                  W = W[,],
                  W_new = W_new[,],  # instead of region index we use the summation matrix
                  x_new1 = xpred.fin[,cov1],
                  x_new1q= xpred.fin[, cov1q]^2,
                  x_new2 = xpred.fin[,cov2],
                  x_new2q = xpred.fin[, cov2q]^2,
                  weight = weight.v, # kg 
                  Weggs = Weggs[,],
                  f=f.ratio.v, 
                  nonzero1g = which(R3.ven.gil[,1] !=0),nonzero1f = which(R3.ven.fyke[,1] !=0),nonzero1t = which(R3.ven.trawl[,1] !=0),
                  nonzero2g = which(R3.ven.gil[,2] !=0),nonzero2f = which(R3.ven.fyke[,2] !=0),nonzero2t = which(R3.ven.trawl[,2] !=0),
                  nonzero3g = which(R3.ven.gil[,3] !=0),nonzero3f = which(R3.ven.fyke[,3] !=0),nonzero3t = which(R3.ven.trawl[,3] !=0)
)
# Remove all unnecessary objects
rm(list = setdiff(ls(), "dat.time.q")) ; gc()


init1 = list("sigma_g"=runif(1,0.1,0.3),"sigma_f"=runif(1,0.1,0.3),"sigma_t"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3), "lambda_t"=runif(1,0.01,0.05),"lambda_g"=runif(1,0.01,0.05),"lambda_f"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))
init2 = list("sigma_g"=runif(1,0.1,0.3),"sigma_f"=runif(1,0.1,0.3),"sigma_t"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3),  "lambda_t"=runif(1,0.01,0.05),"lambda_g"=runif(1,0.01,0.05),"lambda_f"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))
init3 = list("sigma_g"=runif(1,0.1,0.3),"sigma_f"=runif(1,0.1,0.3),"sigma_t"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3),  "lambda_t"=runif(1,0.01,0.05),"lambda_g"=runif(1,0.01,0.05),"lambda_f"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))
init4 = list("sigma_g"=runif(1,0.1,0.3),"sigma_f"=runif(1,0.1,0.3),"sigma_t"=runif(1,0.1,0.3), "sigma_delta"=runif(1,0.1,0.3),  "lambda_t"=runif(1,0.01,0.05),"lambda_g"=runif(1,0.01,0.05),"lambda_f"=runif(1,0.01,0.05),  "r"=runif(1,0.5,2),  "b"=runif(1,0.1,0.5), "l"=runif(1,3,5))

setwd("~Stan_models")
modelname ="bevholt_efftheta_zi" 
fit.time.quad = stan(file =  paste0("sdm_ve_",modelname,".stan"), data = dat.time.q, warmup=300, iter = 1000, thin = 3, chains = 4, init = list(init1, init2, init3, init4))
setwd("~Stan_models") # fit posteriors, save location
saveRDS(fit.time.quad, file = paste0("post_ve_",modelname, ".Rda"))
