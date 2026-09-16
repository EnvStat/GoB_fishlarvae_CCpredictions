#######################
# Outputs and Figures for paper: 
# "Effect of climate change on reproduction of sea-spawning coregonids in the Baltic Sea"
#######################

# Load libraries & data ####
library("rstan")
library("VGAM")
library(matrixStats)
library("boot")
library("raster")
library("vioplot")
library("viridis")
library(LaplacesDemon)
library(RColorBrewer)

# Function
source(functions.R)
# Data
setwd("~Data")
load("FishGoBData")

# weight matrix wf, ven
weight.mtr = cbind(weight.j, 1) # vendace: accounted for weight in the prior (since weight.v constant over ICES regions)
colnames(weight.mtr) = c("weight whitefish", "weight vendace")

# Graphical settings for raster maps
coldiv = colorspace::diverging_hcl(111,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
colcon = viridis(111)
covnames <- c("intercept",covnames1)

e <- extent(as.matrix(spred.fin))
rr <- raster(e, ncol=length(unique(spred.fin[,1])), nrow=length(unique(spred.fin[,2])))


# Load Posterior samples ####
# load wf & ve posteriors
setwd("~Stan_models") 

fit.list1 = vector(mode = "list", length = 2)
file.names =  c("post_wf_bevholt_efftheta_zi.Rda",
                "post_ve_bevholt_efftheta_zi.Rda")
for (f0 in 1:2) {
  fit.list1[[f0]] = readRDS(file.names[f0]); gc()
}

names(fit.list1) = c("Whitefish", "Vendace")


# Posterior check

par.wf = c( "sigma_epsilon", "sigma_delta", "lambda","alpha_bar","r","a_bar","b" ,"l" , "sigma_tau")
par.ve = c( "sigma_g","sigma_f","sigma_t", "sigma_delta", "lambda_g", "lambda_f", "lambda_t","alpha_bar","r","a_bar","b" ,"l" , "sigma_tau")
par1 = c("tau", "gamma1","gamma2","gamma3","lambda11","lambda21","lambda31","lambda12","lambda22","lambda32" )
par2 = c(apply(matrix(c( "beta", "beta_q", "beta_bar", "beta_bar_q", "b_bar", "b_bar_q"), ncol=1),1, function(x) paste0(x,"[",1:7,"]") ))
allpar = list("whitefish" = c(par.wf, par1, par2), "vendace" =  c(par.ve, par1, par2))
# load posterior samples & check
for(f0 in 1:2){
  smr = summary(fit.list1[[f0]], pars =allpar[[f0]], use_cache=F, probs =c(0.025,0.975))
  smrall= smr$summary
  print(summary(smrall[,6:7]))
  print(kableExtra::kable(round(smrall,4), "markdown", booktabs =TRUE))
}

# Figs S9, S10: Prior vs Posterior #####

prior.list = list("prior.wf" = data.frame(
  "sigma_epsilon"=abs(rst(1000,4, mu=0, sigma=sqrt(1))),
  "sigma_delta"=abs(rst(1000,4, mu=0, sigma=sqrt(1))),
  "lambda"=rgamma(1000,2,10),
  "alpha_bar"=rnorm(1000,2.475,0.25),
  "r"=rgamma(1000,1,.1),
  "a_bar"=rnorm(1000,9.145,1), 
  "b"=abs(rst(1000,4, mu=0, sigma=20)),
  "l" = exp(rnorm(1000,4.2,1)),
  "alpha" = rnorm(1000,2.475,2)), 
  "prior.ve" = data.frame(
    "sigma_g"=abs(rst(1000,4, mu=0, sigma=sqrt(1))),
    "sigma_f"=abs(rst(1000,4, mu=0, sigma=sqrt(1))),
    "sigma_t"=abs(rst(1000,4, mu=0, sigma=sqrt(1))),
    "sigma_delta"=abs(rst(1000,4, mu=0, sigma=sqrt(1))),
    "lambda_g"=rgamma(1000,2,10),
    "lambda_f"=rgamma(1000,2,10),
    "lambda_t"=rgamma(1000,2,10),
    "alpha_bar"=rnorm(1000,2.475,.25),
    "r"=rgamma(1000,1,.1),
    "a_bar"=rnorm(1000,5.975,1), 
    "b"=abs(rst(1000,4, mu=0, sigma=20)),
    "l" = exp(rnorm(1000,4.2,1)),
    "alpha" = rnorm(1000,2.475,2)))

# prior vs posterior plots 
for(f in 1:2){
  # alpha
  posts.delsig = as.matrix(fit.list1[[f]], pars = c("delta","Sigma") )
  delta = posts.delsig[,1:54]
  Sigma = posts.delsig[,55:ncol(posts.delsig)]
  alpha = cbind(getalpha((delta[,1:18]), (Sigma)),getalpha((delta[,19:36]), (Sigma)),getalpha((delta[,37:54]), (Sigma)))
  alpha = alpha[,1] +2.475 # 1 sample! apply(alpha,1, mean) !!!!add on 24.9.24, to account for correction due to quad.eff negative constrain 
  # parameters
  par.m = names(prior.list[[f]])
  # posterior
  m = as.matrix(fit.list1[[f]], pars = par.m[-length(par.m)] )  
  m = cbind(m, alpha)
  
  # prior
  prior = prior.list[[f]]
  mycol = c(rgb(.5,.5,.5,.5), rgb(1,.2,0,.5))
  q.par = rep(1, ncol(prior)); 
  if(f==2) q.par[c(9,11,12)] = c(.8,.95,.9) else  q.par[c(3,5)] = c(.9,.9)# xlim
  # image
  if(f==1) par(mfrow=c(3,3), mar = c(4,2,2,1)) else par(mfrow=c(4,4), mar = c(4,2,2,1))
  for (p in 1:ncol(prior)) {
    hist(prior[,p], freq =F, col=mycol[1], border=NA, xlab=paste(colnames(m)[p]), xlim=c(min(prior[,p], m[,p]),max(quantile(prior[,p],q.par[p]), quantile(m[,p], .7))), main ="", cex.lab=1.5, breaks=50)
    hist(m[,p], freq=F, col=mycol[2], border=NA, add=T, breaks=30)
    if(p==1) legend("topright", legend = c("prior", "posterior"), col=mycol, bty = "n", pch=15, cex=1.5)
    
  }
  
  par(mfrow=c(1,1), mar= c(5,4,4,2)+.1)
  
  # Table mean 50CI
  posterior.summary = cbind(apply(m,2, mean), apply(m,2,quantile,.25), apply(m,2, quantile, .75))
  colnames(posterior.summary) = c("mean", "0.25 quantile", "0.75 quantile")
  cat(kableExtra::kable(round(posterior.summary,4), "markdown", booktabs =TRUE))
}

# Fig S1: Vendace fishery data #####
# img.ratio: 681 x 696

#  barplot 
par(mfrow=c(2,1), mar= c(2, 4, 2, 2) + 0.1, mgp= c(3,1,0)) #c(5, 4, 4, 2) + 0.1
mycex =1.4

colors.bp = cbind(c(rgb(.8,.2,0,.3), rgb(1,.6,0,.3), rgb(0.4,.8,0.4,.3),rgb(0,.2,0.8,.3)),
                  c(rgb(.8,.2,0,.5), rgb(1,.6,0,.5), rgb(0.4,.8,0.4,.5),rgb(0,.2,0.8,.5)),
                  c(rgb(.8,.2,0,1), rgb(1,.6,0,1), rgb(0.4,.8,0.4,1),rgb(0,.2,0.8,1)))

Effort = list(E3.ven.gil+ E3.ven.fyke+ E3.ven.trawl, E3.ven.gil +E3.ven.fyke, E3.ven.gil)
Catch = list(R3.ven.gil + R3.ven.fyke + R3.ven.trawl, R3.ven.gil+ R3.ven.fyke, R3.ven.gil)
names(Effort) = names(Catch) = c("Trawl", "Fyke", "Gillnet")
addit = c(F,T,T)
eff.max = ceiling(max(log(Effort[[1]],10)))
catch.max = ceiling(max(Catch[[1]]))

#Effort 
for (e0 in 1:3) {
  E3 = Effort[[e0]]
  colnames(E3)=2008:2010
  logE3 = log(E3,10)
  logE3[which(!is.finite(logE3)) ]=0
  barplot(t(logE3), beside = TRUE,  border = NA,col=colors.bp[1:3, e0],las=1,cex.axis=mycex,las=3,ylim = c(0, eff.max),
          ylab="Effort (net-days)", cex.lab=mycex, yaxt ="n", add = addit[e0])
  axis(2, at = c(0, 1,2,3,4),  labels= c(0,10,100,1000,10000), cex.lab=mycex,cex.axis=mycex,)
}

#Catch
par(mar= c(4, 4, 0, 2)+.1)
for (r0 in 1:3) {
  R3 = Catch[[r0]]
  rownames(R3)=reg.index
  colnames(R3)=2008:2010
  barplot(t(R3/1000), beside = TRUE,  border = NA, col=colors.bp[1:3, r0],las=3,ylim = c(0, catch.max/1000),cex.axis=mycex,cex.names = mycex,
          xlab="region", ylab="Catch (Ton)", axis.lty=1, cex.lab=mycex,  add = addit[r0])
  #legend( "topright", bty="n", fill = colors.bp,legend = 2008:2011, title = "year:", cex = mycex)
  
}

# legend ###
# img.ratio: 944 x 635 
library(scales)
par(mfrow=c(2,1), mar = c(5, 4, 4, 2) + 0.1, mgp = c(2,2,2))
plot.new()
text(.35,-.5,2008, cex =1.5); text(.47,-.5,2009, cex =1.5); text(.59,-.5,2010, cex =1.5)

show_col(colors.bp[1:3,], F, borders = NA)
text(-.4,-.5,"Trawl", cex =1.5); text(-.4,-1.5,"Fyke", cex =1.5); text(-.4,-2.5,"Gillnet", cex =1.5)
#text(.5,.06,2008); text(1.5,.06,2009); text(2.5,.06,2010)


# Fig 1.D: Covariates maps #####
# img.ratio: 1775 x 521 

covnames1 = c("Distance to sand", "River influence" , "Distance to deep", "Exposure","Chlorophyll-a", "Winter ice","Salinity" )
units = c("index", "index", "m", "index", "index", "weeks","psu")
# Plot covariate maps
par(mfrow=c(1,7),mar= c(1, 1, 2.1, 6), oma = c(1,1,.1,.1), bty="L")

for (j0 in 1:7) {
  z.cov = rasterize(spred.fin, rr, xpred.nost[,j0], fun=mean)
  plot(z.cov,  xaxt='n', yaxt='n',legend.width=4,axes = FALSE, box=T, col=colcon, legend=T,  legend.args = list(text=units[j0], 3, cex = 1))
  title(main=covnames1[j0], cex.main=1.8)
  plot.coords()
}
par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))
rm(list = c("xpred.nost"))

# Fig S6: predictive maps absolute change in ice and salt in 2050/2090s-RCP 4.5/RCP 8.5 #########
# img.ratio: 1330 x 732 (for _50s85)  (1284 x 510 for each scenarios)

# Future data, year 2059, 2097: avg salt and ice over models A,B,D
setwd("~Data/saltice_gob_abdfut")
xpred.fut.fin4 = xpred.fut.fin8 <- xpred.fin
year.fut = c(2059, 2097)
rc4.ice = rc8.ice = rc4.salt = rc8.salt = matrix(ncol = 2, nrow = nrow(xpred.fin)) # 2 columns: 50s,90s

for(y0 in 1:2){
  models=c("_A002","_A005", "_B002","_B005", "_D002","_D005" )
  # Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
  # load Ice , salt
  icesalt4A = read.table(paste("gob_si_abdfut_raster", models[1],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt4B = read.table(paste("gob_si_abdfut_raster", models[3],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt4D = read.table(paste("gob_si_abdfut_raster", models[5],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  # Future ice and salt
  icesalt4 = (icesalt4A + icesalt4B + icesalt4D)/3
  icesalt8 = (icesalt8A + icesalt8B + icesalt8D)/3
  # update future ice and salt
  xpred.fut.fin4[,7:8] = t( apply( t(apply(icesalt4,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
  xpred.fut.fin8[,7:8] = t( apply( t(apply(icesalt8,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
  rm(list = c("icesalt4", "icesalt4A", "icesalt4B", "icesalt4D",
              "icesalt8", "icesalt8A", "icesalt8B", "icesalt8D"))
  
  # Ice Salt absolute change
  
  # Unstandardize
  xpred.fut.fin4.nost = xpred.fut.fin8.nost = xpred.fin.nost <- xpred.fin
  xpred.fut.fin4.nost[,2:8] = t(apply(t(apply(xpred.fut.fin4[,2:8], 1, '*', stdxcont)),1,'+' , mxcont ))
  xpred.fut.fin8.nost[,2:8] = t(apply(t(apply(xpred.fut.fin8[,2:8], 1, '*', stdxcont)),1,'+' , mxcont ))
  xpred.fin.nost[,2:8] = t(apply(t(apply(xpred.fin[,2:8], 1, '*', stdxcont)),1,'+' , mxcont ))
  
  # ice and salt absolute change: (x.fut-x.hist)
  rc4.ice[,y0] = (xpred.fut.fin4.nost[,7] - xpred.fin.nost[,7])#/xpred.fin.nost[,7]
  rc8.ice[,y0] = (xpred.fut.fin8.nost[,7] - xpred.fin.nost[,7])#/xpred.fin.nost[,7]
  rc4.salt[,y0] = (xpred.fut.fin4.nost[,8] - xpred.fin.nost[,8])#/xpred.fin.nost[,8]
  rc8.salt[,y0] = (xpred.fut.fin8.nost[,8] - xpred.fin.nost[,8])#/xpred.fin.nost[,8]
}

abs.ch.mtr = cbind(rc4.ice, rc8.ice, rc4.salt, rc8.salt)
colnames(abs.ch.mtr) = paste0("Abs. change ",c("ice RCP4.5 50s","ice RCP4.5 90s","ice RCP8.5 50s","ice RCP8.5 90s", "salt RCP4.5 50s", "salt RCP4.5 90s", "salt RCP8.5 50s", "salt RCP8.5 90s"))
kableExtra::kable(summary(abs.ch.mtr), format = "markdown", caption = "Absolute change in ice-cover and salinity in RCP4.5 and RCP8.5")


## Absolute change predictive maps

# Rasters 50s
z.rc4.ice50 <- rasterize(spred.fin, rr, rc4.ice[,1], fun=mean)
z.rc8.ice50 <- rasterize(spred.fin, rr, rc8.ice[,1], fun=mean)
z.rc4.salt50 <- rasterize(spred.fin, rr, rc4.salt[,1], fun=mean)
z.rc8.salt50 <- rasterize(spred.fin, rr, rc8.salt[,1], fun=mean)
# Rasters 90s
z.rc4.ice90 <- rasterize(spred.fin, rr, rc4.ice[,2], fun=mean)
z.rc8.ice90 <- rasterize(spred.fin, rr, rc8.ice[,2], fun=mean)
z.rc4.salt90 <- rasterize(spred.fin, rr, rc4.salt[,2], fun=mean)
z.rc8.salt90 <- rasterize(spred.fin, rr, rc8.salt[,2], fun=mean)
invisible(gc())

# Plot raster maps Absolute change ###
par(mfrow=c(1,4),mar= c(1, 1, 3.1, 6), oma = c(1,1,.1,.1), bty="L")

# Ice ##
# Range raster map colors
raster.range = max(abs(range(c(rc4.ice, rc8.ice))))
raster.breaks = seq(-raster.range[1], raster.range[1], length = length(coldiv)+1)
# 50s
plot(z.rc4.ice50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 4.5", sub = year.fut[1],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv,breaks = raster.breaks, legend=F, ylab = "Ice-cover, absolute change"); plot.coords()
#title(ylab =  "Ice-cover, absolute change", line=0)
plot(z.rc8.ice50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 8.5", sub = year.fut[1],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     breaks = raster.breaks,col=coldiv, legend=F); plot.coords()
# 90s
plot(z.rc4.ice90, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 4.5", sub = year.fut[2],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv,breaks = raster.breaks, legend=F); plot.coords()#, ylab = "Ice-cover, absolute change")
plot(z.rc8.ice90, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 8.5", sub = year.fut[2],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     breaks = raster.breaks,col=coldiv, legend=F); plot.coords()
# Legend
plot((z.rc4.ice50), breaks = raster.breaks, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=3,axes = FALSE,box=T,  
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend.only=T, legend.args = list(text="Change,\n weeks", cex = .5),
     axis.args = list(at = pretty(seq(-raster.range[1], raster.range[1], l=5)), labels=pretty(seq(-raster.range[1], raster.range[1], l=5))))

# Salt ##
# Range raster map colors
raster.range = max(abs(range(c(rc4.salt, rc8.salt))))
raster.breaks = seq(-raster.range[1], raster.range[1], length = length(coldiv)+1)
# 50s
plot(z.rc4.salt50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 4.5", sub = year.fut[1],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv,breaks = raster.breaks, legend=F, ylab = "Salinity, absolute change"); plot.coords()
#title(ylab =  "Salinity, absolute change", line=0)
plot(z.rc8.salt50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 8.5", sub = year.fut[1],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     breaks = raster.breaks,col=coldiv, legend=F); plot.coords()
# 90s
plot(z.rc4.salt90, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 4.5", sub = year.fut[2],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv,breaks = raster.breaks, legend=F); plot.coords()#, ylab = "Salinity, absolute change")
plot(z.rc8.salt90, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="RCP 8.5", sub = year.fut[2],
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     breaks = raster.breaks,col=coldiv, legend=F); plot.coords()
# Legend
plot((z.rc4.salt50), breaks = raster.breaks, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=3,axes = FALSE,box=T,  
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend.only=T, legend.args = list(text="Change,\n psu", cex = .5),
     axis.args = list(at = pretty(seq(-raster.range[1], raster.range[1], l=5)), labels=pretty(seq(-raster.range[1], raster.range[1], l=5))))
par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))

rm(list = c("xpred.fut.fin4.nost", "xpred.fut.fin8.nost", "xpred.fin.nost"))


# Fig 1.E, 1.F: Future values and absolute change for RCP 8.5 2050s #############
# img.ratio: 1775 x 629 (for _50s85)  (1284 x 510 for each scenarios)

# Future data, year 2059: avg salt and ice over models A,B,D
setwd("~Data/saltice_gob_abdfut")
xpred.fut.fin4 = xpred.fut.fin8 <- xpred.fin
year.fut = c(2059)
#rc4.ice = rc8.ice =  rc8.salt = ac8.salt = rc8.salt = matrix(ncol = 1, nrow = nrow(xpred.fin)) # 2 columns: 50s,90s

for(y0 in 1){
  models=c("_A002","_A005", "_B002","_B005", "_D002","_D005" )
  # Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
  # load Ice , salt
  icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut[y0],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
  # Future ice and salt
  icesalt8 = (icesalt8A + icesalt8B + icesalt8D)/3
  # update future ice and salt
  xpred.fut.fin8[,7:8] = t( apply( t(apply(icesalt8,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
  rm(list = c("icesalt8", "icesalt8A", "icesalt8B", "icesalt8D"))
  
  # Ice Salt absolute change
  
  # Unstandardize
  xpred.fut.fin8.nost = xpred.fin.nost <- xpred.fin
  
  xpred.fut.fin8.nost[,2:8] = t(apply(t(apply(xpred.fut.fin8[,2:8], 1, '*', stdxcont)),1,'+' , mxcont ))
  xpred.fin.nost[,2:8] = t(apply(t(apply(xpred.fin[,2:8], 1, '*', stdxcont)),1,'+' , mxcont ))
  
  # ice and salt fut values, relative and absolute change: (x.fut-x.hist)
  fut8.ice = xpred.fut.fin8.nost[,7]
  ac8.ice = (xpred.fut.fin8.nost[,7] - xpred.fin.nost[,7])
  rc8.ice = (xpred.fut.fin8.nost[,7] - xpred.fin.nost[,7])/xpred.fin.nost[,7]
  
  fut8.salt = xpred.fut.fin8.nost[,8]
  ac8.salt = (xpred.fut.fin8.nost[,8] - xpred.fin.nost[,8])
  rc8.salt = (xpred.fut.fin8.nost[,8] - xpred.fin.nost[,8])/xpred.fin.nost[,8]
}


# Rasters 50s
z.fut8.ice50 <- rasterize(spred.fin, rr, fut8.ice, fun=mean)
z.ac8.ice50 <- rasterize(spred.fin, rr, ac8.ice, fun=mean)
z.rc8.ice50 <- rasterize(spred.fin, rr, rc8.ice, fun=mean)
z.fut8.salt50 <- rasterize(spred.fin, rr, fut8.salt, fun=mean)
z.ac8.salt50 <- rasterize(spred.fin, rr, ac8.salt, fun=mean)
z.rc8.salt50 <- rasterize(spred.fin, rr, rc8.salt, fun=mean)

# Image
par(mfrow=c(1,6),mar= c(1, 1, 3.1, 6), oma = c(1,1,.1,.1), bty="L")
my.leg.wdth = 3
## Salinity
plot(z.fut8.salt50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=my.leg.wdth,axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=colcon, legend=T, main = "\n future values ")
plot.coords()

raster.range = max(abs(range(ac8.salt)))
raster.breaks = seq(-raster.range[1], raster.range[1], length = length(coldiv)+1)
plot(z.ac8.salt50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv, breaks = raster.breaks, main = "Salinity (psu)\n absolute change ", legend=F)
plot.coords()
# Legend
plot((z.ac8.salt50), breaks = raster.breaks, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=my.leg.wdth,axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend.only=T, 
     axis.args = list(at = pretty(seq(-raster.range[1], raster.range[1], l=5)), labels=pretty(seq(-raster.range[1], raster.range[1], l=5))))

raster.range = max(abs(range(rc8.salt)-1))
raster.breaks = seq(1-raster.range[1], 1+raster.range[1], length = length(coldiv)+1)
plot(z.rc8.salt50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv, breaks = raster.breaks, main = "\n relative change ", legend=F)
plot.coords()
# Legend
plot((z.rc8.salt50), breaks = raster.breaks, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=my.leg.wdth,axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend.only=T, 
     axis.args = list(at = pretty(seq(1-raster.range[1], 1+raster.range[1], l=5)), labels=pretty(seq(1-raster.range[1], 1+raster.range[1], l=5))))

## Ice
plot(z.fut8.ice50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=my.leg.wdth,axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=colcon, legend=T, main = "\n future values ")
plot.coords()

raster.range = max(abs(range(ac8.ice)))
raster.breaks = seq(-raster.range[1], raster.range[1], length = length(coldiv)+1)
plot(z.ac8.ice50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv, breaks = raster.breaks, main = "Last winter ice (weeks)\n absolute change ", legend=F)
plot.coords()
# Legend
plot((z.ac8.ice50), breaks = raster.breaks, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=my.leg.wdth,axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend.only=T, 
     axis.args = list(at = pretty(seq(-raster.range[1], raster.range[1], l=5)), labels=pretty(seq(-raster.range[1], raster.range[1], l=5))))

raster.range = max(abs(range(rc8.ice)-1))
raster.breaks = seq(1-raster.range[1], 1+raster.range[1], length = length(coldiv)+1)
plot(z.rc8.ice50, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), 
     col=coldiv, breaks = raster.breaks, main = "\n relative change ", legend=F)
plot.coords()
# Legend
plot((z.rc8.ice50), breaks = raster.breaks, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=my.leg.wdth,axes = FALSE,box=T, ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend.only=T, 
     axis.args = list(at = pretty(seq(1-raster.range[1], 1+raster.range[1], l=5)), labels=pretty(seq(1-raster.range[1], 1+raster.range[1], l=5))))

title(ylab = "RCP 8.5    2050s", outer = T, line = -1)



# Fig S7: Absolute change across climate models - RCP 8.5,2090s ################
# img.ratio: 1318 x 7751

## 2090s
# Future data, year 2097:  salt and ice over models A,B,D
# Data created in : futuredataset325.Rmd
setwd("~Data/saltice_gob_abdfut")
xpred.fut.fin8.90sA = xpred.fut.fin8.90sB = xpred.fut.fin8.90sD <- xpred.fin
year.fut = 2097

models=c("_A002","_A005", "_B002","_B005", "_D002","_D005" )
# Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
# load Ice , salt
icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]

# Unstandardize
xpred.fin.nost <- xpred.fin
xpred.fin.nost[,2:8] = t(apply(t(apply(xpred.fin[,2:8], 1, '*', stdxcont)),1,'+' , mxcont ))

# Ice
# ice absolute change for each climate model: (x.fut-x.hist)
ice8.90sA.absch = (icesalt8A[,1] - xpred.fin.nost[,7])
ice8.90sB.absch = (icesalt8B[,1] - xpred.fin.nost[,7])
ice8.90sD.absch = (icesalt8D[,1] - xpred.fin.nost[,7])

z.ice8.90sA.absch <- rasterize(spred.fin, rr, ice8.90sA.absch, fun=mean)
z.ice8.90sB.absch <- rasterize(spred.fin, rr, ice8.90sB.absch, fun=mean)
z.ice8.90sD.absch <- rasterize(spred.fin, rr, ice8.90sD.absch, fun=mean)

# Salinity
# salt absolute change for each climate model: (x.fut-x.hist)
salt8.90sA.absch = (icesalt8A[,2] - xpred.fin.nost[,8])
salt8.90sB.absch = (icesalt8B[,2] - xpred.fin.nost[,8])
salt8.90sD.absch = (icesalt8D[,2] - xpred.fin.nost[,8])

z.salt8.90sA.absch <- rasterize(spred.fin, rr, salt8.90sA.absch, fun=mean)
z.salt8.90sB.absch <- rasterize(spred.fin, rr, salt8.90sB.absch, fun=mean)
z.salt8.90sD.absch <- rasterize(spred.fin, rr, salt8.90sD.absch, fun=mean)


# Image
par(mfrow=c(1,4),mar= c(1, 4, 3, 3), oma = c(1,1,.1,.1), bty="L", mgp = c(2,1,0)) # , mar = c(2,5,1,.1)

# ice
# color scale fixed 
min.zetat = quantile(c(values((z.ice8.90sA.absch)),values((z.ice8.90sB.absch)),values((z.ice8.90sD.absch))),.005, na.rm=T)
max.zetat = quantile(c(values((z.ice8.90sA.absch)),values((z.ice8.90sB.absch)),values((z.ice8.90sD.absch))),1, na.rm=T) # q. changed to 1 on 151225
lim.z = max(max.zetat, abs(min.zetat))
breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
plot(z.ice8.90sA.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F, main = "Climate model A"); plot.coords()
title(ylab = "Winter ice RCP 8.5 2090s", line = 2.5, cex.lab= 1.5)
plot(z.ice8.90sB.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F, main = "Climate model B"); plot.coords()
plot(z.ice8.90sD.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F, main = "Climate model D"); plot.coords()
plot(z.brk, legend.only=TRUE, col=coldiv,
     legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
     axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                    labels=pretty.round(min(breaks), max(breaks)), 
                    cex.axis=1))

# salt
# color scale fixed 
min.zetat = quantile(c(values((z.salt8.90sA.absch)),values((z.salt8.90sB.absch)),values((z.salt8.90sD.absch))),.005, na.rm=T)
max.zetat = quantile(c(values((z.salt8.90sA.absch)),values((z.salt8.90sB.absch)),values((z.salt8.90sD.absch))),1, na.rm=T) # q. changed to 1 on 151225
lim.z = max(max.zetat, abs(min.zetat))
breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
plot(z.salt8.90sA.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F, main = "Climate model A"); plot.coords()
title(ylab = "Salinity RCP 8.5 2090s", line = 2.5, cex.lab= 1.5)
plot(z.salt8.90sB.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F, main = "Climate model B"); plot.coords()
plot(z.salt8.90sD.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
     ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F, main = "Climate model D"); plot.coords()
plot(z.brk, legend.only=TRUE, col=coldiv,
     legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
     axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                    labels=pretty.round(min(breaks), max(breaks)), 
                    cex.axis=1))

invisible(gc())


# Fig 3.B: Covariate violinplots #############
# img.ratio: ??

## Violin plots covariates: load future data

setwd("~Data/saltice_gob_abdfut")
xpred.fut50.fin4 = xpred.fut50.fin8 = xpred.fut90.fin4 = xpred.fut90.fin8 <- xpred.fin
year.fut = c(2059, 2097)
models=c("_A002","_A005", "_B002","_B005", "_D002","_D005" )

# Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
# load Ice , salt 50s
icesalt4A = read.table(paste("gob_si_abdfut_raster", models[1],"_",year.fut[1],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4B = read.table(paste("gob_si_abdfut_raster", models[3],"_",year.fut[1],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4D = read.table(paste("gob_si_abdfut_raster", models[5],"_",year.fut[1],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut[1],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut[1],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut[1],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
# Future ice and salt
icesalt4 = (icesalt4A + icesalt4B + icesalt4D)/3
icesalt8 = (icesalt8A + icesalt8B + icesalt8D)/3
# update future ice and salt
xpred.fut50.fin4[,7:8] = t( apply( t(apply(icesalt4,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
xpred.fut50.fin8[,7:8] = t( apply( t(apply(icesalt8,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
rm(list = c("icesalt4", "icesalt4A", "icesalt4B", "icesalt4D",
            "icesalt8", "icesalt8A", "icesalt8B", "icesalt8D"))

# Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
# load Ice , salt 90s
icesalt4A = read.table(paste("gob_si_abdfut_raster", models[1],"_",year.fut[2],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4B = read.table(paste("gob_si_abdfut_raster", models[3],"_",year.fut[2],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4D = read.table(paste("gob_si_abdfut_raster", models[5],"_",year.fut[2],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut[2],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut[2],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut[2],".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
# Future ice and salt
icesalt4 = (icesalt4A + icesalt4B + icesalt4D)/3
icesalt8 = (icesalt8A + icesalt8B + icesalt8D)/3
# update future ice and salt
xpred.fut90.fin4[,7:8] = t( apply( t(apply(icesalt4,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
xpred.fut90.fin8[,7:8] = t( apply( t(apply(icesalt8,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
rm(list = c("icesalt4", "icesalt4A", "icesalt4B", "icesalt4D",
            "icesalt8", "icesalt8A", "icesalt8B", "icesalt8D"))

## Image
# for each process (3 p0), for each covariate (7 i),
nbr.proc=3;nbr.cov=7; nbr.rangeval=100; n.mc.samples = 100
responses = array(dim= c(nbr.proc, nbr.cov, nbr.rangeval,n.mc.samples+1)) 

covnames1 = c("Distance to sand", "River influence" , "Distance to deep", "Exposure","Chlorophyll-a", "Winter ice","Salinity" )
XtrainingCat = x.fin[,1]
XtrainingCont = x.fin[,2:8]
XpredCont=xpred.fin[,2:8]
# Response curves range
Xminmax = t(apply(XpredCont,2,"range")) 
Xminmax[4,2] =2.533138 #q.9 for exposure
# Update ice and salt range to future
Xminmax[6,] = range(c(xpred.fin[,7], xpred.fut50.fin4[,7], xpred.fut50.fin8[,7], xpred.fut90.fin4[,7], xpred.fut90.fin8[,7])) # ice
Xminmax[7,] = range(c(xpred.fin[,8], xpred.fut50.fin4[,8], xpred.fut50.fin8[,8], xpred.fut90.fin4[,8], xpred.fut90.fin8[,8])) # salt

# Image ###
proces = c("change in\n suitability", "rel. change in\n spawners dens.", "rel. change in\n max prolif. rate")
rangey = matrix(c(-.5,0,.5, -2 , 0 , 2, -3,0,3), ncol=3, byrow = T) 
mycolor = c("grey95","grey", "grey30") # c("grey", "orange", "blue3")

par(oma = c(2,4,1,1), mfrow = c(1,7), mar = c(3, 1, 1, 1))
# Vioplots
for (i in 1:length(covnames1)) {
  xtemp = as.vector(seq(Xminmax[i,1],Xminmax[i,2],length.out = 100))
  # Current 
  vioplot(cbind(XtrainingCont[,i]),  col=mycolor[1], horizontal=T, xaxt='n', yaxt='n', ylim = Xminmax[i,], xlim=c(0.5,4.5),axes=F, side = "left", lty =2)
  vioplot(cbind(XpredCont[,i]), col=mycolor[1], horizontal=T,add=T,  xaxt='n', yaxt='n',  xlim=c(0.5,2.5),axes=F, side = "right")
  # Future ice and salt
  if(i==6|i==7) {
    # 50s
    vioplot(cbind(-(100:101),xpred.fut50.fin4[,i+1]), col=mycolor[2], horizontal=T,add=T,  xaxt='n', yaxt='n',  xlim=c(0.5,2.5),axes=F, side = "left")
    vioplot(cbind(-(100:101),xpred.fut50.fin8[,i+1]), col=mycolor[3], horizontal=T,add=T,  xaxt='n', yaxt='n',  xlim=c(0.5,2.5),axes=F, side = "right")
    # 90s
    vioplot(cbind(-(100:101),-(100:101),xpred.fut90.fin4[,i+1]), col=mycolor[2], horizontal=T,add=T,  xaxt='n', yaxt='n',  xlim=c(0.5,2.5),axes=F, side = "left")
    vioplot(cbind(-(100:101),-(100:101),xpred.fut90.fin8[,i+1]), col=mycolor[3], horizontal=T,add=T,  xaxt='n', yaxt='n',  xlim=c(0.5,2.5),axes=F, side = "right")
  }
  mtext(covnames1[i], side = 1,line = 3, cex=1)
  axis(1,at=xtemp[c(5,50,95)],labels=round(xtemp[c(5,50,95)]*stdxcont.fin[i]+mxcont.fin[i],0),cex.axis=1)
  if(i==1) axis(2, at =1:3,labels=c("2010s", "2050s", "2090s"), ylab = "Distribution of \n covariates value", las=2, cex.axis=1.4)
}
# legend
par(fig = c(0, 1, 0, 1), oma = c(1, 4, 7, 1), mar = c(3, 1, 0.01, 1), new = TRUE)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("topleft", fill = mycolor, c("Current", "Future, RCP4.5", "Future, RCP8.5"), bty = "n", cex=1.4)
legend("top", lty = 2:1, c("Sampling sites", "Full GoB"), bty = "n", cex=1.4)




# Fig: 3.A: response curves #######
# img.ratio: ?? 

# Covariates effect ###
# Covariates
covnames1 = c("Distance to sand", "River influence" , "Distance to deep", "Exposure","Chlorophyll-a", "Winter ice","Salinity" )
cov1 = cov2 = cov3 = cov1q = cov2q = cov3q = 2:8

prob=c(T,F,F)

# Response curve: computation
# for each model (4 f), for each process (3 p0), for each covariate (7 i),
nbr.fit=length(fit.list1);nbr.proc=3;nbr.cov=7; nbr.rangeval=100; n.mc.samples = 100
responses = array(dim= c(nbr.fit,nbr.proc, nbr.cov, nbr.rangeval,n.mc.samples+1)) #(f x p x j x x x n)
for(f in 1:length(fit.list1)){
  fit.time.quad = fit.list1[[f]]
  # linear coeff posteriors
  # load post param
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q] # prolif. rate
  
  B1.all = list(betabar1, beta1, bbar1)
  B2.all = list(betabar2, beta2, bbar2)
  for(p0 in 1:3) {
    B1 = B1.all[[p0]]
    B2 = B2.all[[p0]]
    for (i in 1:ncol(B1)){
      # select 100 increasing values over covariate range
      xtemp = as.vector(seq(Xminmax[i,1],Xminmax[i,2],length.out = 100))
      ftemp = matrix(xtemp, ncol =1)%*%matrix(B1[,i], nrow=1) + matrix(xtemp^2, ncol =1)%*%matrix(B2[,i], nrow=1)
      ftemp.mean = xtemp * mean(B1[,i]) + xtemp^2 * mean(B2[,i])
      # select one random MCMC sample
      if(f==1&i==7){good.samples = which(B1[,i]<quantile(B1[,i],.975)&B1[,i]>quantile(B1[,i],.025))
      mc.samp = sample(good.samples, size = n.mc.samples)} else mc.samp = sample(1:nrow(B1), size = n.mc.samples)
      ftemp.select = cbind(ftemp[, mc.samp], ftemp.mean)
      # Inverse link-function
      # theta = resp curves matrix, dim = length.cov.range x nbr.mc.sample + 1 (last row = mean)
      if(prob[p0]) theta = inv.logit(ftemp.select) else theta = exp(ftemp.select)
      # centering
      for(j1 in 1:ncol(theta)) theta[,j1] = theta[,j1] - mean(theta[,j1])
      
      responses[f,p0, i,,] = theta
    }
  }
}


# Image ###
proces = c("change in\n suitability ", "rel. change in\n spawners dens. ", "rel. change in\n max prolif. rate ")
rangey = matrix(c(-.5,0,.5, -2 , 0 , 2, -5,0,5), ncol=3, byrow = T) 
mycolor = cbind(c(rgb(.6,.6,.6,.3), rgb(.9,.6,.2,.3)), c("black", "red")) # gray & orange

par(oma = c(2,4,1,1), mfrow = c(3,7), mar = c(1, 1, 1, 1))
# For each process
for(p0 in 1:3) {
  for (i in 1:length(covnames1)){
    xtemp = as.vector(seq(Xminmax[i,1],Xminmax[i,2],length.out = 100))
    # plot frame
    mylim = rangey[p0,c(1,3)] 
    plot(xtemp, xtemp, col="white", ylim=mylim, yaxt="n", xaxt="n", xlab = "", ylab ="", bty ='n')
    if(i==1) mtext(proces[p0], side = 2, line = 2, cex=1)
    ticks = rangey[p0,]
    if(i==1) axis(2,at=ticks,labels=ticks,las=2, cex.axis=1.4) 
    # curves
    for(f in 1:length(fit.list1)){
      theta = responses[f,p0, i,,]  
      for (j1 in 1:ncol(theta)){
        lines(xtemp, theta[,j1], col=mycolor[f,1])
      }}
    for(f in 1:length(fit.list1)){
      theta = responses[f,p0, i,,]  
      # theta from post.expected reg.coeff
      theta.mean = theta[,ncol(theta)]
      lines(xtemp, theta.mean, col = mycolor[f,2], lwd = 1.5)}; 
    abline(h=0, col=1, lty=2)
  }}
# legend
par(fig = c(0, 1, 0, 1), oma = c(1, 1, .1, 1), mar = c(3, 1, 0.01, 1), new = TRUE)
plot.new()
legend("top", names(fit.list1), lty = 1, col = mycolor[,2], horiz = T, bty ='n', cex=1.4)

par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))


# Figs. 4,5,S2,S3, S4: predictive maps of suitability, proliferation rate, spawner density, and larval density in current and future (abs.change) scenarios ######
# img.ratio: 1318 x 739

## 2050s ###
# Future data, year 2059: avg salt and ice over models A,B,D
setwd("/home/piailari/Documents/saltice_gob_abdfut")
xpred.fut.fin4.50s = xpred.fut.fin8.50s <- xpred.fin
year.fut = 2059

models=c("_A002","_A005", "_B002","_B005", "_D002","_D005" )
# Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
# load Ice , salt
icesalt4A = read.table(paste("gob_si_abdfut_raster", models[1],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4B = read.table(paste("gob_si_abdfut_raster", models[3],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4D = read.table(paste("gob_si_abdfut_raster", models[5],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
# Future ice and salt
icesalt4 = (icesalt4A + icesalt4B + icesalt4D)/3
icesalt8 = (icesalt8A + icesalt8B + icesalt8D)/3
# update future ice and salt
xpred.fut.fin4.50s[,7:8] = t( apply( t(apply(icesalt4,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
xpred.fut.fin8.50s[,7:8] = t( apply( t(apply(icesalt8,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
rm(list = c("icesalt4", "icesalt4A", "icesalt4B", "icesalt4D",
            "icesalt8", "icesalt8A", "icesalt8B", "icesalt8D"))

## 2090s ###
# Future data, year 2097: avg salt and ice over models A,B,D
setwd("/home/piailari/Documents/saltice_gob_abdfut")
xpred.fut.fin4.90s = xpred.fut.fin8.90s <- xpred.fin
year.fut = 2097

models=c("_A002","_A005", "_B002","_B005", "_D002","_D005" )
# Future Datasets : E_etrs89    N_etrs89      SALT_FUT      ICELAST_FUT
# load Ice , salt
icesalt4A = read.table(paste("gob_si_abdfut_raster", models[1],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4B = read.table(paste("gob_si_abdfut_raster", models[3],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt4D = read.table(paste("gob_si_abdfut_raster", models[5],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8A = read.table(paste("gob_si_abdfut_raster", models[2],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8B = read.table(paste("gob_si_abdfut_raster", models[4],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
icesalt8D = read.table(paste("gob_si_abdfut_raster", models[6],"_",year.fut,".txt", sep = ""), header = T)[finnish.ind.raster,4:3]
# Future ice and salt
icesalt4 = (icesalt4A + icesalt4B + icesalt4D)/3
icesalt8 = (icesalt8A + icesalt8B + icesalt8D)/3
# update future ice and salt
xpred.fut.fin4.90s[,7:8] = t( apply( t(apply(icesalt4,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
xpred.fut.fin8.90s[,7:8] = t( apply( t(apply(icesalt8,1,'-',mxcont[6:7])),1,'/',stdxcont[6:7] ) )
rm(list = c("icesalt4", "icesalt4A", "icesalt4B", "icesalt4D",
            "icesalt8", "icesalt8A", "icesalt8B", "icesalt8D"))

leg.pos.viop= c("topleft", "bottomleft")
coldiv = colorspace::diverging_hcl(111,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
cov1 = cov2 = cov3 = cov1q = cov2q = cov3q = 2:8

## Predictive maps of log10(Median(P)) ###
for(f in 1:2){
  fit.time.quad = fit.list1[[f]]
  ## load post param ###
  
  # linear coeff posteriors
  alphabar = as.matrix(fit.time.quad, "alpha_bar")
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  abar = as.matrix(fit.time.quad, "a_bar")
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q]
  
  delta = as.matrix(fit.time.quad, "delta")
  eta.pred1 = as.matrix(fit.time.quad, "eta_new")
  theta.pred1 = as.matrix(fit.time.quad, "theta_new")
  # add delta+alpha to eta
  eta.pred2008= exp(delta[,1:18]%*%W_new) * eta.pred1
  eta.pred2009= exp(delta[,19:36]%*%W_new) * eta.pred1
  eta.pred2010= exp(delta[,37:54]%*%W_new) * eta.pred1
  # average over time
  eta.pred1.m = (eta.pred2008+ eta.pred2009+eta.pred2010)/3#*round(theta.pred1) # suitable sites
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new) # [N.samp x n.loc]
  a1 = t(W_weight)*exp(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(as.matrix(abar)) +
                         xpred.fin[,cov3]%*%t(bbar1) +  xpred.fin[,cov3q]^2%*%t(bbar2))
  eggsurv1 = exp(xpred.fin[,cov3]%*%t(bbar1) + xpred.fin[,cov3q]^2%*%t(bbar2)) # eggs survival=x*b.bar
  b1= matrix(rep(as.matrix(fit.time.quad, "b"), nrow(xpred.fin)), nrow=nrow(xpred.fin), byrow=T) 
  eta.tilde.pred1 = t(a1)*eta.pred1.m/(1+W_weight*t(b1)*eta.pred1.m)
  
  # Median values - current, future 50s & 90s + RCP 4.5 & 8.5 ###
  
  # Suitability (theta)
  theta_pred_m1=apply(theta.pred1, 2, "median")
  z.t1 <- rasterize(spred.fin, rr, theta_pred_m1, fun=mean) # or theta.pred
  ## 50s
  theta.pred4.50s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin4.50s[,cov1]%*%t(betabar1) +  xpred.fut.fin4.50s[,cov1q]^2%*%t(betabar2)) )
  theta.pred8.50s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin8.50s[,cov1]%*%t(betabar1) +  xpred.fut.fin8.50s[,cov1q]^2%*%t(betabar2)) )
  theta_pred_m4=apply(theta.pred4.50s, 2, "median")
  theta_pred_m8=apply(theta.pred8.50s, 2, "median")
  z.t4.50s.absch <- rasterize(spred.fin, rr, theta_pred_m4 - theta_pred_m1, fun=mean)
  z.t8.50s.absch <- rasterize(spred.fin, rr, theta_pred_m8 - theta_pred_m1, fun=mean)
  ## 90s
  theta.pred4.90s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin4.90s[,cov1]%*%t(betabar1) +  xpred.fut.fin4.90s[,cov1q]^2%*%t(betabar2)) )
  theta.pred8.90s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin8.90s[,cov1]%*%t(betabar1) +  xpred.fut.fin8.90s[,cov1q]^2%*%t(betabar2)) )
  theta_pred_m4=apply(theta.pred4.90s, 2, "median")
  theta_pred_m8=apply(theta.pred8.90s, 2, "median")
  z.t4.90s.absch <- rasterize(spred.fin, rr, theta_pred_m4 - theta_pred_m1, fun=mean)
  z.t8.90s.absch <- rasterize(spred.fin, rr, theta_pred_m8 - theta_pred_m1, fun=mean)
  
  ## fraction of suitable sites
  theta_pred_m1_frac = apply(round(theta.pred1), 1, "sum")/ncol(theta.pred1)
  theta_pred_m4_frac_50 = apply(round(theta.pred4.50s), 1, "sum")/ncol(theta.pred1)
  theta_pred_m8_frac_50 = apply(round(theta.pred8.50s), 1, "sum")/ncol(theta.pred1)
  theta_pred_m4_frac_90 = apply(round(theta.pred4.90s), 1, "sum")/ncol(theta.pred1)
  theta_pred_m8_frac_90 = apply(round(theta.pred8.90s), 1, "sum")/ncol(theta.pred1)
  
  # Relative variation in proliferation rate (x*b.bar )
  ## 50s
  logeggsurv.pred4.50s = t(xpred.fut.fin4.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(bbar2)) 
  logeggsurv.pred8.50s = t(xpred.fut.fin8.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin8.50s[,cov3q]^2%*%t(bbar2)) 
  ## 90s
  logeggsurv.pred4.90s = t(xpred.fut.fin4.90s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.90s[,cov3q]^2%*%t(bbar2)) 
  logeggsurv.pred8.90s = t(xpred.fut.fin8.90s[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90s[,cov3q]^2%*%t(bbar2)) 
  
  a_pred_m1 = apply(eggsurv1, 1, "median") # transpose: n.loc x n.samples
  z.a1 <- rasterize(spred.fin, rr, log(a_pred_m1,10), fun=mean)
  a_pred_m4_50s=apply(exp(logeggsurv.pred4.50s), 2, "median")
  a_pred_m8_50s=apply(exp(logeggsurv.pred8.50s), 2, "median")
  a_pred_m4_90s=apply(exp(logeggsurv.pred4.90s), 2, "median")
  a_pred_m8_90s=apply(exp(logeggsurv.pred8.90s), 2, "median")
  z.a4.50s.absch <- rasterize(spred.fin, rr, a_pred_m4_50s - a_pred_m1, fun=mean)
  z.a8.50s.absch <- rasterize(spred.fin, rr, a_pred_m8_50s - a_pred_m1, fun=mean)
  z.a4.90s.absch <- rasterize(spred.fin, rr, a_pred_m4_90s - a_pred_m1, fun=mean)
  z.a8.90s.absch <- rasterize(spred.fin, rr, a_pred_m8_90s - a_pred_m1, fun=mean)
  
  # Proliferation rate (log(a) )
  ## 50s
  loga.pred4.50s = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                         xpred.fut.fin4.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(bbar2)) )
  loga.pred8.50s = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                         xpred.fut.fin8.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin8.50s[,cov3q]^2%*%t(bbar2)) )
  ## 90s
  loga.pred4.90s = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                         xpred.fut.fin4.90s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.90s[,cov3q]^2%*%t(bbar2)) )
  loga.pred8.90s = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                         xpred.fut.fin8.90s[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90s[,cov3q]^2%*%t(bbar2)) )
  
  # spawners density (log(eta) )
  # avg delta over 3 years, for the region to which each cell belongs
  deltaalpha= t(2.475 + log((exp(delta[,1:18]%*%W_new) + exp(delta[,19:36]%*%W_new)+ exp(delta[,37:54]%*%W_new))/3) ) # [ n.loc x n.samples ] 
  ## 50s
  # constant change rate over each grid-cell
  eta.pred4 = exp(t( deltaalpha+ xpred.fut.fin4.50s[,cov3]%*%t(beta1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(beta2)) )#*round(theta.pred4.50s)
  eta.pred8 = exp(t( deltaalpha+ xpred.fut.fin8.50s[,cov3]%*%t(beta1) +  xpred.fut.fin8.50s[,cov3q]^2%*%t(beta2)) )#*round(theta.pred8.50s)
  c4 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred4, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  c8 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred4.50s = matrix(log(c4), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred4)) + log(eta.pred4)
  logeta.pred8.50s = matrix(log(c8), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8)) + log(eta.pred8)
  ## 90s
  # constant change rate over each grid-cell
  eta.pred4 = exp(t( deltaalpha+ xpred.fut.fin4.90s[,cov3]%*%t(beta1) +  xpred.fut.fin4.90s[,cov3q]^2%*%t(beta2)) )#*round(theta.pred4.90s)
  eta.pred8 = exp(t( deltaalpha+ xpred.fut.fin8.90s[,cov3]%*%t(beta1) +  xpred.fut.fin8.90s[,cov3q]^2%*%t(beta2)) )#*round(theta.pred8.90s)
  c4 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred4, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  c8 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred4.90s = matrix(log(c4), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred4)) + log(eta.pred4)
  logeta.pred8.90s = matrix(log(c8), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8)) + log(eta.pred8)
  
  eta_pred_m1 = apply(eta.pred1.m, 2, "median")
  z.eta1 <- rasterize(spred.fin, rr, log(eta_pred_m1,10), fun=mean)
  eta_pred_m4_50s=apply(exp(logeta.pred4.50s), 2, "median")
  eta_pred_m8_50s=apply(exp(logeta.pred8.50s), 2, "median")
  eta_pred_m4_90s=apply(exp(logeta.pred4.90s), 2, "median")
  eta_pred_m8_90s=apply(exp(logeta.pred8.90s), 2, "median")
  z.eta4.50s.absch <- rasterize(spred.fin, rr, eta_pred_m4_50s - eta_pred_m1, fun=mean)
  z.eta8.50s.absch <- rasterize(spred.fin, rr, eta_pred_m8_50s - eta_pred_m1, fun=mean)
  z.eta4.90s.absch <- rasterize(spred.fin, rr, eta_pred_m4_90s - eta_pred_m1, fun=mean)
  z.eta8.90s.absch <- rasterize(spred.fin, rr, eta_pred_m8_90s - eta_pred_m1, fun=mean)
  
  
  # larvae density eta.tilde
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new)
  ## 50s
  eta.tilde.pred4.50s =exp(loga.pred4.50s+logeta.pred4.50s-log(1+W_weight*t(b1)*exp(logeta.pred4.50s)))
  eta.tilde.pred8.50s =exp(loga.pred8.50s+logeta.pred8.50s-log(1+W_weight*t(b1)*exp(logeta.pred8.50s)))
  ## 90s 
  eta.tilde.pred4.90s =exp(loga.pred4.90s+logeta.pred4.90s-log(1+W_weight*t(b1)*exp(logeta.pred4.90s)))
  eta.tilde.pred8.90s =exp(loga.pred8.90s+logeta.pred8.90s-log(1+W_weight*t(b1)*exp(logeta.pred8.90s)))
  
  
  etatilde_pred_m1 = apply(eta.tilde.pred1, 2, "median")
  z.etat1 <- rasterize(spred.fin, rr, log(etatilde_pred_m1,10), fun=mean)
  etatilde_pred_m4_50s=apply(eta.tilde.pred4.50s, 2, "median")
  etatilde_pred_m8_50s=apply(eta.tilde.pred8.50s, 2, "median")
  etatilde_pred_m4_90s=apply(eta.tilde.pred4.90s, 2, "median")
  etatilde_pred_m8_90s=apply(eta.tilde.pred8.90s, 2, "median")
  z.etat4.50s.absch <- rasterize(spred.fin, rr, etatilde_pred_m4_50s - etatilde_pred_m1, fun=mean)
  z.etat8.50s.absch <- rasterize(spred.fin, rr, etatilde_pred_m8_50s - etatilde_pred_m1, fun=mean)
  z.etat4.90s.absch <- rasterize(spred.fin, rr, etatilde_pred_m4_90s - etatilde_pred_m1, fun=mean)
  z.etat8.90s.absch <- rasterize(spred.fin, rr, etatilde_pred_m8_90s - etatilde_pred_m1, fun=mean)
  
  
  # Plots only ###
  invisible(gc())
  # New maps
  par(mfrow=c(1,4),mar= c(1, 4, 3, 3), oma = c(1,1,.1,.1), bty="L", mgp = c(2,1,0)) # , mar = c(2,5,1,.1)
  
  # RCP 4.5 ##
  
  # Suitability
  # color scale fixed 
  breaks <- seq(0,1, length.out = length(coldiv)+1)
  
  plot(z.t1, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="Current", ylab="",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F, cex.main = 1.5)
  title(ylab = "Sites suitability", line = 2.5, cex.lab= 1.5, outer = F)

  plot((z.t1), legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,
       axis.args=list(at=c(.2,.4,.6,.8),
                      labels=c(.2,.4,.6,.8),
                      cex.axis=1)); plot.coords()
  # Suitability abs.change
  # color scale fixed 
  min.ztheta = quantile(c(values((z.t4.50s.absch)),values((z.t8.50s.absch)), values((z.t4.90s.absch)),values((z.t8.90s.absch))),.005, na.rm=T)
  max.ztheta = max(c(values((z.t4.50s.absch)),values((z.t8.50s.absch)), values((z.t4.90s.absch)),values((z.t8.90s.absch))), na.rm=T)
  lim.z = max(max.ztheta, abs(min.ztheta))
  breaks <- seq(-lim.z,lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot(z.t4.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="2050s",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F, cex.main = 1.5); plot.coords()
  plot(z.t4.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="2090s",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F, cex.main = 1.5); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  vioplot.absch(theta_pred_m1_frac,theta_pred_m4_frac_50, theta_pred_m4_frac_90, log10=F, leg.pos = leg.pos.viop[f])
  
  
  # Relative variation in proliferation rate (log(a) )
  # color scale fixed 
  min.za = quantile(c(values((z.a1))),.000001, na.rm=T)
  max.za = max(c(values((z.a1))), na.rm=T)
  breaks <- seq(min.za, max.za, length.out = length(colcon)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot((z.a1), xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,  ylab = "",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),  breaks= breaks,col=colcon, legend=F); plot.coords()
  title(ylab = "Relative variation in proliferation rate", line = 2.5, cex.lab= 1.5,  outer = F)

  plot(z.brk, legend.only=TRUE, col=colcon,
       legend.width = 2,breaks= breaks,legend.args = list(text=expression("10"^x), 3),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)),  
                      cex.axis=1))
  # Rel var. proliferation rate abs.change
  # color scale fixed 
  min.za = quantile(c(values((z.a4.50s.absch)),values((z.a8.50s.absch)),values((z.a4.90s.absch)),values((z.a8.90s.absch))),.0001, na.rm=T)
  max.za = max(c(values((z.a4.50s.absch)),values((z.a8.50s.absch)),values((z.a4.90s.absch)),values((z.a8.90s.absch))), na.rm=T)
  lim.z = max(max.za, abs(min.za))
  breaks <- seq(-lim.z,lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot(z.a4.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.a4.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)),  
                      cex.axis=1))
  vioplot.absch(a_pred_m1,a_pred_m4_50s, a_pred_m4_90s, leg.pos = leg.pos.viop[f])
  
  
  # Spawners density (log(eta) )
  # color scale fixed 
  min.zeta = quantile(c(values((z.eta1))),.000001, na.rm=T)
  max.zeta = max(c(values((z.eta1))), na.rm=T)
  breaks <- seq(min.zeta, max.zeta, length.out = length(colcon)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot((z.eta1), xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,  ylab = "",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),  breaks= breaks,col=colcon, legend=F); plot.coords()
  title(ylab = "Spawners density", line = 2.5, cex.lab= 1.5,  outer = F)

  plot(z.brk, legend.only=TRUE, col=colcon,
       legend.width = 2,breaks= breaks,legend.args = list(text=expression("10"^x), 3),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  # Spawners density abs.change
  # color scale fixed 
  min.zeta = quantile(c(values((z.eta4.50s.absch)),values((z.eta8.50s.absch)), values((z.eta4.90s.absch)),values((z.eta8.90s.absch))),.0001, na.rm=T)
  max.zeta = max(c(values((z.eta4.50s.absch)),values((z.eta8.50s.absch)), values((z.eta4.90s.absch)),values((z.eta8.90s.absch))), na.rm=T)
  lim.z = max(max.zeta, abs(min.zeta))
  breaks <- seq(-lim.z,lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot(z.eta4.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.eta4.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  vioplot.absch(eta_pred_m1,eta_pred_m4_50s, eta_pred_m4_90s,leg.pos = leg.pos.viop[f])
  
  
  # Larvae density (log(eta.tilde) )
  # color scale fixed 
  min.zetat = quantile(c(values((z.etat1))),.001, na.rm=T)
  max.zetat = quantile(c(values((z.etat1))),.9995, na.rm=T)
  breaks <- seq(min.zetat, max.zetat, length.out = length(colcon)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot((z.etat1), xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,  ylab = "",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),  breaks= breaks,col=colcon, legend=F); plot.coords()
  title(ylab = "Larvae density", line = 2.5, cex.lab= 1.5, outer = F)

  plot(z.brk, legend.only=TRUE, col=colcon,
       legend.width = 2,breaks= breaks,legend.args = list(text=expression("10"^x), 3),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  # Larvae density, abs.change 
  # color scale fixed 
  min.zetat = quantile(c(values((z.etat4.50s.absch)),values((z.etat8.50s.absch)), values((z.etat4.90s.absch)),values((z.etat8.90s.absch))),.001, na.rm=T)
  max.zetat = quantile(c(values((z.etat4.50s.absch)),values((z.etat8.50s.absch)),values((z.etat4.90s.absch)),values((z.etat8.90s.absch))),.9995, na.rm=T)
  lim.z = max(max.zetat, abs(min.zetat))
  breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot(z.etat4.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.etat4.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  vioplot.absch(etatilde_pred_m1,etatilde_pred_m4_50s, etatilde_pred_m4_90s, leg.pos =leg.pos.viop[f])
  
  invisible(gc())
  
  #par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))
  
  
  # RCP 8.5 ## (Supplement)
  
  # Suitability
  # color scale fixed 
  breaks <- seq(0,1, length.out = length(coldiv)+1)
  
  plot(z.t1, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,main="Current", ylab="",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F, cex.main = 1.5); plot.coords()
  title(ylab = "Sites suitability", line = 2.5, cex.lab= 1.5)

  plot((z.t1), legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,
       axis.args=list(at=c(.2,.4,.6,.8),
                      labels=c(.2,.4,.6,.8),
                      cex.axis=1))
  # Suitability abs.change
  # color scale fixed 
  min.ztheta = quantile(c(values((z.t4.50s.absch)),values((z.t8.50s.absch)), values((z.t4.90s.absch)),values((z.t8.90s.absch))),.0001, na.rm=T)
  max.ztheta = max(c(values((z.t4.50s.absch)),values((z.t8.50s.absch)), values((z.t4.90s.absch)),values((z.t8.90s.absch))), na.rm=T)
  lim.z = max(max.ztheta, abs(min.ztheta))
  breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot(z.t8.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="2050s", cex.main = 1.5,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F); plot.coords()
  plot(z.t8.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="2090s", cex.main = 1.5,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  vioplot.absch(theta_pred_m1_frac,theta_pred_m8_frac_50, theta_pred_m8_frac_90, log10=F, rcp = "RCP 8.5", leg.pos =leg.pos.viop[f])
  
  
  # Relative variation in proliferation rate (log(a) )
  # color scale fixed 
  min.za = quantile(c(values((z.a1))),.000001, na.rm=T)
  max.za = max(c(values((z.a1))), na.rm=T)
  breaks <- seq(min.za, max.za, length.out = length(colcon)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot((z.a1), xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,  ylab = "",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),  breaks= breaks,col=colcon, legend=F); plot.coords()
  title(ylab = "Relative variation in proliferation rate", line = 2.5, cex.lab= 1.5)

  plot(z.brk, legend.only=TRUE, col=colcon,
       legend.width = 2,breaks= breaks,legend.args = list(text=expression("10"^x), 3),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)),  
                      cex.axis=1))
  # Relative variation in proliferation rate abs.change
  # color scale fixed 
  min.za = quantile(c(values((z.a4.50s.absch)),values((z.a8.50s.absch)),values((z.a4.90s.absch)),values((z.a8.90s.absch))),.0001, na.rm=T)
  max.za = max(c(values((z.a4.50s.absch)),values((z.a8.50s.absch)),values((z.a4.90s.absch)),values((z.a8.90s.absch))), na.rm=T)
  lim.z = max(max.za, abs(min.za))
  breaks <- seq(-lim.z,lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot(z.a8.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  
  plot(z.a8.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)),  
                      cex.axis=1))
  vioplot.absch(a_pred_m1,a_pred_m8_50s, a_pred_m8_90s, rcp = "RCP 8.5",leg.pos = leg.pos.viop[f])
  
  
  # Spawners density (log(eta) )
  # color scale fixed 
  min.zeta = quantile(c(values((z.eta1))),.000001, na.rm=T)
  max.zeta = max(c(values((z.eta1))), na.rm=T)
  breaks <- seq(min.zeta, max.zeta, length.out = length(colcon)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot((z.eta1), xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,  ylab = "",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),  breaks= breaks,col=colcon, legend=F); plot.coords()
  title(ylab = "Spawners density", line = 2.5, cex.lab= 1.5)

  plot(z.brk, legend.only=TRUE, col=colcon,
       legend.width = 2,breaks= breaks,legend.args = list(text=expression("10"^x), 3),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  # Spawners density abs.change
  # color scale fixed 
  min.zeta = quantile(c(values((z.eta4.50s.absch)),values((z.eta8.50s.absch)), values((z.eta4.90s.absch)),values((z.eta8.90s.absch))),.0001, na.rm=T)
  max.zeta = max(c(values((z.eta4.50s.absch)),values((z.eta8.50s.absch)), values((z.eta4.90s.absch)),values((z.eta8.90s.absch))), na.rm=T)
  lim.z = max(max.zeta, abs(min.zeta))
  breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot(z.eta8.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.eta8.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  vioplot.absch(eta_pred_m1,eta_pred_m8_50s, eta_pred_m8_90s,  rcp = "RCP 8.5", leg.pos =leg.pos.viop[f])
  
  
  # Larvae density (log(eta.tilde) )
  # color scale fixed 
  min.zetat = quantile(c(values((z.etat1))),.005, na.rm=T)
  max.zetat = quantile(c(values((z.etat1))),1, na.rm=T) # q. changed to 1 on 151225
  breaks <- seq(min.zetat, max.zetat, length.out = length(colcon)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot((z.etat1), xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,  ylab = "",
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),  breaks= breaks,col=colcon, legend=F); plot.coords()
  title(ylab = "Larvae density", line = 2.5, cex.lab= 1.5)

  plot(z.brk, legend.only=TRUE, col=colcon,
       legend.width = 2,breaks= breaks,legend.args = list(text=expression("10"^x), 3),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  # Larvae density, abs.change 
  # color scale fixed 
  min.zetat = quantile(c(values((z.etat4.50s.absch)),values((z.etat8.50s.absch)), values((z.etat4.90s.absch)),values((z.etat8.90s.absch))),.005, na.rm=T)
  max.zetat = quantile(c(values((z.etat4.50s.absch)),values((z.etat8.50s.absch)),values((z.etat4.90s.absch)),values((z.etat8.90s.absch))),1, na.rm=T) # q. changed to 1 on 151225
  lim.z = max(max.zetat, abs(min.zetat))
  breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot(z.etat8.50s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.etat8.90s.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  vioplot.absch(etatilde_pred_m1,etatilde_pred_m8_50s, etatilde_pred_m8_90s,  rcp = "RCP 8.5", leg.pos = leg.pos.viop[f])
  
  invisible(gc())
  par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))
}



# Fig. S8: predictive maps absolute change in future suitability, proliferation rate, and larval density, under climate models A,B,D and RCP 8.5-2090s ########################
# img.ratio: 1318 x 776

leg.pos.viop= c("topleft", "bottomleft")
coldiv = colorspace::diverging_hcl(111,h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7))
cov1 = cov2 = cov3 = cov1q = cov2q = cov3q = 2:8

# predictive maps of log10(Median(P))
for(f in 1:2){
  fit.time.quad = fit.list1[[f]]
  ## load post param ####
  
  # linear coeff posteriors
  alphabar = as.matrix(fit.time.quad, "alpha_bar")
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  abar = as.matrix(fit.time.quad, "a_bar")
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q]
  
  delta = as.matrix(fit.time.quad, "delta")
  eta.pred1 = as.matrix(fit.time.quad, "eta_new")
  theta.pred1 = as.matrix(fit.time.quad, "theta_new")
  # add delta+alpha to eta
  eta.pred2008= exp(delta[,1:18]%*%W_new) * eta.pred1
  eta.pred2009= exp(delta[,19:36]%*%W_new) * eta.pred1
  eta.pred2010= exp(delta[,37:54]%*%W_new) * eta.pred1
  # average over time
  eta.pred1.m = (eta.pred2008+ eta.pred2009+eta.pred2010)/3#*round(theta.pred1) # suitable sites
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new) # [N.samp x n.loc]
  a1 = t(W_weight)*exp(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(as.matrix(abar)) +
                         xpred.fin[,cov3]%*%t(bbar1) +  xpred.fin[,cov3q]^2%*%t(bbar2))
  eggsurv1 = exp(xpred.fin[,cov3]%*%t(bbar1) + xpred.fin[,cov3q]^2%*%t(bbar2)) # eggs survival=x*b.bar
  b1= matrix(rep(as.matrix(fit.time.quad, "b"), nrow(xpred.fin)), nrow=nrow(xpred.fin), byrow=T) 
  eta.tilde.pred1 = t(a1)*eta.pred1.m/(1+W_weight*t(b1)*eta.pred1.m)
  
  # Median values - current, future 90s + RCP 8.5 + clim.mod A,B,D ###
  
  # Suitability (theta)
  theta_pred_m1=apply(theta.pred1, 2, "median")
  z.t1 <- rasterize(spred.fin, rr, theta_pred_m1, fun=mean) # or theta.pred
  
  ## 90s
  theta.pred8.90sA = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                    xpred.fut.fin8.90sA[,cov1]%*%t(betabar1) +  xpred.fut.fin8.90sA[,cov1q]^2%*%t(betabar2)) )
  theta.pred8.90sB = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                    xpred.fut.fin8.90sB[,cov1]%*%t(betabar1) +  xpred.fut.fin8.90sB[,cov1q]^2%*%t(betabar2)) )
  theta.pred8.90sD = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                    xpred.fut.fin8.90sD[,cov1]%*%t(betabar1) +  xpred.fut.fin8.90sD[,cov1q]^2%*%t(betabar2)) )
  theta_pred_m8A=apply(theta.pred8.90sA, 2, "median")
  theta_pred_m8B=apply(theta.pred8.90sB, 2, "median")
  theta_pred_m8D=apply(theta.pred8.90sD, 2, "median")
  z.t8.90sA.absch <- rasterize(spred.fin, rr, theta_pred_m8A - theta_pred_m1, fun=mean)
  z.t8.90sB.absch <- rasterize(spred.fin, rr, theta_pred_m8B - theta_pred_m1, fun=mean)
  z.t8.90sD.absch <- rasterize(spred.fin, rr, theta_pred_m8D - theta_pred_m1, fun=mean)
  
  
  # Eggs survival (x*b.bar )
  ## 90s
  logeggsurv.pred8.90sA = t(xpred.fut.fin8.90sA[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90sA[,cov3q]^2%*%t(bbar2)) 
  logeggsurv.pred8.90sB = t(xpred.fut.fin8.90sB[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90sB[,cov3q]^2%*%t(bbar2)) 
  logeggsurv.pred8.90sD = t(xpred.fut.fin8.90sD[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90sD[,cov3q]^2%*%t(bbar2)) 
  
  a_pred_m1 = apply(eggsurv1, 1, "median") # transpose: n.loc x n.samples
  z.a1 <- rasterize(spred.fin, rr, log(a_pred_m1,10), fun=mean)
  a_pred_m8_90sA=apply(exp(logeggsurv.pred8.90sA), 2, "median")
  a_pred_m8_90sB=apply(exp(logeggsurv.pred8.90sB), 2, "median")
  a_pred_m8_90sD=apply(exp(logeggsurv.pred8.90sD), 2, "median")
  z.a8.90sA.absch <- rasterize(spred.fin, rr, a_pred_m8_90sA - a_pred_m1, fun=mean)
  z.a8.90sB.absch <- rasterize(spred.fin, rr, a_pred_m8_90sB - a_pred_m1, fun=mean)
  z.a8.90sD.absch <- rasterize(spred.fin, rr, a_pred_m8_90sD - a_pred_m1, fun=mean)
  
  # Proliferation rate (log(a) )
  ## 90s
  loga.pred8.90sA = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                          xpred.fut.fin8.90sA[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90sA[,cov3q]^2%*%t(bbar2)) )
  loga.pred8.90sB = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                          xpred.fut.fin8.90sB[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90sB[,cov3q]^2%*%t(bbar2)) )
  loga.pred8.90sD = log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                          xpred.fut.fin8.90sD[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90sD[,cov3q]^2%*%t(bbar2)) )
  
  # spawners density (log(eta) )
  # avg delta over 3 years, for the region to which each cell belongs
  deltaalpha= t(2.475 + log((exp(delta[,1:18]%*%W_new) + exp(delta[,19:36]%*%W_new)+ exp(delta[,37:54]%*%W_new))/3) ) # [ n.loc x n.samples ] 
  ## 90s
  # constant change rate over each grid-cell
  eta.pred8A = exp(t( deltaalpha+ xpred.fut.fin8.90sA[,cov3]%*%t(beta1) +  xpred.fut.fin8.90sA[,cov3q]^2%*%t(beta2)) )#*round(theta.pred4.90s)
  eta.pred8B = exp(t( deltaalpha+ xpred.fut.fin8.90sB[,cov3]%*%t(beta1) +  xpred.fut.fin8.90sB[,cov3q]^2%*%t(beta2)) )#*round(theta.pred8.90s)
  eta.pred8D = exp(t( deltaalpha+ xpred.fut.fin8.90sD[,cov3]%*%t(beta1) +  xpred.fut.fin8.90sD[,cov3q]^2%*%t(beta2)) )#*round(theta.pred8.90s)
  c8A = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8A, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  c8B = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8B, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  c8D = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8D, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred8.90sA = matrix(log(c8A), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8A)) + log(eta.pred8A)
  logeta.pred8.90sB = matrix(log(c8B), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8B)) + log(eta.pred8B)
  logeta.pred8.90sD = matrix(log(c8D), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8D)) + log(eta.pred8D)

  
  
  # larvae density eta.tilde
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new)
  ## 90s 
  eta.tilde.pred8.90sA =exp(loga.pred8.90sA+logeta.pred8.90sA-log(1+W_weight*t(b1)*exp(logeta.pred8.90sA)))
  eta.tilde.pred8.90sB =exp(loga.pred8.90sB+logeta.pred8.90sB-log(1+W_weight*t(b1)*exp(logeta.pred8.90sB)))
  eta.tilde.pred8.90sD =exp(loga.pred8.90sD+logeta.pred8.90sD-log(1+W_weight*t(b1)*exp(logeta.pred8.90sD)))
  
  
  etatilde_pred_m1 = apply(eta.tilde.pred1, 2, "median")
  z.etat1 <- rasterize(spred.fin, rr, log(etatilde_pred_m1,10), fun=mean)
  etatilde_pred_m8_90sA=apply(eta.tilde.pred8.90sA, 2, "median")
  etatilde_pred_m8_90sB=apply(eta.tilde.pred8.90sB, 2, "median")
  etatilde_pred_m8_90sD=apply(eta.tilde.pred8.90sD, 2, "median")
  z.etat8.90sA.absch <- rasterize(spred.fin, rr, etatilde_pred_m8_90sA - etatilde_pred_m1, fun=mean)
  z.etat8.90sB.absch <- rasterize(spred.fin, rr, etatilde_pred_m8_90sB - etatilde_pred_m1, fun=mean)
  z.etat8.90sD.absch <- rasterize(spred.fin, rr, etatilde_pred_m8_90sD - etatilde_pred_m1, fun=mean)
  
  
  # Plots only ###
  invisible(gc())
  # New maps
  
  
  # RCP 8.5 90s Clim Mod A,B,D ## (Supplement)
  
  
  # Suitability abs.change
  par(mfrow=c(1,4),mar= c(1, 4, 3, 3), oma = c(1,1,.1,.1), bty="L", mgp = c(2,1,0)) # , mar = c(2,5,1,.1)
  
  # color scale fixed 
  min.ztheta = quantile(c(values((z.t8.90sA.absch)),values((z.t8.90sB.absch)),values((z.t8.90sD.absch))),.0001, na.rm=T)
  max.ztheta = max(c(values((z.t8.90sA.absch)),values((z.t8.90sB.absch)),values((z.t8.90sD.absch))), na.rm=T)
  lim.z = max(max.ztheta, abs(min.ztheta))
  breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot(z.t8.90sA.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="Climate model A", cex.main = 1.5,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F); plot.coords()
  title(ylab = "Sites suitability", line = 2.5, cex.lab= 1.5)
  plot(z.t8.90sB.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="Climate model B", cex.main = 1.5,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F); plot.coords()
  plot(z.t8.90sD.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main="Climate model D", cex.main = 1.5,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), col=coldiv, breaks = breaks, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  
  
  # Relative variation in proliferation rate abs.change
  par(mfrow=c(1,4),mar= c(1, 4, 3, 3), oma = c(1,1,.1,.1), bty="L", mgp = c(2,1,0)) # , mar = c(2,5,1,.1)
  
  # color scale fixed 
  min.za = quantile(c(values((z.a8.90sA.absch)),values((z.a8.90sB.absch)),values((z.a8.90sD.absch))),.0001, na.rm=T)
  max.za = max(c(values((z.a8.90sA.absch)),values((z.a8.90sB.absch)),values((z.a8.90sD.absch))), na.rm=T)
  lim.z = max(max.za, abs(min.za))
  breaks <- seq(-lim.z,lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  
  plot(z.a8.90sA.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  title(ylab = "Relative variation in proliferation rate", line = 2.5, cex.lab= 1.5)
  plot(z.a8.90sB.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.a8.90sD.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)),  
                      cex.axis=1))
  
  
  # Larvae density (log(eta.tilde) )
  par(mfrow=c(1,4),mar= c(1, 4, 3, 3), oma = c(1,1,.1,.1), bty="L", mgp = c(2,1,0)) # , mar = c(2,5,1,.1)
  
  # color scale fixed 
  min.zetat = quantile(c(values((z.etat8.90sA.absch)),values((z.etat8.90sB.absch)),values((z.etat8.90sD.absch))),.005, na.rm=T)
  max.zetat = quantile(c(values((z.etat8.90sA.absch)),values((z.etat8.90sB.absch)),values((z.etat8.90sD.absch))),1, na.rm=T) # q. changed to 1 on 151225
  lim.z = max(max.zetat, abs(min.zetat))
  breaks <- seq(-lim.z, lim.z, length.out = length(coldiv)+1)
  z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
  plot(z.etat8.90sA.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  title(ylab = "Larvae density", line = 2.5, cex.lab= 1.5)
  plot(z.etat8.90sB.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.etat8.90sD.absch, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T,
       ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])), breaks= breaks,col=coldiv, legend=F); plot.coords()
  plot(z.brk, legend.only=TRUE, col=coldiv,
       legend.width = 2,breaks= breaks,legend.args = list(text="Abs. change", 3, cex=1),
       axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                      labels=pretty.round(min(breaks), max(breaks)), 
                      cex.axis=1))
  
  invisible(gc())
  par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))
}


# Fig. 2: Variance partition #################
# img.ratio: 

# Varpart part barplots
covnames1 = c("Distance to sand", "River influence" , "Distance to deep", "Exposure","Chlorophyll-a", "Winter ice","Salinity" )

for(f in 1:2){
  fit.time.quad = fit.list1[[f]]
  ## load post param
  # linear coeff posteriors
  alphabar = as.matrix(fit.time.quad, "alpha_bar")
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  abar = as.matrix(fit.time.quad, "a_bar")
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q]

  delta = as.matrix(fit.time.quad, "delta")
  eta.pred1 = as.matrix(fit.time.quad, "eta_new")
  theta.pred1 = as.matrix(fit.time.quad, "theta_new")
  # add delta+alpha to eta
  eta.pred2008= exp(delta[,1:18]%*%W_new) * eta.pred1
  eta.pred2009= exp(delta[,19:36]%*%W_new) * eta.pred1
  eta.pred2010= exp(delta[,37:54]%*%W_new) * eta.pred1
  # average over time
  eta.pred1.m = (eta.pred2008+ eta.pred2009+eta.pred2010)/3
  rm(list=c("eta.pred1", "eta.pred2008", "eta.pred2009", "eta.pred2010"))
  a1 = exp(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(as.matrix(abar)) +
             xpred.fin[,cov3]%*%t(bbar1) +  xpred.fin[,cov3q]^2%*%t(bbar2)) # [n.loc x n.samp]
  
  
  # Posterior values - future 50s RCP 4.5 ###
  
  # Suitability (theta)
  ## 50s
  logit.theta.pred4.50s = t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                xpred.fut.fin4.50s[,cov1]%*%t(betabar1) +  xpred.fut.fin4.50s[,cov1q]^2%*%t(betabar2)) )

  # Proliferation rate (log(a) )
  ## 50s
  loga.pred4.50s = t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                         xpred.fut.fin4.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(bbar2)) )

  # spawners density (log(eta) )
  # avg delta over 3 years, for the region to which each cell belongs
  deltaalpha= t(2.475 + log((exp(delta[,1:18]%*%W_new) + exp(delta[,19:36]%*%W_new)+ exp(delta[,37:54]%*%W_new))/3) ) # [ n.loc x n.samples ] 
  ## 50s
  # constant change rate over each grid-cell
  eta.pred4 = exp(t( deltaalpha+ xpred.fut.fin4.50s[,cov3]%*%t(beta1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(beta2)) )
  c4 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred4, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred4.50s = matrix(log(c4), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred4)) + log(eta.pred4)

  invisible(gc())
  
  ## Variance partition current scenarios
  vpar.eta = varpart1(xpred.fin[,2:8], beta = beta1, beta2 = beta2, log(eta.pred1.m)); invisible(gc())
  vpar.theta = varpart1(xpred.fin[,2:8], beta = betabar1, beta2 = betabar2,logit(theta.pred1)); invisible(gc())
  vpar.a = varpart1(xpred.fin[,2:8], beta = bbar1, beta2 = bbar2,log(t(a1))); invisible(gc())
  
  var.eta.mean = mean(apply(log(eta.pred1.m),1, var) )
  var.theta.mean = mean(apply(logit(theta.pred1),1, var) , na.rm = T)
  var.a.mean = mean(apply(log(t(a1)),1, var) )
  tot.vars = c( var.theta.mean,var.eta.mean,var.a.mean); invisible(gc())
  
  ## 50s
  ## Variance partition RCP 4.5
  vpar.eta.4.50s = varpart1(xpred.fut.fin4.50s[,2:8], beta = beta1, beta2 = beta2, logeta.pred4.50s); invisible(gc())
  vpar.theta.4.50s = varpart1(xpred.fut.fin4.50s[,2:8], beta = betabar1, beta2 = betabar2,(logit.theta.pred4.50s)); invisible(gc())
  vpar.a.4.50s = varpart1(xpred.fut.fin4.50s[,2:8], beta = bbar1, beta2 = bbar2,loga.pred4.50s); invisible(gc())
  
  var.eta.mean.4.50s = mean(apply(logeta.pred4.50s,1, var) )
  var.theta.mean.4.50s = mean(apply((logit.theta.pred4.50s),1, var) , na.rm = T)
  var.a.mean.4.50s = mean(apply(loga.pred4.50s,1, var) )
  tot.vars.4.50s = c( var.theta.mean.4.50s,var.eta.mean.4.50s,var.a.mean.4.50s); invisible(gc())
  
  # save average values in one matrix [6 x 7] to plot them in order
  # Mean values
  vpar.mtr = rbind(apply(vpar.theta, 2, mean, na.rm=T), apply(vpar.theta.4.50s, 2, mean, na.rm=T), # suit
                   apply(vpar.eta, 2, mean, na.rm=T), apply(vpar.eta.4.50s, 2, mean, na.rm=T),     # spawner
                   apply(vpar.a, 2, mean, na.rm=T), apply(vpar.a.4.50s, 2, mean, na.rm=T))         # prolif.rate
  # CI 50%
  q1=0.1;q2=0.9
  upper.v = rbind(apply(vpar.theta, 2, quantile,q2, na.rm=T), apply(vpar.theta.4.50s, 2, quantile,q2, na.rm=T), # suit
                  apply(vpar.eta, 2, quantile,q2, na.rm=T), apply(vpar.eta.4.50s, 2, quantile,q2, na.rm=T),     # spawner
                  apply(vpar.a, 2, quantile,q2, na.rm=T), apply(vpar.a.4.50s, 2, quantile,q2, na.rm=T))         # prolif.rate
  lower.v = rbind(apply(vpar.theta, 2, quantile,q1, na.rm=T), apply(vpar.theta.4.50s, 2, quantile,q1, na.rm=T), # suit
                  apply(vpar.eta, 2, quantile,q1, na.rm=T), apply(vpar.eta.4.50s, 2, quantile,q1, na.rm=T),     # spawner
                  apply(vpar.a, 2, quantile,q1, na.rm=T), apply(vpar.a.4.50s, 2, quantile,q1, na.rm=T))         # prolif.rate
  
  rm(list = c("vpar.theta", "vpar.eta", "vpar.a","vpar.theta.4.50s", "vpar.eta.4.50s", "vpar.a.4.50s"))
  
  
  # barplots ##
  covnames1 = c("Distance to sand", "River influence" , "Distance to deep", "Exposure","Chlorophyll-a", "Winter ice","Salinity" )
  mycex = .8
  mycol = c(rgb(1,.1,.2),rgb(1,.1,.2, .3), rgb(.2,.8,.3), rgb(.2,.8,.3, .3), rgb(.3,.1,1), rgb(.3,.1,1, .3))
  par(mfrow = c(1, 1), mar = c(2, 2, 2, 1))
  range.all = c(0,.8) #range(vpar.mtr,  na.rm = T) 
  plot("", xlim = c(0.5, 7+0.5), ylim = range.all, axes= F, xaxt='n', bty = 'n', ylab="", xlab="", cex.lab=.9)
  # Add barplot from matrix
  brp=barplot(vpar.mtr,beside=T, ylab ="",xlab = "", cex.axis =mycex, yaxt="n", col=mycol,border=NA,las = 1, cex.names =mycex, names.arg = covnames1, ylim =range.all) #
  arrows(brp,upper.v, brp, lower.v, angle=90, code=3, length=.05) # col="gray, lty=3
  legend("topleft", legend = paste0(c("suitability (","spawner density (", "max proliferation rate (" ),round(tot.vars,2),")"), fill=mycol[c(1,3,5)], bty="n", cex =mycex, title = "2010s (total variance)")
  legend("top", legend = paste0(c("suitability (","spawner density (", "max proliferation rate (" ),round(tot.vars.4.50s,2),")"), fill=mycol[c(2,4,6)], bty="n", cex =mycex, title = "RCP 4.5, 2050s (total variance)")
  axis(2, at =round(seq(range.all[1],(range.all[2]), l=7),2), cex.axis=mycex, las =2, line =-1.5, ylab = "")
  # title(ylab ="Proportion of variance", line =1, cex.lab=mycex)
  title(main = names(fit.list1)[f], line = 3)
  par(mfrow = c(1, 1), mar = c(5, 4, 4, 2)+.1)
  
}

# Fig S5: Posterior of abs. change ############
# img.ratio: 815 x 425

subreg.list = list(c(27,28,32, 37,42,47), c(23,24), c(19,20,15,16,11,12,6,7,2,3))
names(subreg.list) = c("Bothnian Sea", "Quark", "Bothnian Bay")
subreg.nbrsites = c(sum(region.fin%in% subreg.list[[1]]), sum(region.fin%in% subreg.list[[2]]), sum(region.fin%in% subreg.list[[3]])) # wrong length: must be nbr MC samples

for(f in 1:2){
  # load posteriors
  fit.time.quad = fit.list1[[f]]
  ## load post param
  # linear coeff posteriors
  alphabar = as.matrix(fit.time.quad, "alpha_bar")
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  abar = as.matrix(fit.time.quad, "a_bar")
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q]

  delta = as.matrix(fit.time.quad, "delta")
  eta.pred1 = as.matrix(fit.time.quad, "eta_new")
  theta.pred1 = as.matrix(fit.time.quad, "theta_new")
  # add delta+alpha to eta
  eta.pred2008= exp(delta[,1:18]%*%W_new) * eta.pred1
  eta.pred2009= exp(delta[,19:36]%*%W_new) * eta.pred1
  eta.pred2010= exp(delta[,37:54]%*%W_new) * eta.pred1
  # average over time
  eta.pred1.m = (eta.pred2008+ eta.pred2009+eta.pred2010)/3
  rm(list=c("eta.pred1", "eta.pred2008", "eta.pred2009", "eta.pred2010"))
  a1 = exp(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(as.matrix(abar)) +
             xpred.fin[,cov3]%*%t(bbar1) +  xpred.fin[,cov3q]^2%*%t(bbar2)) 
  b1= matrix(rep(as.matrix(fit.time.quad, "b"), nrow(xpred.fin)), nrow=nrow(xpred.fin), byrow=T) 
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new) # [N.samp x n.loc]
  eta.tilde.pred1 = W_weight*exp(log(t(a1)) + log(eta.pred1.m)-log(1+W_weight*t(b1)*eta.pred1.m))
  
  # Median values - current, future 50s & 90s + RCP 4.5 & 8.5 ###
  
  # Suitability (theta)
  ## 50s
  theta.pred4.50s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin4.50s[,cov1]%*%t(betabar1) +  xpred.fut.fin4.50s[,cov1q]^2%*%t(betabar2)) )
  theta.pred8.50s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin8.50s[,cov1]%*%t(betabar1) +  xpred.fut.fin8.50s[,cov1q]^2%*%t(betabar2)) )
  
  ## 90s
  theta.pred4.90s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin4.90s[,cov1]%*%t(betabar1) +  xpred.fut.fin4.90s[,cov1q]^2%*%t(betabar2)) )
  theta.pred8.90s = t( inv.logit(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(alphabar) +
                                   xpred.fut.fin8.90s[,cov1]%*%t(betabar1) +  xpred.fut.fin8.90s[,cov1q]^2%*%t(betabar2)) )
  
  
  # Proliferation rate (log(a) )
  ## 50s
  loga.pred4.50s = t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                         xpred.fut.fin4.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(bbar2)) )
  loga.pred8.50s = t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                         xpred.fut.fin8.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin8.50s[,cov3q]^2%*%t(bbar2)) )
  ## 90s
  loga.pred4.90s = t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                         xpred.fut.fin4.90s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.90s[,cov3q]^2%*%t(bbar2)) )
  loga.pred8.90s = t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                         xpred.fut.fin8.90s[,cov3]%*%t(bbar1) +  xpred.fut.fin8.90s[,cov3q]^2%*%t(bbar2)) )
  
  
  
  # spawners density (log(eta) )
  # avg delta over 3 years, for the region to which each cell belongs
  deltaalpha= t(2.475 + log((exp(delta[,1:18]%*%W_new) + exp(delta[,19:36]%*%W_new)+ exp(delta[,37:54]%*%W_new))/3) ) # [ n.loc x n.samples ] 
  ## 50s
  # constant change rate over each grid-cell
  eta.pred4 = exp(t( deltaalpha+ xpred.fut.fin4.50s[,cov3]%*%t(beta1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(beta2)) )
  eta.pred8 = exp(t( deltaalpha+ xpred.fut.fin8.50s[,cov3]%*%t(beta1) +  xpred.fut.fin8.50s[,cov3q]^2%*%t(beta2)) )
  c4 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred4, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  c8 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred4.50s = matrix(log(c4), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred4)) + log(eta.pred4)
  logeta.pred8.50s = matrix(log(c8), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8)) + log(eta.pred8)
  ## 90s
  # constant change rate over each grid-cell
  eta.pred4 = exp(t( deltaalpha+ xpred.fut.fin4.90s[,cov3]%*%t(beta1) +  xpred.fut.fin4.90s[,cov3q]^2%*%t(beta2)) )
  eta.pred8 = exp(t( deltaalpha+ xpred.fut.fin8.90s[,cov3]%*%t(beta1) +  xpred.fut.fin8.90s[,cov3q]^2%*%t(beta2)) )
  c4 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred4, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  c8 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred8, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred4.90s = matrix(log(c4), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred4)) + log(eta.pred4)
  logeta.pred8.90s = matrix(log(c8), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred8)) + log(eta.pred8)
  
  
  # larvae density eta.tilde
  # from density to number? *90000m2
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new)
  ## 50s
  eta.tilde.pred4.50s =W_weight*exp(loga.pred4.50s+logeta.pred4.50s-log(1+W_weight*t(b1)*exp(logeta.pred4.50s)))
  eta.tilde.pred8.50s =W_weight*exp(loga.pred8.50s+logeta.pred8.50s-log(1+W_weight*t(b1)*exp(logeta.pred8.50s)))
  ## 90s 
  eta.tilde.pred4.90s =W_weight*exp(loga.pred4.90s+logeta.pred4.90s-log(1+W_weight*t(b1)*exp(logeta.pred4.90s)))
  eta.tilde.pred8.90s =W_weight*exp(loga.pred8.90s+logeta.pred8.90s-log(1+W_weight*t(b1)*exp(logeta.pred8.90s)))
  
  invisible(gc())
  
  # Suitability index
  suit.pred1 = matrix(rbinom(length(theta.pred1), size = 1,theta.pred1 ) , ncol = ncol(theta.pred1))
  suit.pred4.50s = matrix(rbinom(length(theta.pred4.50s), size = 1,theta.pred4.50s ) , ncol = ncol(theta.pred4.50s))
  suit.pred8.50s = matrix(rbinom(length(theta.pred8.50s), size = 1,theta.pred8.50s ) , ncol = ncol(theta.pred8.50s))
  suit.pred4.90s = matrix(rbinom(length(theta.pred4.90s), size = 1,theta.pred4.90s ) , ncol = ncol(theta.pred4.90s))
  suit.pred8.90s = matrix(rbinom(length(theta.pred8.90s), size = 1,theta.pred8.90s ) , ncol = ncol(theta.pred8.90s))
  
  
  # Vioplots quantities
  proc.list1 = list(suit.pred1, eta.pred1.m, eta.tilde.pred1)
  proc.list4.50s = list(suit.pred4.50s,  eta.tilde.pred4.50s)
  proc.list8.50s = list(suit.pred8.50s, eta.tilde.pred8.50s)
  proc.list4.90s = list(suit.pred4.90s, eta.tilde.pred4.90s)
  proc.list8.90s = list(suit.pred8.90s,eta.tilde.pred8.90s)
  
  names(proc.list4.50s) = names(proc.list8.50s) = names(proc.list4.90s) = names(proc.list8.90s) =  c("absolute change\n in fraction of suitable sites","log relative change\n in total larvae production") 
  invisible(gc())
  
  # Vioplots
  mycol = rep(blues9[4],3) 
  mycol1 = rep(blues9[8],3)  

  par(mfrow =c(1,2), mar =c(4,2,1,1), oma = c(1,5,1,1), mgp = c(3,1,0))
  
  # 50s : suitable sites
  for (p1 in 1 ) {
    # for each process = suitable site
    #proces1 = proc.list1[[p1]]
    proces4 = proc.list4.50s[[p1]] - proc.list1[[p1]]
    proces8 = proc.list8.50s[[p1]] - proc.list1[[p1]]
    # RCP 4.5 
    # posterior distribution of absolute change of p1... per region
    l4=list((apply(proces4[,which(region.fin%in% subreg.list[[1]])],1, sum))/subreg.nbrsites[1], 
            (apply(proces4[,which(region.fin%in% subreg.list[[2]])],1, sum))/subreg.nbrsites[2],  
            (apply(proces4[,which(region.fin%in% subreg.list[[3]])],1, sum))/subreg.nbrsites[3])
    # RCP 8.5 
    # posterior distribution of absolute change of p1... per region
    l8=list((apply(proces8[,which(region.fin%in% subreg.list[[1]])],1, sum))/subreg.nbrsites[1], 
            (apply(proces8[,which(region.fin%in% subreg.list[[2]])],1, sum))/subreg.nbrsites[2],  
            (apply(proces8[,which(region.fin%in% subreg.list[[3]])],1, sum))/subreg.nbrsites[3])
    
    
    # Set up plot without violins
    plot("", ylim = c(0.5, length(l4)+0.5), xlim = range.q.ll(list(l4,l8)),  yaxt='n', bty = 'n', ylab="", xlab=names(proc.list4.50s)[p1], cex.lab=.9)
    # Add violins from list
    invisible(lapply(seq_along(l4), function(xx)
      vioplot(l4[[xx]], at = xx, col = mycol[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "left")))
    invisible(lapply(seq_along(l8), function(xx)
      vioplot(l8[[xx]], at = xx, col = mycol1[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "right")))
    
    if(p1==1) axis(2, at =1:3, labels = names(subreg.list), cex.axis =1, las =2) #else mtext(expression("10"^x), side = 1, line = 4, cex = .7)
    if(p1==1) title(main = "2050s") # Title
  }
  
 
  # 50s : larvae (log(10))
  for (p1 in 2 ) {
    # for each process =  tot larvae
    proces1 = proc.list1[[p1]]
    proces4 = proc.list4.50s[[p1]] 
    proces8 = proc.list8.50s[[p1]] 
    # RCP 4.5 
    # posterior distribution of absolute change of p1... per region
    l4=list(log(apply(proces4[,which(region.fin%in% subreg.list[[1]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[1]])],1, sum), 10), 
            log(apply(proces4[,which(region.fin%in% subreg.list[[2]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[2]])],1, sum), 10),  
            log(apply(proces4[,which(region.fin%in% subreg.list[[3]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[3]])],1, sum), 10))
    # RCP 8.5 
    # posterior distribution of log relative change of p1... per region
    l8=list(log(apply(proces8[,which(region.fin%in% subreg.list[[1]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[1]])],1, sum), 10), 
            log(apply(proces8[,which(region.fin%in% subreg.list[[2]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[2]])],1, sum), 10),  
            log(apply(proces8[,which(region.fin%in% subreg.list[[3]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[3]])],1, sum), 10))
    # From Number to Fraction of suitable sites
    if(p1==1){
      subreg.nbrsites1 = c(rep(subreg.nbrsites[1], length(l4[[1]])), rep(subreg.nbrsites[2], length(l4[[1]])), rep(subreg.nbrsites[3], length(l4[[1]])))
      l4 = lapply(l4, '/', subreg.nbrsites1)
      l8 = lapply(l8, '/', subreg.nbrsites1)      }
    
    # Set up plot without violins
    plot("", ylim = c(0.5, length(l4)+0.5), xlim = range.q.ll(list(l4,l8)),  yaxt='n', bty = 'n', ylab="", xlab=names(proc.list4.50s)[p1], cex.lab=.9)
    # Add violins from list
    invisible(lapply(seq_along(l4), function(xx)
      vioplot(l4[[xx]], at = xx, col = mycol[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "left")))
    invisible(lapply(seq_along(l8), function(xx)
      vioplot(l8[[xx]], at = xx, col = mycol1[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "right")))
    
    if(p1==1) axis(2, at =1:3, labels = names(subreg.list), cex.axis =1, las =2) #else mtext(expression("10"^x), side = 1, line = 4, cex = .7)
    if(p1==1) title(main = "2050s") # Title
  }
  
  # 90s : suitable sites
  for (p1 in 1 ) {
    # for each process = suitable site 
    #proces1 = proc.list1[[p1]]
    proces4 = proc.list4.90s[[p1]] - proc.list1[[p1]]
    proces8 = proc.list8.90s[[p1]] - proc.list1[[p1]]
    # RCP 4.5 
    # posterior distribution of absolute change of p1... per region
    l4=list((apply(proces4[,which(region.fin%in% subreg.list[[1]])],1, sum))/subreg.nbrsites[1], 
            (apply(proces4[,which(region.fin%in% subreg.list[[2]])],1, sum))/subreg.nbrsites[2],  
            (apply(proces4[,which(region.fin%in% subreg.list[[3]])],1, sum))/subreg.nbrsites[3])
    # RCP 8.5 
    # posterior distribution of absolute change of p1... per region
    l8=list((apply(proces8[,which(region.fin%in% subreg.list[[1]])],1, sum))/subreg.nbrsites[1], 
            (apply(proces8[,which(region.fin%in% subreg.list[[2]])],1, sum))/subreg.nbrsites[2],  
            (apply(proces8[,which(region.fin%in% subreg.list[[3]])],1, sum))/subreg.nbrsites[3])
    
    
    # Set up plot without violins
    plot("", ylim = c(0.5, length(l4)+0.5), xlim = range.q.ll(list(l4,l8)),  yaxt='n', bty = 'n', ylab="", xlab=names(proc.list4.50s)[p1], cex.lab=.9)
    # Add violins from list
    invisible(lapply(seq_along(l4), function(xx)
      vioplot(l4[[xx]], at = xx, col = mycol[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "left")))
    invisible(lapply(seq_along(l8), function(xx)
      vioplot(l8[[xx]], at = xx, col = mycol1[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "right")))
    
    if(p1==1) axis(2, at =1:3, labels = names(subreg.list), cex.axis =1, las =2) #else mtext(expression("10"^x), side = 1, line = 4, cex = .7)
    if(p1==1) title(main = "2090s") # Title
    invisible(gc())
  }
  
  # 90s : Larvae (log10)
  for (p1 in 2 ) {
    # for each process = tot larvae
    proces1 = proc.list1[[p1]]
    proces4 = proc.list4.90s[[p1]] 
    proces8 = proc.list8.90s[[p1]] 
    # RCP 4.5 
    # posterior distribution of log relative change of p1... per region
    l4=list(log(apply(proces4[,which(region.fin%in% subreg.list[[1]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[1]])],1, sum), 10), 
            log(apply(proces4[,which(region.fin%in% subreg.list[[2]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[2]])],1, sum), 10),  
            log(apply(proces4[,which(region.fin%in% subreg.list[[3]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[3]])],1, sum), 10))
    # RCP 8.5 
    # posterior distribution of absolute change of p1... per region
    l8=list(log(apply(proces8[,which(region.fin%in% subreg.list[[1]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[1]])],1, sum), 10), 
            log(apply(proces8[,which(region.fin%in% subreg.list[[2]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[2]])],1, sum), 10),  
            log(apply(proces8[,which(region.fin%in% subreg.list[[3]])],1, sum), 10)- log(apply(proces1[,which(region.fin%in% subreg.list[[3]])],1, sum), 10))

    
    # Set up plot without violins
    plot("", ylim = c(0.5, length(l4)+0.5), xlim = range.q.ll(list(l4,l8)),  yaxt='n', bty = 'n', ylab="", xlab=names(proc.list4.50s)[p1], cex.lab=.9)
    # Add violins from list
    invisible(lapply(seq_along(l4), function(xx)
      vioplot(l4[[xx]], at = xx, col = mycol[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "left")))
    invisible(lapply(seq_along(l8), function(xx)
      vioplot(l8[[xx]], at = xx, col = mycol1[xx], add = T, axes=F, bty = 'n',horizontal=T, side = "right")))
    
    if(p1==1) axis(2, at =1:3, labels = names(subreg.list), cex.axis =1, las =2) #else mtext(expression("10"^x), side = 1, line = 4, cex = .7)
    if(p1==1) title(main = "2090s") # Title
  }
}
par(mfrow=c(1,1), mar= c(5,4,4,2)+.1)
plot.new()
legend("center", c("RCP 4.5", "RCP 8.5"), fill = c(mycol[1], mycol1[1])) # RCPs legend



# Fig. 7: Abs. change in impact on larvae ###########
# img.ratio: 771 x 821

prc.names = c("Proliferarion rate","Spawner density")
prc.names1 = c("Larvae rate of change \n wrt proliferarion rate","Larvae rate of change \n wrt  spawner density")
title.pos = c(0,-32)
par(mfrow=c(1,4),mar= c(1, 1, 3.1, 6), oma = c(1,1,.1,.1), bty="L")

for(f in 1:2){
  # load posteriors
  fit.time.quad = fit.list1[[f]]
  ## load post param
  # linear coeff posteriors
  alphabar = as.matrix(fit.time.quad, "alpha_bar")
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  abar = as.matrix(fit.time.quad, "a_bar")
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q]

  delta = as.matrix(fit.time.quad, "delta")
  eta.pred1 = as.matrix(fit.time.quad, "eta_new")
  # add delta+alpha to eta
  eta.pred2008= exp(delta[,1:18]%*%W_new) * eta.pred1
  eta.pred2009= exp(delta[,19:36]%*%W_new) * eta.pred1
  eta.pred2010= exp(delta[,37:54]%*%W_new) * eta.pred1
  # average over time
  eta.pred1.m = (eta.pred2008+ eta.pred2009+eta.pred2010)/3
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new) # [N.samp x n.loc]
  a1 = t(W_weight)*exp(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(as.matrix(abar)) +
                         xpred.fin[,cov3]%*%t(bbar1) +  xpred.fin[,cov3q]^2%*%t(bbar2)) 
  b1= matrix(rep(as.matrix(fit.time.quad, "b"), nrow(xpred.fin)), nrow=nrow(xpred.fin), byrow=T) 
  
  # Median values - current, future 50s + RCP 4.5###
  
  # Proliferation rate (log(a) )
  ## 50s
  loga.pred4.50s =  log(W_weight) + t( (matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(abar) +
                                          xpred.fut.fin4.50s[,cov3]%*%t(bbar1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(bbar2)) )
  
  
  # spawners density (log(eta) )
  # avg delta over 3 years, for the region to which each cell belongs
  deltaalpha= t(2.475 + log((exp(delta[,1:18]%*%W_new) + exp(delta[,19:36]%*%W_new)+ exp(delta[,37:54]%*%W_new))/3) ) # [ n.loc x n.samples ] 
  ## 50s
  # constant change rate over each grid-cell
  eta.pred4 = exp(t( deltaalpha+ xpred.fut.fin4.50s[,cov3]%*%t(beta1) +  xpred.fut.fin4.50s[,cov3q]^2%*%t(beta2)) )
  c4 = apply(eta.pred1.m, 1, sum)/ apply(eta.pred4, 1, sum) #[s x 1] tot.spawn.biomass change coeff
  logeta.pred4.50s = matrix(log(c4), ncol = 1)%*%matrix(1, nrow = 1,ncol = ncol(eta.pred4)) + log(eta.pred4)
  
  # Impacts
  b0 = as.matrix(fit.time.quad, "b")
  Imp0 = impact.larvae(t(a1), eta.pred1.m,b0 )
  Imp.cl4.50s = impact.larvae(exp(loga.pred4.50s), exp(logeta.pred4.50s),b0 ) 
  invisible(gc())
  
  # Imp: plot of abs. change d(eta.tilde)/d(p0) ###
  # Raster maps
  # p0: prolif.rate & spawner dens
  for(p0 in 1:2){
    imp_cl_m0 = (apply(Imp0[[p0]], 2, "median"))
    #z.imp0 <- rasterize(spred.fin, rr, imp_cl_m0, fun=mean)
    imp_cl_m4_50s=(apply(Imp.cl4.50s[[p0]] , 2, "median"))
    z.imp4.50s <- rasterize(spred.fin, rr, imp_cl_m4_50s - imp_cl_m0, fun=mean)
    
    # color scale fixed
    max.zimp = max(abs(imp_cl_m4_50s - imp_cl_m0)) 
    breaks <- seq(-max.zimp, max.zimp, length.out = length(colcon)+1)
    z.brk <- rasterize(spred.fin, rr, breaks, fun=mean)
    
    # plot
    plot(z.imp4.50s, xlim=cbind(min(spred.fin[,1]),max(spred.fin[,1])), xaxt='n', yaxt='n',legend.width=2,axes = FALSE,box=T, main=prc.names1[p0],
         ylim=cbind(min(spred.fin[,2]),max(spred.fin[,2])),col=coldiv, legend=F, breaks= breaks) # 
    plot.coords()
    plot(z.brk, legend.only=TRUE, col=coldiv,
         legend.width = 2,breaks= breaks,legend.args = list(text="absoulte change", 3, cex = 1),
         axis.args=list(at=pretty.round(min(breaks), max(breaks)),
                        labels=pretty.round(min(breaks), max(breaks)),
                        cex.axis=1))
  }
  title(main = names(fit.list1)[f], outer = T, line =title.pos[f])
  invisible(gc())
}
par(mfrow=c(1,1), mar= c(5,4,4,2)+.1)



# Fig. 6: reproduction curves ###########
# img.ratio: 725 x 620

subreg.list = list(c(27,28,32, 37,42,47), c(23,24), c(19,20,15,16,11,12,6,7,2,3))
names(subreg.list) = c("Bothnian Sea", "Quark", "Bothnian Bay")
cov1 = cov2 = cov3 = cov1q = cov2q = cov3q = 2:8

n.loc = 100 # location to sample, for each subreg

## predictive maps of Median log10(P)
for(f in 1:2){
  fit.time.quad = fit.list1[[f]]
  ## load post param
  # linear coeff posteriors
  alphabar = as.matrix(fit.time.quad, "alpha_bar")
  betabar1 = as.matrix(fit.time.quad,  paste0("beta_bar","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  betabar2 = as.matrix(fit.time.quad,  paste0("beta_bar_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # sites suitability
  beta1 = as.matrix(fit.time.quad,  paste0("beta","[",1:7,"]")); colnames(betabar1) = covnames[cov1]
  beta2 = as.matrix(fit.time.quad,  paste0("beta_q","[",1:7,"]")); colnames(betabar2) = covnames[cov1q] # spawner density
  abar = as.matrix(fit.time.quad, "a_bar")
  bbar1 = as.matrix(fit.time.quad,  paste0("b_bar","[",1:7,"]")); colnames(bbar1) = covnames[cov3]
  bbar2 = as.matrix(fit.time.quad,  paste0("b_bar_q","[",1:7,"]")); colnames(bbar2) = covnames[cov3q]

  delta = as.matrix(fit.time.quad, "delta")
  eta.pred1 = as.matrix(fit.time.quad, "eta_new")
  theta.pred1 = as.matrix(fit.time.quad, "theta_new")
  # add delta+alpha to eta
  eta.pred2008= exp(delta[,1:18]%*%W_new) * eta.pred1
  eta.pred2009= exp(delta[,19:36]%*%W_new) * eta.pred1
  eta.pred2010= exp(delta[,37:54]%*%W_new) * eta.pred1
  # average over time
  eta.pred1.m = (eta.pred2008+ eta.pred2009+eta.pred2010)/3
  rm(list=c("eta.pred1", "eta.pred2008", "eta.pred2009", "eta.pred2010"))
  a1 = exp(matrix(1,nrow=dim(xpred.fin)[1],ncol=1)%*%t(as.matrix(abar)) +
             xpred.fin[,cov3]%*%t(bbar1) +  xpred.fin[,cov3q]^2%*%t(bbar2)) 
  b1= matrix(rep(as.matrix(fit.time.quad, "b"), nrow(xpred.fin)), nrow=nrow(xpred.fin), byrow=T) 
  W_weight = (matrix(1,nrow=dim(eta.pred1.m)[1],ncol=1))%*%(weight.mtr[,f]%*%W_new) # [N.samp x n.loc]
  eta.tilde.pred1 = W_weight*t(a1)*eta.pred1.m/(1+W_weight*t(b1)*eta.pred1.m)
  
  b = as.matrix(fit.time.quad, "b")
  # avg delta over 3 years, for the region to which each cell belongs
  deltaalpha= t(2.475 + log((exp(delta[,1:18]%*%W_new) + exp(delta[,19:36]%*%W_new)+ exp(delta[,37:54]%*%W_new))/3) ) # [ n.loc x n.samples ] 
  
  
  # mean values
  theta_pred_m1=apply(theta.pred1, 2, "mean")
  eta_pred_m1=apply(eta.pred1.m, 2, "mean")
  eta_tilde_pred_m1=apply(eta.tilde.pred1, 2, "mean")
  
  
  ## Reproduction model
  set.seed(123)
  
  # * sample z for each location from Bernoulli distribution with E[\theta]
  z =rbinom(length(theta_pred_m1), size = 1, theta_pred_m1)
  
  # * select 100 locations: z = 1
  # * - for each selected location, visualize the expected spawners density, reproduction curve, and expected larvae density 
  for (reg in 1:3) {
    # * sample 100 locations from set of locations with z=1 , in the subregion
    suit.loc = which((z==1)&(region.fin%in% subreg.list[[reg]]))
    suit.loc.samp = sample(suit.loc, n.loc)
    if(f==1) maxeta = 0.05 else maxeta = 0.15 # quantile(eta_pred_m1[suit.loc],.997)
    eta.temp =seq(0,maxeta, l =1000)
    #eta.tilde.temp = matrix(nrow=length(eta.temp), ncol = n.loc) # [1000 eta values, 100 locations]
    eta.mc = eta_pred_m1[suit.loc.samp] 
    eta.tilde.mc1 = eta_tilde_pred_m1[suit.loc.samp]
    
    a1.temp = apply(W_weight*t(a1), 2, mean)[suit.loc.samp] # [100 locations] expected value
    b1.temp = apply(W_weight*t(b1), 2, mean)[suit.loc.samp] # [100 locations]
    eta.tilde.temp = matrix(eta.temp, ncol = 1)%*% matrix(a1.temp, nrow = 1)/(1+matrix(eta.temp, ncol = 1)%*% matrix(b1.temp, nrow = 1)) # [1000 eta values, 100 locations]
    eta.tilde.mc = eta.mc*a1.temp/(1+ eta.mc*b1.temp)
    
    maxetatilde = quantile(eta.tilde.temp, .95)

    
    ## Plot
    mycol =viridis(n.loc)
    
    
    # layout for marginal density violinplots
    # Layout
    layout(matrix(c(2, 3, 0, 1),
                  nrow = 2, ncol = 2,
                  byrow = TRUE),
           widths = c(1,6),
           heights  = c(6, 1), respect = TRUE)
    
    # Bottom vioplot
    par(mar = c(2, 4.1, 0, 0), bty = "n", mgp =c(1,1,0))
    vioplot(eta_pred_m1[suit.loc], xaxt='n', yaxt='n', horizontal = TRUE,ylim= c(0,maxeta),
            col = "white", xlab= "spawner density (1/m2)", colMed = mycol[50])
    # Left vioplot
    par(mar = c(5.1, 2, 0, 0), bty = "n")
    vioplot(eta_tilde_pred_m1[suit.loc], xaxt='n', yaxt='n', ylim =c(0, maxetatilde),
            col = "white", ylab= "larval density (1/m3)", colMed = mycol[50])
    
    # Top and right margin of the main plot
    par(mar = c(5.1, 4.1, 0, 0))
    # plot
    mycol =viridis(n.loc)
    
    plot(eta.temp, eta.temp*0, col="white", "l", xlab= "", ylab="",axes=F,
         main ="", xlim= c(0,maxeta), ylim =c(0, maxetatilde)) #quantile(eta.tilde.temp,.9)
    axis(1, at = pretty(c(0,maxeta), n = 6))
    axis(2, at = pretty(c(0, maxetatilde), n = 6))
    
    c0=1
    for (mc in order(eta.mc)) {
      lines(eta.temp, eta.tilde.temp[,mc], col=mycol[c0])
      points((eta.mc[mc]), eta.tilde.mc[mc], col=mycol[c0], pch =16, cex=1)
      # points on axes
      points((eta.mc[mc]), par("usr")[3], col=mycol[c0], pch ="l", cex=1.5, xpd=T)
      points(par("usr")[1], eta.tilde.mc[mc], col=mycol[c0], pch ="—", cex=1.5, xpd=T)
      
      c0=c0+1
    }
    points(sort(eta.mc), eta.tilde.mc[order(eta.mc)], col=mycol, pch =16, cex=1)
    title(names(subreg.list)[reg], outer = T, cex.main = 2)   
  }
  invisible(gc())
}

par(mfrow=c(1,1),mar= c(5.1, 4.1, 4.1, 2.1))



