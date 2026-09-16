# load library and data ####
library("rstan")
library("VGAM")
library(matrixStats)
library("boot")
library("raster")
library("vioplot")
library("viridis")
library(LaplacesDemon)
library(RColorBrewer)

setwd("~Data")
load("FishGoBData")


# Add lat/lon grid to maps ####
# From etrs89-LAEA (UTM) to long/lat (WGS84)
# default raster: crs= etrs89
e <- extent(as.matrix(spred.fin))
rr <- raster(e, ncol=length(unique(spred.fin[,1])), nrow=length(unique(spred.fin[,2])))

z0 = rasterize(spred.fin, rr, xpred.fin[,1], fun=mean)

# Add parallels and meridians (crs=lat/lon) to plot of raster (crs=etrs89)
# to run only after plotting the map of interest: plot(z.xxx,...); plot.coords()
plot.coords = function(z.salt = z0){
  library(raster)
  library(sp)
  crs(z.salt) = CRS("+init=epsg:3035") # CRS etrs89
  # Add lat/lon lines
  # wgs coord to etrs coord: 
  # Parallels (wgs): lat= (60,62,64), lon=(19,...,26)
  lat.wgs = c(60,62,64)
  lon.wgs = seq(19,26, length= 100)
  for (l0 in lat.wgs) {
    wgs.coord = cbind( lon.wgs, rep(l0, 100) )
    pts = SpatialPoints(wgs.coord, proj4string=CRS("+proj=longlat +datum=WGS84"))
    pts_projected = spTransform(pts, CRS(projection(z.salt)) )
    etrs.coord = coordinates(pts_projected)
    lines( etrs.coord , lty=3, col="gray30" )
  }
  # Meridians (wgs): lat= (60,...,66), lon = (19,20,21)
  lat.wgs = seq(60,66,length = 100)
  lon.wgs = c(22,23,24)
  for (l0 in lon.wgs) {
    wgs.coord = cbind(rep(l0, 100), lat.wgs)
    pts = SpatialPoints(wgs.coord, proj4string=CRS("+proj=longlat +datum=WGS84"))
    pts_projected = spTransform(pts, CRS(projection(z.salt)) )
    etrs.coord = coordinates(pts_projected)
    lines( etrs.coord , lty=3, col="gray30" )
  }
  # coords1 <- locator(n = 2) # find location of tikcs
  axis(1, at = c(4970003, 5025658), labels = paste0(c(22,23), "°E"), las=1, cex.axis = 1.3)
  axis(2, at = c(4361493, 4586535), labels = paste0(c(62,64), "°N"), las=3, cex.axis = 1.3)
  
}

# Vioplots abs.change Fig. 4, 5, S2, S3 ####
vioplot.absch = function(v_pred_m1,v_pred_50s,  v_pred_90s, leg = T, log10 =F, rcp= "RCP 4.5", leg.pos="topleft"){
  if(log10){
    v_pred_m1 = 10^v_pred_m1
    v_pred_50s = 10^v_pred_50s
    v_pred_90s = 10^v_pred_90s
  }
  library(vioplot)
  par(mgp = c(2,1,0), bty = "L") # , mar = c(2,5,1,.1)
  y.min = min(v_pred_90s - v_pred_m1)-.1
  y.max =  max(v_pred_50s - v_pred_m1)+.1
  plot(1,1, col = "white", xlim = c(0.5,1.5), ylim = c(y.min, y.max), xaxt= 'n', ylab = "Absolute change", xlab = "", cex.lab = 1.5)
  # 
  vioplot(v_pred_50s - v_pred_m1, side = "left", col = "grey", add = T)
  vioplot(v_pred_90s - v_pred_m1, side = "right", col = "grey30", add = T)
  axis(1, at = 1, rcp, cex.axis =1.5)
  if(leg) legend(leg.pos, fill = c("gray", "gray30"), c("2050s", "2090s"), bty = 'n', cex = 1.5)
}

# Variance partition, Fig.2 ####
# Return: tot.var, posterior for each var. component
varpart1 =function(xpred, beta, beta2, etapred, covnames = covnames1){
  # variance of latent var, at each iteration
  var.eta = apply(etapred,1, var) 
  var.par = matrix(ncol = ncol(beta), nrow = nrow(beta)) # [ n.samples x n.covar ]
  colnames(var.par) = covnames
  for (j0 in 1:ncol(beta)) {
    var.par[,j0] = apply(beta[,j0]%*%t(xpred[,j0]) +beta2[,j0]%*%t(xpred[,j0]^2),1, var) / var.eta
  }
  return(var.par)
}



# Impact larvae Fig.7 ####
impact.larvae = function(a.pred, eta.pred, b, W_weihght = W_weight){
  denom.tmp = 1 + eta.pred*(matrix(b, ncol=1)%*%matrix(1, nrow=1, ncol = ncol(eta.pred)))*W_weight
  # eta
  imp.eta = a.pred*W_weight/ denom.tmp^2
  # a
  imp.a = eta.pred/denom.tmp
  
  return(list("Imp.a"=imp.a, "Imp.eta"=imp.eta))
}

# Other functions ####
# Pretty rounding of sequence 
pretty.round = function(min.sq, max.sq, l.sq = 6){
  sq = seq(min.sq, max.sq, l = 100)
  sq.round = pretty(sq, l.sq)
  return(sq.round)
}
# Get quantiles from list of lists, to set vioplots xlim (Fig.S5) 
range.q.ll = function(ll, qnt = c(.00005,.99995)){
  minmax = matrix(nrow = length(ll), ncol = 2)
  for(i0 in 1:length(ll)) minmax[i0,] = range(lapply(ll[[i0]], quantile, qnt, na.rm=T))
  return(range(minmax))
}

# Get spawner density intercept (from random effect) 
getalpha = function(delta.p, Sigma.delta.p, sigma2.alpha = 4){
  alpha = matrix(ncol = 100, nrow = nrow(delta.p))
  
  for(jj in 1:nrow(delta.p)){
    Sigma.delta.p.inv = solve(matrix(Sigma.delta.p[jj,], nrow = ncol(delta.p) ))
    unit.v = matrix(1, ncol = ncol(delta.p))
    mean.alpha = sigma2.alpha *unit.v %*% (Sigma.delta.p.inv)%*%delta.p[jj,]
    var.alpha = sigma2.alpha - sigma2.alpha*unit.v %*% Sigma.delta.p.inv %*% t(unit.v) *sigma2.alpha
    var.alpha[which(var.alpha<0)]=0 
    alpha[jj,] =rnorm(100, mean.alpha, sqrt(var.alpha))
  }
  
  return(alpha)}


