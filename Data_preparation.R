#######################
# Get Data for results in paper: 
# "Effect of climate change on reproduction of sea-spawning coregonids in the Baltic Sea"
#######################

# Load data from files ####
setwd("~Data")

data.train = read.table(file = "data_training.txt", col.names = T, row.names = F)
data.fishery = read.table(file = "data_fishery.txt", col.names = T, row.names = F)
data.gob = read.table(file = "data_fin_gob.txt", col.names = T, row.names = F)

# Assign data to R objects ####

# training data
t = data.train$Year
V = data.train$Water.volume
y = data.train$Whitefih
y.v = data.train$Vendace
reg = data.train$ICES.Reg
s.fin = cbind(data.train$E_etrs89, data.train$N_etrs89)
x.nost = data.train[, 5:11]

# predictive data
xpred.nost = data.gob[, 1:7]
reg.gob = data.gob$ICES.Reg
spred.fin = cbind(data.gob$E_etrs89, data.gob$N_etrs89)
covnames1 = colnames(data.gob[, 1:7])

mxcont = apply(x.nost, 2, mean)
stdxcont = apply(x.nost, 2, sd)
x.fin = t( apply( t(apply(x.nost,1,'-',mxcont)),1,'/',stdxcont ) )
xpred.fin = t( apply( t(apply(xpred.nost,1,'-',mxcont)),1,'/',stdxcont ) )
x.fin = cbind(rep(1, length(y)),x.fin) # add intercept
xpred.fin = cbind(rep(1, length(y)), xpred.fin) # add intercept

# fishery data
reg.inedx = data.fishery$ICES.Reg
E3 = data.fishery[,2:4]; R3 = data.fishery[,5:7]
E3.ven.gil = data.fishery[,8:10]; R3.ven.gil = data.fishery[,11:13]
E3.ven.fyke = data.fishery[,14:16]; R3.ven.fyke = data.fishery[,17:19]
E3.ven.trawl = data.fishery[,20:22]; R3.ven.trawl = data.fishery[,23:25]

weight.j = data.fishery$Weight_wf
weight.v = data.fishery$weight_ve
s.reg = cbind(data.fishery$E_ETRS89, data.fishery$N_ETRS89)

# Female ratio
f.ratio = 0.5487805 # wf
f.ratio.v = 0.58 # ve

# Summation matrices that picks up the region and time specific grid ####

# W [54 x 319] :  regions*yaers x sampling sites 
# rows 1:18 = region j in year 2009, rows 19:36 = region j in year 2010,  rows 37:54 = region j in year 2011
W = matrix(0,nrow=length(reg.index)*3,ncol=length(reg))
count = 1
for (i1 in reg.index*3){
  if(t[count]==2009) W[count,reg==i1] = 1 
  if(t[count]==2010) W[count+18,reg==i1] = 1 
  if(t[count]==2011) W[count+36,reg==i1] = 1 
  count = count+1;
}
# Weggs [18 x 319]: regions x sampling sites
Weggs = matrix(0,nrow=length(reg.index),ncol=length(reg))
count = 1
for (i1 in reg.index){
  Weggs[count,reg==i1] = 1 
  count = count+1;
}
# W_new [18 x 97078] : regions x raster.cells
W_new = matrix(0,nrow=length(reg.index),ncol=length(region))
count = 1
for (i1 in reg.index){
  W_new[count,region==i1] = 1 
  count = count+1;
}



# Save all data as one R object ####
rm(list = c("data.train", "data.gob", "data.fishery"))
save(list = ls(), file = "FishGoBData")
  