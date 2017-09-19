library(adegenet)
library(adegenet)

# install the dev version of diveRsity
library(devtools)
install_github("parallel")
library(diveRsity)

#MPBres<-divMigrate(infile="MPBpregenepopFULL.txt",outfile = "RUN2",boots=5,plot_network=TRUE,stat="d",para=TRUE,filter_threshold = 0.2)
fuck<-readLines("CICGENEPOP.txt")
tail(gsub("VirMAL]([:digits:]+)",fuck[54900],replacement="VirMAL"))
source("divMigrate_JW.R")
MPBres_JW<-divMigrate_JW(infile="MPBpregenepopFULL.txt",outfile = NULL,boots=100,plot_network=TRUE,stat="d",para=F,filter_threshold = 0.055, coordinates=as.matrix(coords),VS=7)
divMigrate("CICGENEPOP.txt", outfile=NULL, boots=10, plot_network=TRUE)
# Which Fst measurement?
# What bias is interesting?
coords_meta<-read.csv("coords.txt",header=TRUE)
coords<-read.csv("coords.txt",header=TRUE)[6:7]
coordsSQ<-read.csv("projXY_sq_sites.csv",header=TRUE)
coordsrd<-round(coords,0)
coordsSQrd<-round(coordsSQ,0)

coords_deg<-coords_meta[,c(4,3)][coordsrd[,1]%in%coordsSQrd[,1],]
ldscp_sel<-coords_meta$as_charact[coordsrd[,1]%in%coordsSQrd[,1]]
ldscp_sel<-droplevels(ldscp_sel)

MPBres_JW$dRelMig
gen_con<-MPBres_JW$dRelMig[which(unique(MPBoriginal00[,6])%in%ldscp_sel),which(unique(MPBoriginal00[,6])%in%ldscp_sel)]
diag(gen_con)<-0
####################################################################
gup<-gen_con
glo<-gen_con
gup[lower.tri(gen_con)]<-t(gen_con)[lower.tri(gen_con)]
glo[upper.tri(gen_con)]<-t(gen_con)[upper.tri(gen_con)]
diag(gup)<-0
diag(glo)<-0


wi<-read.csv("Ugeo.txt_ones.txt_projXY_sq_sites.csv_asymmetric_WindCost.cdmatrix.csv",header=FALSE)
wup<-wi
wlo<-wi
wup[lower.tri(wi)]<-t(wi)[lower.tri(wi)]
wlo[upper.tri(wi)]<-t(wi)[upper.tri(wi)]
diag(wup)<-0
diag(wlo)<-0

source("MMRR.R");

geo.dist<-dist(coordsSQrd)
glowlo<-MMRR(glo,list(as.matrix(wlo),as.matrix(geo.dist)),nperm=1000)
gupwup<-MMRR(gup,list(as.matrix(wup),as.matrix(geo.dist)),nperm=1000)
glowup<-MMRR(glo,list(as.matrix(wup),as.matrix(geo.dist)),nperm=100000)
gupwlo<-MMRR(gup,list(as.matrix(wlo),as.matrix(geo.dist)),nperm=1000)
MMRR(glo,list(as.matrix(geo.dist)),nperm=1000)
MMRR(glo,list(as.matrix(wlo)),nperm=1000)
MMRR(gup,list(as.matrix(geo.dist)),nperm=1000)
MMRR(gup,list(as.matrix(wup)),nperm=1000)

### Wind
wi<-read.csv("Ugeo.txt_ones.txt_projXY_sq_sites.csv_asymmetric_WindCost.cdmatrix.csv",header=FALSE)
# Minimum distances
wi_min<-matrix(data=NA,32,32)
g_min<-matrix(data=NA,32,32)
for (i in 1:32){
  for (j in 1:32){
    wi_min[i,j]<-min(wi[i,j],wi[j,i])
    wi_min[j,i]<-min(wi[i,j],wi[j,i])
    g_min[i,j]<-min(gen_con[i,j],gen_con[j,i])
    g_min[j,i]<-min(gen_con[i,j],gen_con[j,i])
  }
}
# Maximum distance
wi_max<-matrix(data=NA,32,32)
g_max<-matrix(data=NA,32,32)
for (i in 1:32){
  for (j in 1:32){
    wi_max[i,j]<-max(wi[i,j],wi[j,i])
    wi_max[j,i]<-max(wi[i,j],wi[j,i])
    g_max[i,j]<-max(gen_con[i,j],gen_con[j,i])
    g_max[j,i]<-max(gen_con[i,j],gen_con[j,i])
  }
}

wi_max
wi_min

gw_min<-MMRR(g_min, list(wi_max, as.matrix(geo.dist)), nperm=100)
gw_max<-MMRR(g_max ,list(wi_min, as.matrix(geo.dist)), nperm=100)

ggeo_min<-MMRR(g_min, list(as.matrix(geo.dist)), nperm=100)
ggeo_min
ggeo_max<-MMRR(g_max, list(as.matrix(geo.dist)), nperm=100)
ggeo_max

#####################################################################
##########################################################################
#################################################################################
########################################################################################





wi_std<-(as.matrix(wi)-mean(as.matrix(wi)))/sd(as.matrix(wi))
geo_std<-(as.matrix(geo.dist)-mean(as.matrix(geo.dist)))/sd(as.matrix(geo.dist))
wup_std<-wi_std
wlo_std<-wi_std
wup_std[lower.tri(wi_std)]<-t(wi_std)[lower.tri(wi_std)]
wlo_std[upper.tri(wi_std)]<-t(wi_std)[upper.tri(wi_std)]
diag(wup_std)<-0
diag(wlo_std)<-0
gupwlo_std<-MMRR(gup,list(as.matrix(wlo_std),as.matrix(geo_std)),nperm=1000)



mantel.rtest(dist(glo),dist(wlo))
mantel.rtest(dist(gup),dist(wup))

mantel.rtest(dist(glo),dist(wup))
mantel.rtest(dist(gup),dist(wlo))

mantel.rtest(dist(gup),geo.dist)
mantel.rtest(dist(glo),geo.dist)

plot(mantel.rtest(dist(gup),dist(wup)))
plot(mantel.rtest(dist(glo),dist(wlo)))
plot(mantel.rtest(dist(glo),dist(wlo)))
plot(mantel.rtest(dist(glo),geo.dist))

mantel <- function(x, y, n.iter=999, stat=function(a,b) sum(a*b)) {
  permute <- function(z) {
    i <- sample.int(nrow(z), nrow(z))
    return (z[i, i])
  }
  sapply(1:n.iter, function(i) stat(x, permute(y)))
}

hist(mantel(gen_con,wi),xlim=c(1.8*10^7,3.2*10^7),col="pink",main="Permutation results (wind)",xlab="Statistic")
abline(v=sum(gen_con*wi),col="red",lwd=4)

hist(mantel(gen_con,as.matrix(geo.dist)),xlim=c(1.4*10^7,2.2*10^7),col="pink",main="Permutation results (geographic)",xlab="Statistic")
abline(v=sum(gen_con*as.matrix(geo.dist)),col="red",lwd=4)

wi_std<-(as.matrix(wi)-mean(as.matrix(wi)))/sd(as.matrix(wi))
geo_std<-(as.matrix(geo.dist)-mean(as.matrix(geo.dist)))/sd(as.matrix(geo.dist))

wi_std<-(as.matrix(wi)-mean(as.matrix(wi)))/sd(as.matrix(wi))
geo_std<-(as.matrix(geo.dist)-mean(as.matrix(geo.dist)))/sd(as.matrix(geo.dist))

hist(wi_std)
hist(geo_std)

subs<-wi_std-geo_std
diag(subs)<-0
hist(mantel(gen_con,subs), col="pink", main="Permutation results (sub)",xlab="Statistic", breaks=100, yaxs = "i",  xaxs = "i")
abline(v=sum(gen_con*subs),col="red",lwd=4)

div<-wi/as.matrix(geo.dist)
diag(div)<-0
hist(mantel(gen_con,div), col="pink", main="Permutation results (div)",xlab="Statistic", breaks=100, yaxs = "i",  xaxs = "i")
abline(v=sum(gen_con*div),col="red",lwd=4)

gup.gdm<-data.frame(ldscp_sel, as.matrix(gup), stringsAsFactors=FALSE)

names(gup.gdm)<-c("ldscp_sel",1:32)

wup.gdm<-data.frame(ldscp_sel,as.matrix(wup),stringsAsFactors=FALSE)
names(wup.gdm)<-c("ldscp_sel",1:32)

wlo.gdm<-data.frame(ldscp_sel,as.matrix(wlo),stringsAsFactors=FALSE)
names(wlo.gdm)<-c("ldscp_sel",1:32)

coords.gdm_form<-data.frame(ldscp_sel,coordsSQ,stringsAsFactors=FALSE)
names(coords.gdm_form)<-c("ldscp_sel","POINT_X","POINT_Y")

#IMPORTANT: drop unused factor levels
frstpr<-formatsitepair(bioData=gup.gdm,bioFormat=3,predData=coords.gdm_form,distPreds = list(wup.gdm=wup.gdm), XColumn="POINT_X", YColumn="POINT_Y",siteColumn="ldscp_sel")
frstpr[,7]<-frstpr[,8]+rnorm(496)
head(frstpr)
gdm_all<-gdm(frstpr,geo=TRUE)
gdm.varImp(frstpr,geo=TRUE)


