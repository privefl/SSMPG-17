x <- readRDS("sim1a.rds")
library(vegan)
library(adespatial)
library(adegenet)
library(bigsnpr)
# Find genetic clusters
clus<-find.clusters(t(x$G), n.pca = 300, n.clust = 4)
# DAPC
dapc<-dapc(t(x$G),  grp=clus$grp, n.pca=299, n.da=3)
#PLots
scatter.dapc(dapc)

# Combine coordinates
spa<-cbind(x$x,x$y)
# Create Moran's I Eigenvector Maps
mem_pos<-dbmem(spa,MEM.autocor = "positive")
# RDA with all MEM
mem_null.rda<-rda(t(x$G),mem_pos)
(mem_null.R2a<-RsquareAdj(mem_null.rda)$adj.r.squared)
(mem_null.R2<-RsquareAdj(mem_null.rda)$r.squared)
# Forward selection
mem_null.fwd<-forward.sel(t(x$G),mem_pos, adjR2thresh = mem_null.R2, alpha=0.05, nperm=99)
mem_null.sign<-sort(functional_group1.MEMpos.fwd$order)
mem_null.order<-as.data.frame(mem_pos[c(MEMpos.sign)])