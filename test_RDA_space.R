library(vegan)
library(adespatial)
library(adegenet)
library(bigsnpr)
x <- readRDS("sim1a.rds")
x1<-snp_attach("back/simu1.rds")
G <- x1$genotypes[]

# Find genetic clusters
clus<-find.clusters(G, n.pca = 250, n.clust=3)
# DAPC
dapc<-dapc(G,  grp=clus$grp, n.pca=250)
#PLots
scatter.dapc(dapc)

saveRDS(clus$grp, file="gen_clusters.rds")

# Combine coordinates
spa<-cbind(x1$fam$x,x1$fam$y)
# Create Moran's I Eigenvector Maps
mem_pos<-dbmem(spa,MEM.autocor = "positive")
# RDA with all MEM
mem_null.rda<-rda(G,mem_pos)
(mem_null.R2a<-RsquareAdj(mem_null.rda)$adj.r.squared)
(mem_null.R2<-RsquareAdj(mem_null.rda)$r.squared)
# Forward selection
mem_null.fwd<-forward.sel(G,mem_pos, adjR2thresh = 5, alpha=0.05, nperm=99)
mem_null.sign<-sort(functional_group1.MEMpos.fwd$order)
mem_null.order<-as.data.frame(mem_pos[c(MEMpos.sign)])