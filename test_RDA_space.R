x <- readRDS("sim1a.rds")
library(vegan)
library(adespatial)

mod<-rda(t(x$G)~x$x+x$y+x$envi)
x$G
x$pop

spa<-cbind(x$x,x$y)
mem_pos<-dbmem(spa,MEM.autocor = "positive")

modsumm<-summary(mod)
sort(abs(modsumm$sites[,1]),decreasing=TRUE)[1:13]
which(modsumm)


mem_null.rda<-rda(t(x$G),mem_pos)
(mem_null.R2a<-RsquareAdj(mem_null.rda)$adj.r.squared)
(mem_null.R2<-RsquareAdj(mem_null.rda)$r.squared)
mem_null.fwd<-forward.sel(t(x$G),mem_pos, adjR2thresh = mem_null.R2, alpha=0.05, nperm=9999)
mem_null.sign<-sort(functional_group1.MEMpos.fwd$order)
mem_null.order<-as.data.frame(mem_pos[c(MEMpos.sign)])

