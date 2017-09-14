library(NAM)
library(bigsnpr)
infos <- readRDS("validation/sim2a.rds")
G <- add_code256(big_copy(t(infos$G), type = "raw"), code = bigsnpr:::CODE_012)
maf <- snp_MAF(G)
ind.col <- which(maf > 0.05)
G2 <- big_copy(G, ind.col = ind.col)
G2
CHR <- infos$chromosome[ind.col]
normalize <- function(x) {
  qx <- ppoints(length(x))
  qnorm(qx[rank(x)])
} 
G3 <- scale(G2[])
rownames(G3)<-make.names(1:1000, unique = TRUE)
colnames(G3)<-make.names(1:ncol(G3), unique = TRUE)
mod<-gwas3(y=infos$phenotype1, gen=G3, fam=infos$pop,
           chr=as.vector(table(CHR)), cov=infos$envi)

dim(G3)
