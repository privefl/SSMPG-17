library(NAM)
library(bigsnpr)
gwas3(y=infos$phenotype1, gen=G3, fam=infos$pop, chr=as.vector(table(CHR)), cov=infos$envi)
