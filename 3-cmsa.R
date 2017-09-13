betas <- big_CMSA(
  FUN = big_spLinReg,
  X = G, 
  feval = function(pred, target) -sum((pred - target)^2), 
  y.train = simu1$fam$affection,
  covar.train = cov, 
  alpha = 0.5
)



ind <- which(betas[cols_along(G)] != 0)
tmp <- pval3[ind]
ind[tmp < 0.01]

snps3 <- which(infos$position %in% simu1$map$physical.pos[ind[tmp < 0.01]])
snps4 <- union(snps3, snps2)

write(snps4, "gwas-clumping-union-cmsa.txt", ncolumns = 1)
# 0.36052	0.53846	0.75862	5