
proj <- PCs %*% t(diff(test2$centers[2:3, ]))

qplot(PCs[, 1], PCs[, 2], col = as.vector(proj)) + 
  viridis::scale_color_viridis()


gwas <- big_univLinReg(G2, proj)
snp_qq(gwas)
gwas.gc <- snp_gc(gwas)
snp_qq(gwas.gc)
snp_manhattan(gwas.gc, CHR, POS, dist.sep.chrs = 1e5)

tmp2 <- tmp
diag(tmp2) <- 0
max(abs(tmp2))
