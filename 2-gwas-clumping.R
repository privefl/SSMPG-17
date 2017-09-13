hist(tmp <- simu1$fam$affection, freq = FALSE)
test <- fitdistrplus::fitdist(simu1$fam$affection, "norm")

lines(density(rnorm(1e5, test$estimate)), col = "red")


normalize <- function(x) {
  qx <- ppoints(length(x))
  qnorm(qx[rank(x)])
} 
plot(normalize(tmp), tmp)

tukey <- function(x) {
  iqr <- 1.5 * IQR(x)
  lower <- quantile(x, 0.25) - iqr
  upper <- quantile(x, 0.75) + iqr
  which(x < lower | x > upper)
}
tukey(simu1$fam[["env"]])

env.normed <- normalize(simu1$fam[["env"]])
hist(env.normed)
which(is.infinite(env.normed))
cov3 <- cbind(as.matrix(env.normed))  #, model.matrix(~pop.ju))

gwas3 <- big_univLinReg(G, normalize(simu1$fam$affection), covar.train = cov3)
snp_qq(gwas3)
gwas3.gc <- snp_gc(gwas3)
snp_qq(gwas3.gc)
plot(gwas$score, gwas3$score)
abline(0, 1, col = "red")

snp_manhattan(gwas3.gc, CHR, unlist(tapply(CHR, CHR, seq_along)), 
              dist.sep.chrs = 50)
snp_manhattan(gwas.gc, CHR, unlist(tapply(CHR, CHR, seq_along)), 
              dist.sep.chrs = 50)


ind.keep <- snp_clumping(G, CHR, S = abs(gwas3.gc$score), size = 100)

snp_qq(snp_gc(gwas3.gc[ind.keep, ]))

plot(qvalue::qvalue(pval3[ind.keep]))
sort(pval3[ind.keep])

snps2 <- which(infos$position %in% 
                 simu1$map$physical.pos[ind.keep[order(pval3[ind.keep])[1:8]]])
  
write(snps2, "gwas-clumping2.txt", ncolumns = 1)
# 0.49029	0.38462	0.375	

