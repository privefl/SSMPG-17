
yware_scale <- function(y.row) {
  function(X, 
           ind.row = rows_along(X), 
           ind.col = cols_along(X)) {
    
    data.frame(
      center = big_scale()(X, ind.row, ind.col)$center,
      scale = 1 / big_univLinReg(X, y.row, ind.row, ind.col)$estim
    )
  }
}

test <- yware_scale(simu1$fam$affection)(G)
head(test)
obj.svd <- snp_autoSVD(
  G, CHR, POS,
  fun.scaling = yware_scale(simu1$fam$affection),
  roll.size = 10,
  int.min.size = 10
)

ind.keep <- snp_clumping(G, CHR)
obj.svd <- big_randomSVD(G, yware_scale(simu1$fam$affection), k = 20)
plot(obj.svd)
plot(obj.svd, type = "scores") + 
  aes(color = simu1$fam$env)
plot(obj.svd, type = "scores", scores = 5:6)

gwas <- snp_pcadapt(G, obj.svd$u[, 1:6])
snp_qq(gwas)
snp_manhattan(gwas, CHR, unlist(tapply(CHR, CHR, seq_along)), 
              dist.sep.chrs = 50) + 
  scale_y_log10()

plot(qvalue::qvalue(predict(gwas, log10 = FALSE)))

ind <- order(predict(gwas))[1:100]
snps6 <- which(infos$position %in% simu1$map$physical.pos[ind])
snps5 %in% snps6   # only 5 /8
snps7 <- union(snps5[-9], snps6)
plot(sort(predict(gwas))[1:100])

write(snps7, "pcadapt-yware.txt", ncolumns = 1)
