
str(infos <- readRDS("cichlid_data/cichlid.rds"))

# rle(infos$species)
# rle(infos$place)
# rle(infos$place)


library(bigsnpr)
library(ggplot2)

G <- add_code256(big_copy(t(infos$G), type = "raw"), code = bigsnpr:::CODE_012)
summary(maf <- snp_MAF(G))

ind.maf.ok <- which(maf > 0.05)
G <- big_copy(G, ind.col = ind.maf.ok)

CHR <- readr::parse_number(infos$chromosome)[ind.maf.ok]
POS <- infos$position[ind.maf.ok]

ind.keep <- snp_clumping(G, CHR)

obj.svd <- big_SVD(G, snp_scaleBinom(), ind.col = ind.keep, k = 20)
plot(obj.svd, type = "loadings", loadings = 1:6, coeff = 0.5)

plot(obj.svd)
plot(obj.svd, type = "scores", coeff = 1.5, scores = 1:2) + 
  aes(color = infos$fishing_pressure[-outlier]) + 
  labs(color = "Pop")


gwas <- snp_pcadapt(G, obj.svd$u[, 1:3])
snp_qq(gwas)
gwas.gc <- snp_gc(gwas)
snp_qq(gwas.gc)
snp_manhattan(gwas.gc, CHR, POS, dist.sep.chrs = 1e5)

outlier <- c(46, 47)
G2 <- big_copy(G, ind.row = rows_along(G)[-outlier])
G2




ind.keep <- snp_clumping(G2, CHR)

obj.svd <- big_SVD(G2, snp_scaleBinom(), ind.col = ind.keep, k = 20)
plot(obj.svd, type = "loadings", loadings = 1:6, coeff = 0.5)

plot(obj.svd)
plot(obj.svd, type = "scores", coeff = 1.5, scores = 1:2) + 
  aes(color = infos$fishing_pressure[-outlier]) + 
  labs(color = "Pop")


gwas2 <- snp_pcadapt(G2, obj.svd$u[, 1:2])
snp_qq(gwas2)
gwas2.gc <- snp_gc(gwas2)
snp_qq(gwas2.gc)
snp_manhattan(gwas2.gc, CHR, POS, dist.sep.chrs = 1e5)


obj.svd <- snp_autoSVD(G2, CHR, POS)
