dat <- readRDS("sim1a.rds")
plot(dat$phenotype1, dat$phenotype2)
plot(dat$phenotype1, dat$envi)
plot(dat$x, dat$y, col = dat$pop, pch = dat$pop %% 26 + 1, cex = 0.6)

library(bigsnpr)
library(plotly)

G <- add_code256(big_copy(t(dat$G), type = "raw"), code = bigsnpr:::CODE_012)
CHR <- dat$chromosome
POS <- dat$position

obj.svd <- snp_autoSVD(G, CHR, POS, size = 200, roll.size = 10, int.min.size = 2)

plot(obj.svd)
plot(obj.svd, type = "scores", scores = 7:8)

pcadapt1 <- snp_pcadapt(G, obj.svd$u[, 1])
snp_qq(pcadapt1)
pcadapt1 <- snp_gc(pcadapt1)
snp_qq(pcadapt1)
(snp_manhattan(pcadapt1, CHR, unlist(tapply(CHR, CHR, seq_along)), 
              dist.sep.chrs = 200) +
  ggplot2::ylab("-log10(pval)")) %>%
  ggplotly()

write(order(predict(pcadapt1))[1:4], "submit-top4.test")

pcadapt6 <- snp_pcadapt(G, obj.svd$u[, 1:6])
snp_qq(pcadapt6)
pcadapt6 <- snp_gc(pcadapt6)
snp_qq(pcadapt6)
which.max(abs(pcadapt6$score))
(snp_manhattan(pcadapt1, CHR, unlist(tapply(CHR, CHR, seq_along)), 
               dist.sep.chrs = 200) +
    ggplot2::ylab("-log10(pval)")) %>%
  ggplotly()

pval <- predict(pcadapt6, log10 = FALSE)
ind.keep <- snp_clumping(G, CHR, S = pcadapt6$score, size = 50)
write(ind.keep[order(pval[ind.keep])][1:10], "submit-top10-clumping.txt")
