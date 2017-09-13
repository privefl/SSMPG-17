library(bigsnpr)
library(tidyverse)

simu1 <- snp_attach("back/simu1.rds")

infos <- readRDS("training/sim1a.rds")
ind <- match(simu1$fam$sample.ID, paste0("i", infos$id))
infos.ind <- data.frame(infos[c("x", "y", "phenotype1", "envi")])[ind, ]
simu1$fam <- mutate(simu1$fam,
                    x = infos[["x"]][ind],
                    y = infos[["y"]][ind],
                    affection = infos[["phenotype1"]][ind],
                    env = infos[["envi"]][ind],
                    pop = infos[["pop"]][ind])
  
simu1 <- snp_save(simu1)
G <- simu1$genotypes
G
CHR <- simu1$map$chromosome
POS <- simu1$map$physical.pos

obj.svd <- snp_autoSVD(G, CHR, POS)
plot(obj.svd)
plot(obj.svd, type = "scores") + 
  aes(color = simu1$fam$affection)

plot(obj.svd$u[, 1], simu1$fam$affection)

ggplot(simu1$fam) + 
  geom_point(aes(x, y, color = affection)) + 
  viridis::scale_color_viridis()

ggplot(simu1$fam) + 
  geom_point(aes(x, y, color = env)) + 
  viridis::scale_color_viridis()

rgl::plot3d(x = simu1$fam[c("x", "y", "env")], size = 5,
            col = simu1$fam$affection)


cov <- cbind(as.matrix(simu1$fam[c("x", "y", "env")]),
             model.matrix(~as.factor(simu1$fam$pop)))
gwas <- big_univLinReg(G, simu1$fam$affection, covar.train = cov)
snp_qq(gwas)
gwas.gc <- snp_gc(gwas)
snp_qq(gwas.gc)

snp_manhattan(gwas.gc, CHR, unlist(tapply(CHR, CHR, seq_along)), 
              dist.sep.chrs = 50)

pval <- predict(gwas.gc, log10 = FALSE)

plot(qvalue::qvalue(pval))
snps <- which(infos$position %in% simu1$map$physical.pos[order(pval)[1:40]])
write(snps, "gwas-cov42.txt")  # 6/13 good for 41 submitted

## With pop of Julian
pop.ju <- readRDS("gen_clusters.rds")
cov2 <- cbind(as.matrix(simu1$fam[c("env")]),
              model.matrix(~pop.ju))

gwas2 <- big_univLinReg(G, simu1$fam$affection, covar.train = cov2)
snp_qq(gwas2)
gwas2.gc <- snp_gc(gwas2)
snp_qq(gwas2.gc)
plot(gwas$score, gwas2$score)
abline(0, 1, col = "red")

snp_manhattan(gwas2.gc, CHR, unlist(tapply(CHR, CHR, seq_along)), 
              dist.sep.chrs = 50)

pval <- predict(gwas2.gc, log10 = FALSE)

plot(qvalue::qvalue(pval))
snps <- which(infos$position %in% simu1$map$physical.pos[order(pval)[1:40]])
write(snps, "gwas-cov42.txt")  # 6/13 good for 41 submitted
