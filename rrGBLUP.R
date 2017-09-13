# infos <- readRDS("training/sim1a.rds")
infos <- readRDS("validation/sim2a.rds")

maf <- snp_MAF(G)
ind.col <- which(maf > 0.05)
G2 <- big_copy(G, ind.col = ind.col)
G2

CHR <- infos$chromosome[ind.col]

cov_base <- cbind(normalize(infos$envi),
                  model.matrix(~as.factor(infos$pop)))
cov_add <- NULL
ind <- NULL
G3 <- G2[]


cov3 <- cbind(cov_base, cov_add) 

ans <- mixed.solve(normalize(infos$phenotype1), Z = G3, X = cov3, SE = TRUE)

score <- ans$u / ans$u.SE
length(score)
pval <- pchisq(score^2, df = 1, lower.tail = FALSE)
plot(-log10(pval))

# gwas3 <- big_univLinReg(G2, normalize(infos$phenotype1), covar.train = cov3)
# gwas3[is.na(gwas3$score), "score"] <- 0
# gwas3.gc <- snp_gc(gwas3)


# ind.keep <- snp_clumping(G2, CHR, S = -pval, size = 100, 
#                          exclude = ind)

new.ind <- which.min(pval)
ind <- c(ind, new.ind)
cov_add <- G2[, ind, drop = FALSE]
G3[, new.ind] <- rnorm(1000)
ind.col[ind]


snps12 <- ind.col[ind]
write(snps12, "gwas-5perc.txt", ncolumns = 1)
