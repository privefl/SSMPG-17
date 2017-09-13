library(bigsnpr)
# infos <- readRDS("training/sim1a.rds")
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
cov_base <- cbind(normalize(infos$envi),
                  model.matrix(~as.factor(infos$pop)))
cov_add <- NULL
ind <- NULL
G3 <- scale(G2[])

cov3 <- cbind(cov_base, cov_add) 

ans <- rrBLUP::mixed.solve(normalize(infos$phenotype1), Z = G3, X = cov3, SE = TRUE,
                           K = crossprod(G3))

score <- ans$u / ans$u.SE
length(score)
pval <- pchisq(score^2, df = 1, lower.tail = FALSE)
plot(-log10(pval))


new.ind <- which.min(pval)
ind <- c(ind, new.ind)
cov_add <- scale(G2[, ind, drop = FALSE])
G3[, new.ind] <- rnorm(1000)
ind.col[ind]


snps13 <- ind.col[ind]
write(snps13, "gwas-eval.txt", ncolumns = 1)
