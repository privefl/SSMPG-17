
POP <- simu1$fam$pop
PHENO <- simu1$fam$affection

library(foreach)
doParallel::registerDoParallel(3)
test <- foreach(i1 = unique(POP), .combine = 'c') %:%
  foreach(i2 = seq_len(i1 - 1)) %dopar% {
    
    ind <- which(POP == i1 | POP == i2)
    env.normed <- normalize(simu1$fam[["env"]][ind])
    
    gwas3 <- big_univLinReg(G, normalize(PHENO[ind]), ind.train = ind, 
                            covar.train = as.matrix(env.normed))
    gwas3[is.na(gwas3$score), "score"] <- 0
    predict(snp_gc(gwas3))
  }
doParallel::stopImplicitCluster()

tmp <- -colMeans(do.call("rbind", test))
plot(tmp)

