# install.packages("http://mouse.cs.ucla.edu/emma/emma_1.1.2.tar.gz", 
#                  repos = NULL, type = "source")
# devtools::install_github("Gregor-Mendel-Institute/mlmm")
library(mlmm)
# devtools::install_github("privefl/bigsnpr")
library(bigsnpr)

names.fake <- paste0("I", rows_along(G2))
K <- big_tcrossprodSelf(G2, snp_scaleBinom())[]
dimnames(K) <- list(names.fake, names.fake)

K[1:5, 1:5]
X <- G2[]
dimnames(X) <- list(names.fake, ind.col)
names(proj) <- names.fake

mygwas <- mlmm::mlmm(Y = proj, X = X, K = K,
                   nbchunks = 2,
                   maxsteps = 100)

snp_infos <- data.frame(
  SNP = ind.col,
  Chr = CHR, 
  Pos = POS
)

plot_step_table(mygwas,'extBIC') # EBIC plot
plot_step_table(mygwas,'maxpval') # mbonf criterion plot
plot_step_RSS(mygwas) # % variance plot
plot_fwd_GWAS(mygwas,1,snp_infos,0.1,main="step 1") # 1st mlmm step plot
plot_fwd_GWAS(mygwas,2,snp_infos,0.1,main="step 2") # 2nd mlmm step plot
plot_fwd_GWAS(mygwas,3,snp_infos,0.1,main="step 3") # 3rd mlmm step plot
plot_fwd_GWAS(mygwas,3,snp_infos,0.1,main="step 4") # 3rd mlmm step plot

plot_fwd_GWAS(mygwas,3,snp_infos,0.1,main="step 5") # 3rd mlmm step plot
plot_opt_GWAS(mygwas,'extBIC',snp_infos,0.1,main="optimal (EBIC)") # optimal step according to eBIC plot
plot_opt_GWAS(mygwas,'mbonf',snp_infos,0.1,main="optimal (mBonf)") # optimal step according to mbonf plot
