library(bigsnpr)

bed1 <- snp_plinkQC(
  prefix.in = "training/sim1a",
  file.type = "--vcf",
  maf = 0.05
)

# bed2 <- snp_plinkIBDQC(
#   bedfile.in = bed1,
#   do.blind.QC = FALSE,
#   pi.hat = 0.55
# )
# 
# library(ggplot2)
# 
# ggplot(df.IBD) + 
#   geom_point(aes(Z0, Z1, color = RT))
# 
# ggplot(df.IBD) + 
#   geom_point(aes(seq_along(PI_HAT), PI_HAT))

bed2 <- snp_plinkIBDQC(
  bedfile.in = bed1,
  pi.hat = 0.55
)

snp_readBed(bed2, backingfile = "back/simu1")


