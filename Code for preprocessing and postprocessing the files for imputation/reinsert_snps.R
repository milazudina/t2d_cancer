gwas_imp_path <- "/rdsgpfs/general/project/eph-prokopenko-lab-silk/live/data/GWAS_summary-stats_imputed/MAGIC_fg_combined/GWAS_FG_imp_qual07.txt"
gwas_path <- "/rdsgpfs/general/user/lz5515/home/WORK/magic_sex/fg_gwas_combined.txt"

gwas_imp <- read.table(gwas_imp_path, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, sep = "\t")
gwas <- read.table(gwas_path, header = TRUE, fill = TRUE, stringsAsFactors = FALSE, sep = "\t")

gwas_imp$chr <- as.numeric(gwas_imp$chr)
gwas_imp$pos <- as.numeric(gwas_imp$pos)
gwas_imp$P.imp <- as.numeric(gwas_imp$P.imp)
gwas_imp$N.imp <- as.numeric(gwas_imp$N.imp)
gwas_imp$z_imp <- as.numeric(gwas_imp$z_imp)
gwas_imp$r2.pred <- as.numeric(gwas_imp$r2.pred)
gwas_imp$maf <- as.numeric(gwas_imp$maf)
gwas_imp$Z_reimputed <- NULL
gwas_imp$r2_reimputed <- NULL
gwas_imp$lambda <- NULL
gwas_imp$bst.imp <- NULL
gwas_imp <- na.omit(gwas_imp)

f.z2b <- function(z, af, n){
  # z = imputed z statistics
  # af = allele frequency
  # n = sample size (effective)
  se.b <- 1/sqrt(2* af * (1-af) * n)
  b <- z * se.b
  return(c(b, se.b))
}

output.b.se <- lapply(1:nrow(gwas_imp), function(x) 
  f.z2b(gwas_imp$z_imp[x], gwas_imp$maf[x], gwas_imp$N.imp[x])
)

gwas_imp$beta <- sapply(output.b.se, function(x) x[1])
gwas_imp$se <- sapply(output.b.se, function(x) x[2])

###### ORIGINAL GWAS PART ######
gwas$r2.pred <- 1

# only truly imputed snps will go to gwas_imp2
gwas_SSIMPsnps <- gwas_imp[!gwas_imp$SNP %in% gwas$rsID, ]

# preparing files for merge
colnames(gwas_SSIMPsnps)[colnames(gwas_SSIMPsnps) == c("z_imp")] <- "z"
colnames(gwas_SSIMPsnps)[colnames(gwas_SSIMPsnps) == c("maf")] <- "EAF"
colnames(gwas_SSIMPsnps)[colnames(gwas_SSIMPsnps) == c("N.imp")] <- "N"
colnames(gwas_SSIMPsnps)[colnames(gwas_SSIMPsnps) == c("P.imp")] <- "p"
colnames(gwas_SSIMPsnps)[colnames(gwas_SSIMPsnps) == c("SNP")] <- "rsID"
gwas$source <- "GWAS"

gwas_reinserted <- rbind(gwas_SSIMPsnps, gwas)
colnames(gwas_reinserted)[colnames(gwas_reinserted) == "r2.pred"] <- "imp.qual"

gwas_reinserted_sorted <- gwas_reinserted[order(gwas_reinserted$rsID), ]
gwas_reinserted_sorted$beta <- signif(gwas_reinserted_sorted$beta, digits = 6)
gwas_reinserted_sorted$se <- signif(gwas_reinserted_sorted$se, digits = 6)

gwas_reinserted_sorted <- gwas_reinserted_sorted[gwas_reinserted_sorted$rsID != ".", ]

write.table(gwas_reinserted_sorted, "/project/silk/data/GWAS_summary-stats_imputed/MAGIC_fg_combined/GWAS_FG_imp_post.txt", quote = F, sep = '\t', row.names = F)


