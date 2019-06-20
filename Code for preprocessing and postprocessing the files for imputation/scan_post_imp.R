HDL_imputed <- read.table("/Users/mila/ssimp_web/GWAS_HDL_imp_qual07.txt", header = TRUE, sep = "\t", stringsAsFactors = F, fill = T, na.strings = NA)
HDL_imputed$P.imp <- as.numeric(HDL_imputed$P.imp)

HDL <- read.table("/Users/mila/ssimp_web/GWAS_HDL.txt", header = TRUE, sep = "\t", stringsAsFactors = F, fill = T, na.strings = NA)

HDL$SNP_hg19 <- unlist(lapply(strsplit(HDL$SNP_hg19, "r"), "[", 2))
HDL$chr <- unlist(lapply(strsplit(HDL$SNP_hg19, ":"), "[", 1))
HDL_chr11 <- HDL[HDL$chr == 11, ]
range(HDL_chr6$P.value)

HDL_imputed_chr11 <- HDL_imputed[HDL_imputed$chr == 11, ]
range(HDL_imputed_chr6$P.imp)

for (i in 1:22){
  tmp <- HDL[HDL$chr == i, ]
  tmp_imp <- HDL_imputed[HDL_imputed$chr == i, ]
  cat("chrm", i, "\n")
  cat("Min P-value before imputation:", min(tmp$P.value), "\n")
  cat(tmp[which.min(tmp$P.value), c("rsid")], "\n")
  cat("Min P-value after imputation:", min(tmp_imp$P.imp), "\n")
  cat(tmp_imp[which.min(tmp$P.value), c("source")], "\n")
  cat("\n")
}


