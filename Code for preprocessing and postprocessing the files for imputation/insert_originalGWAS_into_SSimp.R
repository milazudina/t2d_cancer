
#Metabochip Model-1
meta_imp <- read.table("/project/silk/data/GWAS_summary-stats_imputed/RG/16042019_rerun/Metabochip_3datasets_Model_1_imputed_qual04.txt", header=TRUE, sep="\t")
meta <- read.table("/project/silk/data/GWAS_summary-stats_imputed/RG/original_files/Metabochip_3datasets_Model_1_2019Mar_1.tbl", header=TRUE, sep="\t")

meta_imp <- meta_imp[!meta_imp$source == "source",] #remove left in header lines

meta_imp$MarkerName <- paste0("chr", meta_imp$chr, ":", meta_imp$pos) #make chr:pos IDs for SSimp
meta_imp2 <- meta_imp[!meta_imp$SNP %in% meta$MarkerName & !meta_imp$MarkerName %in% meta$MarkerName,] #exclude variants present in the original GWAS from SSImp 
meta_imp2[16] <- NULL

meta$Allele1 <- toupper(as.character(meta$Allele1)) #capitalize alleles
meta$Allele2 <- toupper(as.character(meta$Allele2))
meta$maf <- ifelse(meta$Freq1 <= 0.5, meta$Freq1, 1-meta$Freq1)

meta <- meta[,c(1:3,8:10,17)] #keep MarkerName, allele, effect and P.value columns from original GWAS

meta$z_imp <- meta$Effect / meta$StdErr  #compute Z score for original GWAS P.values

#add misiing columns before binding
meta$chr <- NA
meta$pos <- NA
meta$source <- "GWAS"
meta$r2.pred <- NA
meta$lambda <- NA
meta$Z_reimputed <- NA
meta$r2_reimputed <- NA
meta$N.imp <- NA
meta$bst.imp <- NA
names(meta)[1] <- "SNP"
names(meta)[6] <- "P.imp"

meta <- meta[, colnames(meta_imp2)]

#bind original GWAS and SSimp
meta_all <- rbind(meta, meta_imp2)

write.table(meta_all, file="/project/silk/users/au615/RG/Metabochip_3datasets_Model_1_imputed_qual04_originalGWAS.txt", sep="\t", row.names=FALSE, quote=FALSE)

rm(list=ls())


#Metabochip Model-5





#HapMap Model-1



#HapMap Model-5 P1




#HapMap Model-5 P2
