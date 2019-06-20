source('functions_ztable.R')

setwd("/Users/mila/t2d_cancer/cluster_analysis")

final_loci_list <- read.table("established_loci_combined.txt", header = TRUE, sep = "\t", dec = ".", na.strings = "NA", stringsAsFactors = FALSE, fill = TRUE)
phenotype_table <- readxl::read_xlsx("phenotype_table.xlsx")

# this is a rudiment, but the code will complain if you skip this step. I'll remove it later this week.
phenotype_table_writein <- phenotype_table
z_table <- read.table("/Users/mila/t2d_cancer/cluster_analysis/cluster_analysis/z_table_2019-02-18 18:24:21.txt", header = TRUE, sep = "\t", dec = ".", na.strings = "NA", stringsAsFactors = FALSE, fill = TRUE)


# For the first time only
standard_ss <- normalise_summary_stats('/Users/mila/t2d_cancer/gwas_ss_sources/METAANALYSIS_DIAGRAM_SE1.txt',
                                       sep = NA,
                                       phenotype = 't2d',
                                       chrpos_col = 'chrpos',
                                       ea_col = 'Allele1',
                                       nea_col = 'Allele2',
                                       beta_col = 'Effect',
                                       se_col = 'StdErr')

# things here will need tweaking
z_table <- merge(final_loci_list, standard_ss, by.x = 'chr.pos', by.y = 'chrpos', all.x = TRUE)
not_found <- sum(is.na(z_table[, c('z_t2d')]))
z_table$V2 <- NULL
z_table$original_freq <- NULL
z_table$new_freq <- NULL
colnames(z_table)[6] <- 't2d_effect_allele'
colnames(z_table)[7] <- 't2d_other_allele'

for (i in 1:nrow(z_table)) {
  if (!is.na(z_table$z_t2d[i]) && z_table$z_t2d[i] < 0) {
    temp <- z_table$t2d_other_allele[i] 
    z_table$t2d_other_allele[i] <- z_table$t2d_effect_allele[i]
    z_table$t2d_effect_allele[i] <- temp
    z_table$z_t2d[i] <- -z_table$z_t2d[i]
  }
  else if (!is.na(z_table$z_t2d[i]) && z_table$z_t2d[i] > 0) {
  }
}


for (i in 1:nrow(phenotype_table)) {
  # no need to change anything here
  standard_ss <- normalise_summary_stats(phenotype_table[i, c('pathname')], 
                                         sep = phenotype_table[i, c('sep')], 
                                         phenotype = phenotype_table[i, c('phenotype')], 
                                         chrpos_col = phenotype_table[i, c('chrpos_col')], 
                                         chr_col = phenotype_table[i, c('chr_col')],
                                         pos_col = phenotype_table[i, c('pos_col')],
                                         rsid_col = phenotype_table[i, c('rsid_col')],
                                         ea_col = phenotype_table[i, c('ea_col')], 
                                         nea_col = phenotype_table[i, c('nea_col')], 
                                         beta_col = phenotype_table[i, c('beta_col')], 
                                         se_col = phenotype_table[i, c('se_col')],
                                         z_col = phenotype_table[i, c('z_col')]
  )
  z_table <- add_phenotype(z_table, standard_ss)
}



# copied from the beginning:
# pathname <- '../gwas_ss/Leptin/Leptin_Adjusted_for_BMI.txt'
# pathname <- "../gwas_ss/PAI.res"
# pathname <- '../gwas_ss/METAANALYSIS_DIAGRAM_SE1.txt'
# summary_stats <- read.table(pathname, sep = '\t', header = TRUE, stringsAsFactors = FALSE, na.strings = "NA")

# standard_ss <- normalise_summary_stats('/rds/general/user/lz5515/home/WORK/t2d_cancer/gwas_ss/METAANALYSIS_DIAGRAM_SE1.txt',
#                                        sep = NA,
#                                        phenotype = 't2d',
#                                        chrpos_col = 'chrpos',
#                                        ea_col = 'Allele1',
#                                        nea_col = 'Allele2',
#                                        beta_col = 'Effect',
#                                        se_col = 'StdErr')

# z_table <- add_phenotype(z_table, standard_ss)
# homa_b_imputed <- read.table("../gwas_ss/HOMA_B_0.6.txt", header = TRUE, sep = "\t", dec = ".", na.strings = "NA", stringsAsFactors = FALSE, fill = TRUE)
# homa_b_imputed$chrpos <- paste(homa_b_imputed$chr, homa_b_imputed$pos, sep = ':')
# check_b_06 <- merge(final_loci_list_37_unique, homa_b_imputed, by.x = 'chr.pos', by.y = 'chrpos', all.x = TRUE)
# 
# homa_ir_imputed <- read.table("../gwas_ss/HOMA_IR_0.6.txt", header = TRUE, sep = "\t", dec = ".", na.strings = "NA", stringsAsFactors = FALSE, fill = TRUE)
# homa_ir_imputed$chrpos <- paste(homa_ir_imputed$chr, homa_ir_imputed$pos, sep = ':')
# check_06 <- merge(final_loci_list_37_unique, homa_ir_imputed, by.x = 'chr.pos', by.y = 'chrpos', all.x = TRUE)
