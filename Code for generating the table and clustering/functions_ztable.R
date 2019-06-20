# all_summary_stats <- 

# rsid = TRUE -> rsIDs are present, merge by rs_id
# rsid = FALSE -> merge by chr.pos
# chr:pos = TRUE -> no column concatenation is needed
# chr:pos = FALSE -> concatenate position and chromosome# columns
# Feed it column numbers, not names
#pathname <- "../gwas_ss/breast_cancer/bcac_meta.txt"


# STEP 1 ###########################################################################
rename_columns <- function(summary_stats, ea_col, nea_col, beta_col, se_col, z_col = NA, phenotype) {
  ea_col_number <- which(names(summary_stats) == as.character(ea_col))
  nea_col_number <- which(names(summary_stats) == as.character(nea_col))
  colnames(summary_stats)[ea_col_number] <- 'effect_allele'
  colnames(summary_stats)[nea_col_number] <- 'other_allele'
  summary_stats[, ea_col_number] <- toupper(summary_stats[, ea_col_number])
  summary_stats[, nea_col_number] <- toupper(summary_stats[, nea_col_number])
  
  if (is.na(z_col)) {
    beta_col_number <- which(names(summary_stats) == as.character(beta_col))
    se_col_number <- which(names(summary_stats) == as.character(se_col))
    colnames(summary_stats)[beta_col_number] <- 'beta'
    colnames(summary_stats)[se_col_number] <- 'se'
  }
  else {
    z_column <- paste("z", as.character(phenotype), sep = "_")
    z_col_number <- which(names(summary_stats) == as.character(z_col))
    colnames(summary_stats)[z_col_number] <- z_column
  }
  return(summary_stats)
}

# STEP 2 #########################################################################
reduce_ss <- function(summary_stats_renamed, rsid_col, chr_col, pos_col, chrpos_col, z_col = NA, phenotype) {
  chrpos_col_number <- which(names(summary_stats_renamed) == as.character(chrpos_col))
  colnames(summary_stats_renamed)[chrpos_col_number] <- 'chrpos'
  if (is.na(rsid_col)) {
    # if FALSE, we need to merge by positions; check that they are one column
    if (is.na(chrpos_col)) {
      # if FALSE, we need to concatenate them first
      chr_col_number <- which(names(summary_stats_renamed) == as.character(chr_col))
      pos_col_number <- which(names(summary_stats_renamed) == as.character(pos_col))
      summary_stats_renamed$chrpos <- paste(summary_stats_renamed[ , chr_col_number], summary_stats_renamed[ , pos_col_number], sep = ':')
      if (is.na(z_col)) {
        reduced_summary_stats <- summary_stats_renamed[, c('chrpos', 'effect_allele', 'other_allele', 'beta', 'se')]
      } else {
        z_column <- paste("z", as.character(phenotype), sep = "_")
        reduced_summary_stats <- summary_stats_renamed[, c('chrpos', 'effect_allele', 'other_allele', as.character(z_column))]
        if (class(reduced_summary_stats[, as.character(z_column)]) != "numeric") {
          reduced_summary_stats[, as.character(z_column)] <- as.numeric(reduced_summary_stats[, as.character(z_column)])
        }
      }
    }
    else {
      if (is.na(z_col)){
        chrpos_col_number <- which(names(summary_stats_renamed) == as.character(chrpos_col))
        colnames(summary_stats_renamed)[chrpos_col_number] <- 'chrpos'
        reduced_summary_stats <- summary_stats_renamed[, c('chrpos', 'effect_allele', 'other_allele', 'beta', 'se')]
      } else {
        z_column <- paste("z", as.character(phenotype), sep = "_")
        chrpos_col_number <- which(names(summary_stats_renamed) == as.character(chrpos_col))
        colnames(summary_stats_renamed)[chrpos_col_number] <- 'chrpos'
        reduced_summary_stats <- summary_stats_renamed[, c('chrpos', 'effect_allele', 'other_allele', as.character(z_column))]
        if (class(reduced_summary_stats[, as.character(z_column)]) != "numeric") {
          reduced_summary_stats[, as.character(z_column)] <- as.numeric(reduced_summary_stats[, as.character(z_column)])
        }
      }
    }
  }
  else {
    if (is.na(z_col)){
      rsid_col_number <- which(names(summary_stats_renamed) == as.character(rsid_col))
      colnames(summary_stats_renamed)[rsid_col_number] <- 'rsid'
      reduced_summary_stats <- summary_stats_renamed[ , c('rsid', 'effect_allele', 'other_allele', 'beta', 'se')]
    } else {
      z_column <- paste("z", as.character(phenotype), sep = "_")
      rsid_col_number <- which(names(summary_stats_renamed) == as.character(rsid_col))
      colnames(summary_stats_renamed)[rsid_col_number] <- 'rsid'
      reduced_summary_stats <- summary_stats_renamed[ , c('rsid', 'effect_allele', 'other_allele', as.character(z_column))]
      if (class(reduced_summary_stats[, as.character(z_column)]) != "numeric") {
        reduced_summary_stats[, as.character(z_column)] <- as.numeric(reduced_summary_stats[, as.character(z_column)])
      }
    }
  }
  return(reduced_summary_stats)
}

# STEP 3 #########################################################################
compute_z <- function(reduced_summary_stats, phenotype) {
  beta_col_number <- which(names(reduced_summary_stats) == 'beta')
  se_col_number <- which(names(reduced_summary_stats) == 'se')
  if (class(reduced_summary_stats[, beta_col_number]) != 'numeric') {
    reduced_summary_stats[, beta_col_number] <- as.numeric(reduced_summary_stats[, beta_col_number])
    cat("Number of unavailable BETA values:", sum(is.na(reduced_summary_stats[, beta_col_number])), "\n")
  }
  if (class(reduced_summary_stats[, se_col_number]) != 'numeric') {
    reduced_summary_stats[, se_col_number] <- as.numeric(reduced_summary_stats[, se_col_number])
    cat("Number of unavailable SE values: ", sum(is.na(reduced_summary_stats[, se_col_number])), "\n")
  }
  z_column <- paste("z", as.character(phenotype), sep = "_")
  reduced_summary_stats$z <- reduced_summary_stats[, beta_col_number]/reduced_summary_stats[, se_col_number]
  z_col_number <- which(names(reduced_summary_stats) == 'z')
  colnames(reduced_summary_stats)[z_col_number] <- as.character(z_column)
  
  reduced_summary_stats$beta <- NULL
  reduced_summary_stats$se <- NULL
  
  return(reduced_summary_stats)
}

# #########################################################################
normalise_summary_stats <- function(pathname, phenotype, sep = '\t', rsid_col = NA, chr_col = NA, pos_col = NA, chrpos_col = NA, ea_col = NA, nea_col = NA, beta_col = NA, se_col = NA, z_col = NA, normalize_z = FALSE) {
  if (is.na(sep)){
    summary_stats <- read.table(as.character(pathname), sep = '\t', header = TRUE, stringsAsFactors = FALSE, na.strings = "NA", fill = TRUE)
  }
  else if (sep == ',') {
    summary_stats <- read.csv(as.character(pathname), header = TRUE, stringsAsFactors = FALSE, na.strings = "NA", fill = TRUE)
  }
  else if (sep == 'space'){
    summary_stats <- read.table(as.character(pathname), sep = ' ', header = TRUE, stringsAsFactors = FALSE, na.strings = "NA", fill = TRUE)
  }
  
  cat("Phenotype:", as.character(phenotype), "\n")
  cat("Number of SNPs:", nrow(summary_stats), "\n")
  phenotype_table_writein[i, c('n_snps')] <- nrow(summary_stats)
  summary_stats_renamed <- rename_columns(summary_stats, ea_col, nea_col, beta_col, se_col, z_col, phenotype)
  summary_stats_reduced <- reduce_ss(summary_stats_renamed, rsid_col, chr_col, pos_col, chrpos_col, z_col, phenotype)
  if (is.na(z_col)){
    summary_stats_standard <- compute_z(summary_stats_reduced, phenotype)
  } else {
    summary_stats_standard <- summary_stats_reduced
    cat(colnames(summary_stats_standard), "\n")
  }
  
  if (normalize_z == TRUE) {
    
  }
  
  # always same order: 1) snp_id, 2) ea, 3) nea, 4) z
  return(summary_stats_standard)
}



add_phenotype <- function(z_table, summary_stats_standard) {
  #extract_phenotype
  # head(summary_stats_standard, 10)
  z_phenotype <- colnames(summary_stats_standard)[4]
  phenotype <- gsub('z_', '', z_phenotype)
  cat(phenotype, '\n')
  if (phenotype != 't2d') {
    if (colnames(summary_stats_standard)[1] == 'chrpos'){
      z_table_new <- merge(z_table, summary_stats_standard, by.x = 'chr.pos', by.y = 'chrpos', all.x = TRUE)
      #head(z_table_new, 10)
    }
    if (colnames(summary_stats_standard)[1] == 'rsid'){
      z_table_new <- merge(z_table, summary_stats_standard, by.x = 'rsid', by.y = 'rsid', all.x = TRUE)
      #head(z_table_new, 10)
    }
    z_col_num <- which(names(z_table_new) == as.character(z_phenotype)) 
    not_found <- sum(is.na(z_table_new[ ,z_col_num]))
    # or maybe below it's better to use locilist
    snps_not_found <- z_table_new[is.na(z_table_new[ , z_col_num]), ] 
    cat("Number of missing snps:", not_found, "\n")
    phenotype_table_writein$missing_snps[i] <- not_found
    phenotype_table_writein$snps_found[i] <- nrow(final_loci_list) - not_found
    count <- 0
    # Don't touch anything here
    for (j in 1:nrow(z_table_new)) {
      # a check that this loci was found in this summary statistics
      if (is.na(z_table_new[j, c('t2d_effect_allele')]) == FALSE){
        # a check that we have information for what allele is effective
        if (is.na(z_table_new[j, c('effect_allele')]) == FALSE) {
          # a check that effective allele we operate with is different from effective allele in summary statistics
          if (z_table_new[j, c('effect_allele')] != z_table_new[j, c('t2d_effect_allele')]) {
            # we don't need to change alleles actually, only z-score
            z_table_new[j, z_col_num] <- -z_table_new[j, z_col_num]
            count <- count + 1
          }
        }
      } 
      else {
        z_table_new[j, c('t2d_effect_allele')] <- z_table_new[j, c('effect_allele')]
        z_table_new[j, c('t2d_other_allele')] <- z_table_new[j, c('other_allele')]
      }
    }
  }
  else {
    z_table <- merge(final_loci_list_37_unique, standard_ss, by.x = 'chr.pos', by.y = 'chrpos', all.x = TRUE)
    write.table(z_table, file=paste("z_table", Sys.Date(), sep = "_"), sep = "\t", quote = FALSE)
    not_found <- sum(is.na(z_table[, c('z_t2d')]))
    z_table$V2 <- NULL
    z_table$original_freq <- NULL
    z_table$new_freq <- NULL
    colnames(z_table)[c('effect_allele')] <- 't2d_effect_allele'
    colnames(z_table)[c('other_allele')] <- 't2d_other_allele'
  }
  cat('The number of alleles aligned:', count, '\n')
  phenotype_table_writein$number_of_alleles_aligned[i] <- count
  cat("\n")
  z_table_new$effect_allele <- NULL
  z_table_new$other_allele <- NULL
  write.table(z_table_new, file=paste("z_table_", Sys.time(), ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  return(z_table_new)
}
