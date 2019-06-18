options(scipen = 999)
setwd("/Users/mila/ssimp_web/")
library(data.table)

#hapmap_list <- read.table("/Users/mila/ssimp_web/hapmap_snplist.txt", header = T, stringsAsFactors = F, fill = T, na.strings = NA)

# User inputs
phenotype <- 'Adipogen'
ss_path <- '/Users/mila/ssimp_web/GWAS_Adipogen.txt'
positions_added_manually <- T
input_path <- "/rds/general/user/lz5515/home/WORK/ssimp_web/"
imprangeSaveTo <- paste(phenotype, "_imprange.txt", sep = "")
gwasSaveTo <- paste("GWAS_", phenotype,".txt", sep = "")
scriptSaveTo <- paste("ssimp_", phenotype, "_array.sh", sep = "")
output_path <- paste(input_path, phenotype, "/", sep = "")

ss_file <- read.table(ss_path, header = !positions_added_manually, stringsAsFactors = F, fill = T)
# ss_file <- ss_file[ , c(1:8, 11, 12, 18)]
# ss_file <- ss_file[ss_file$position != -9, ]
# ss_file_hapmap <- merge(hapmap_list, ss_file, by.x = 'rsID', by.y = 'rs_number')
# ss_file <- ss_file_hapmap

#ss_file$SNP_hg18 <- unlist(lapply(strsplit(ss_file$SNP_hg18, "r"), "[", 2))

# colnames(ss_file)[1] <- "rsid"
# colnames(ss_file)[2] <- "rsid"
# colnames(ss_file)[3] <- "rsid"

setDF(ss_file)

if (positions_added_manually == TRUE){
  ss_file$CHR <- unlist(lapply(strsplit(ss_file[,ncol(ss_file)], ":"), "[", 1))
  ss_file$POS <- unlist(lapply(strsplit(ss_file[,ncol(ss_file)-1], ":"), "[", 2))
  ss_file[ ,ncol(ss_file)-2] <- NULL
}

chr_colnum <- 10
pos_colnum <- 11
snp_colnum <- 1
effect_allele_colnum <- 4
other_allele_colnum <- 5
effect_size_colnum <- 6
zscore_colnum <- NA
se_colnum <- 7
p_colnum <- 8
sample_size_colnum <- 9
sample_size <- NA

check_missingness <- function(ss_file){
  cat("Missing entries check: \n")
  
  if (sum(is.na(ss_file)) != 0){
    print(ss_file[which(is.na(ss_file)), ])
    cat("Excluded", nrow(ss_file[which(is.na(ss_file)), ]), "SNPs with NA. \n")
    ss_file_nona <- na.omit(ss_file)
    return(ss_file_nona)
  } 
  else if (sum(ss_file == -9) != 0 || is.na(ss_file)){
    print(ss_file[which(ss_file == -9), ])
    cat("Excluded", nrow(ss_file[which(ss_file == -9), ]), "SNPs with -9. \n")
    ss_file_nona <- ss_file[which(ss_file != -9), ]
    return(ss_file_nona)
  }
  else {
    cat("All good, no missing data.\n")
  }
  
  return(ss_file)
}

remove_positions <- function(ss_file, chr_colnum, pos_colnum){
  ss_file_out <- ss_file
  ss_file_out[ ,pos_colnum] <- NULL
  ss_file_out[ ,chr_colnum] <- NULL
  return(ss_file_out)
}

generate_script <- function(phenotype, input_path, scriptSaveTo, imprangeSaveTo){
  cat("#PBS -lwalltime=24:00:00\n#PBS -lselect=1:ncpus=1:mem=10gb\n#!/bin/bash -x\n\n", 
      file = scriptSaveTo)
  cat("module load gsl\nmodule load anaconda3/personal\n\ni=$PBS_ARRAY_INDEX\n", 
      file = scriptSaveTo,
      append = T)
  cat("input_path=", input_path, sep="", file = scriptSaveTo, append = T)
  cat("\noutput_path=", output_path, sep="", file = scriptSaveTo, append = T)
  cat("\nmkdir $output_path", file = scriptSaveTo, append = T)
  
  cat("\n\nsummary_stats=", gwasSaveTo, sep = "", file = scriptSaveTo, append = T)
  cat("\ndeclare -a imprange\nreadarray -O 1 -t imprange < ", input_path, imprangeSaveTo, sep = '', file = scriptSaveTo, append = T)
  cat("\n\n$WORK/ssimp_software-master/ssimp --gwas ${input_path}${summary_stats} --ref $WORK/reference_panels/1000genomes/ALL.chr{CHRM}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --out res-${i}.txt ${imprange[$i]} --tag.maf 0.01 --impute.maf 0.01",
      file = scriptSaveTo, append = T)
  cat("\n\ncp * ${output_path}", file = scriptSaveTo, append = T)
}

format_columns <- function(ss_file,
                           chr_colnum,
                           pos_colnum,
                           snp_colnum,
                           effect_allele_colnum,
                           other_allele_colnum,
                           effect_size_colnum = NA,
                           se_colnum = NA,
                           p_colnum,
                           z_colnum = NA,
                           sample_size_colnum = NA,
                           sample_size = NA){
  
  if (is.na(sample_size_colnum)){
    ss_file$N <- sample_size
    sample_size_colnum <- ncol(ss_file)
  } else{
    colnames(ss_file)[sample_size_colnum] = "N"
  }
  
  if (is.na(se_colnum) | is.na(effect_size_colnum)){
    num_cols <- c(chr_colnum, pos_colnum, z_colnum, p_colnum, sample_size_colnum)
    colnames(ss_file)[zscore_colnum] <- "ZSCORE"
  } else {
    num_cols <- c(chr_colnum, pos_colnum, effect_size_colnum, se_colnum, p_colnum, sample_size_colnum)
    if (colnames(ss_file)[effect_size_colnum] != "BETA"){
      colnames(ss_file)[effect_size_colnum] = "BETA"
    }
    if (colnames(ss_file)[snp_colnum] != "SNP"){
      colnames(ss_file)[snp_colnum] = "SNP"
    } 
    if (colnames(ss_file)[se_colnum] != "SE"){
      colnames(ss_file)[se_colnum] = "SE"
    } 
  }
  
  colnames(ss_file)[chr_colnum] <- "CHR"
  colnames(ss_file)[pos_colnum] <- "POS"
  
  if (colnames(ss_file)[effect_allele_colnum] != "EFFECT_ALLELE"){
    colnames(ss_file)[effect_allele_colnum] = "EFFECT_ALLELE"
  }
  if (colnames(ss_file)[other_allele_colnum] != "OTHER_ALLELE"){
    colnames(ss_file)[other_allele_colnum] = "OTHER_ALLELE"
  }
  if (colnames(ss_file)[snp_colnum] != "SNP"){
    colnames(ss_file)[snp_colnum] = "SNP"
  } 
  if (colnames(ss_file)[p_colnum] != "P"){
    colnames(ss_file)[p_colnum] = "P"
  } 
  
  ss_file[ , num_cols] <- as.numeric(unlist(ss_file[ , num_cols]))
  
  return(ss_file)
}

remove_othercols <- function(ss_file,
                            chr_colnum,
                            pos_colnum,
                            snp_colnum,
                            effect_allele_colnum,
                            other_allele_colnum,
                            effect_size_colnum = NA,
                            z_colnum = NA,
                            se_colnum = NA,
                            p_colnum,
                            sample_size_colnum=NA){
  
  if (is.na(sample_size_colnum)) {
    sample_size_colnum <- ncol(ss_file1)
  }
  
  if (is.na(se_colnum) | is.na(effect_size_colnum)){
    necessary_cols <- c(chr_colnum,
                        pos_colnum,
                        snp_colnum,
                        effect_allele_colnum,
                        other_allele_colnum,
                        z_colnum,
                        sample_size_colnum)
  } else {
    necessary_cols <- c(chr_colnum,
                        pos_colnum,
                        snp_colnum,
                        effect_allele_colnum,
                        other_allele_colnum,
                        effect_size_colnum,
                        se_colnum,
                        p_colnum,
                        sample_size_colnum)
  }

  ss_file <- ss_file[ , necessary_cols]
  return(ss_file)
}

ss_file1 <- format_columns(ss_file,
                           chr_colnum,
                           pos_colnum,
                           snp_colnum,
                           effect_allele_colnum,
                           other_allele_colnum,
                           effect_size_colnum = effect_size_colnum,
                           se_colnum = se_colnum,
                           p_colnum = p_colnum,
                           z_colnum = NA,
                           sample_size_colnum = sample_size_colnum,
                           sample_size = NA)

ss_file2 <- remove_othercols(ss_file1,
                             chr_colnum,
                             pos_colnum,
                             snp_colnum,
                             effect_allele_colnum,
                             other_allele_colnum,
                             effect_size_colnum = effect_size_colnum,
                             z_colnum = NA,
                             se_colnum = se_colnum,
                             p_colnum = p_colnum,
                             sample_size_colnum=sample_size_colnum)

ss_file3 <- check_missingness(ss_file2)

generate_script(phenotype, input_path, scriptSaveTo, imprangeSaveTo)

ss_file_out <- remove_positions(ss_file3, 1, 2)
write.table(ss_file_out, gwasSaveTo, quote = F, sep = '\t', row.names = F)

i <- 0
count <- 0
imprange_array <- vector()

for (chr in 1:22){
  GwasByChr <- ss_file3[ which(ss_file3$CHR == chr), ]
  NumSnps <- nrow(GwasByChr)
  print(NumSnps)
  count <- count + NumSnps
  
  # we want a size of chunk to be around 10000, but that might give nasty remainders
  Chunk <- 15000
  # adjust the size till we have a better remainder
  while((NumSnps%%Chunk < Chunk/2) == T){
    Chunk <- Chunk + 500
  }
  cat("Final number of snps per job for chr", chr, ":", Chunk, "\n")
  cat("Number of jobs per chr", chr, ":", (NumSnps%/%Chunk)+1, "\n")
  
  GwasByChrSorted <- GwasByChr[order(GwasByChr$POS), ]
  ChunkStart <- 1
  while (ChunkStart < NumSnps){
    ChunkEnd <- ChunkStart + Chunk
    StartPos <- GwasByChrSorted[ChunkStart, c("POS")]
    EndPos <- GwasByChrSorted[ChunkEnd, c("POS")]
    ChunkStart <- ChunkEnd
    i <- i + 1
    
    if (is.na(EndPos) ){
      EndPos <- GwasByChrSorted[NumSnps, c("POS")]
    }
    
    cat("For ", i, ": --impute.range ", chr, ":", StartPos, "-", chr, ":", EndPos, "\n", sep = "")
    
    imputerange <- paste("--impute.range ", chr, ":", StartPos, "-", chr, ":", EndPos, sep = "")
    
    imprange_array[i] <- imputerange
    
    write(imprange_array, imprangeSaveTo, ncolumns = 1)
  }
}

remove(list = ls())


                               