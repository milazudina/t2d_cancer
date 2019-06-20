install.packages("readxl")
library(plyr)
library(readxl)
install.packages("data.table")
library(data.table)
install.packages("stringr")
library(stringr)

# Import all the loci tables from Zhanna
cardiometabolic_loci <- readxl::read_xlsx("Cardiometabolic_6.12.18.xlsx")
depression_loci <- readxl::read_xls("depression_snps.xls")
crc_loci <- readxl::read_xlsx("GWAS_catalog_Cancer_summary_14.02.19.xlsx", sheet = 1)
proc_loci <- readxl::read_xlsx("GWAS_catalog_Cancer_summary_14.02.19.xlsx", sheet = 2)
brc_loci <- readxl::read_xlsx("GWAS_catalog_Cancer_summary_14.02.19.xlsx", sheet = 3)
panc_loci <- readxl::read_xlsx("GWAS_catalog_Cancer_summary_14.02.19.xlsx", sheet = 4)

# Clean up: leave only 4 columns (chr.pos, rsid, phenotype, locus)
normalise_locitables <- function(input, chrpos_col, rsid_col, phenotype_col, locus_col) {
  new_name <- paste(input, 'reduced', sep = '_')
  colnames(input)[chrpos_col] <- 'chr.pos'
  colnames(input)[phenotype_col] <- 'phenotype'
  colnames(input)[locus_col] <- 'locus'
  input_loci_reduced <- input[, c(chrpos_col, rsid_col, phenotype_col, locus_col)]
  return(input_loci_reduced)
}

# function to add all LD details
addLD <- function(data) {
  data$LD <- word(data$Details, 5)
  data$LD <- ifelse(grepl("rs", data$LD), data$LD, NA)
  allLociLD <- merge(allLociLD, data, by.x = 'rsid', by.y = 'RS.Number', all.x = T)
  allLociLD$LD[is.na(allLociLD$LD)] <- allLociLD$rsid[is.na(allLociLD$LD)]
  return(allLociLD)
}

cardiometabolic_loci_reduced <- normalise_locitables(cardiometabolic_loci, 3, 4, 5, 6)
proc_loci_reduced <- normalise_locitables(proc_loci, 3, 4, 5, 6)
panc_loci_reduced <- normalise_locitables(panc_loci, 3, 4, 5, 6)
crc_loci_reduced <- normalise_locitables(crc_loci, 3, 4, 5, 6)
brc_loci_reduced <- normalise_locitables(brc_loci, 3, 4, 5, 6)

depression_loci$Markername <- paste(depression_loci$Chromosome, depression_loci$`Postion (bp)`, sep = ":")
depression_loci_reduced <- depression_loci[, c(2, 46)]
colnames(depression_loci_reduced)[1] <- 'rsid'
colnames(depression_loci_reduced)[2] <- 'chr.pos'
depression_loci_reduced$phenotype <- 'Depression'
depression_loci_reduced <- depression_loci_reduced[, c(2, 1, 3)]

remove(cardiometabolic_loci, proc_loci, panc_loci, crc_loci, brc_loci, depression_loci)

# Concatenate all the tables
locitables_list <- list(cardiometabolic_loci_reduced, depression_loci_reduced, crc_loci_reduced, proc_loci_reduced, brc_loci_reduced, panc_loci_reduced)
combined_loci <- ldply(locitables_list, rbind)
remove(cardiometabolic_loci_reduced, proc_loci_reduced, panc_loci_reduced, crc_loci_reduced, brc_loci_reduced, depression_loci_reduced)
remove(locitables_list)
# Tidy up & check for NA
combined_loci$phenotype <- toupper(combined_loci$phenotype)
#phenotype_counts <- count(combined_loci$phenotype)

# remove entries where both chr.pos and rsid are NA
# removed rs6762208 which was not found in GWAS catalog and had NA for phenotype
# chr.pos == :
# deleted 2 lines woith multiple snps and NAs
combined_loci <- combined_loci[which(combined_loci$chr.pos != 'NA' & combined_loci$rsid != 'NA'), ]
combined_loci <- combined_loci[combined_loci$rsid != 'HLA-DRB1*07:01', ]


# Import the pre-made in shell list of snps with positions in build 37 - see the file commands-smth-smth
# This step is crucial because some of the positions from breast cancer and other new ss are from build 38

rsid37_build <- read.table('rsid_37build.txt', header = FALSE, sep = " ", dec = ".", na.strings = "NA", stringsAsFactors = FALSE, fill = TRUE)

# At this point, start doing the merge with positions for db141
# Done in bash #
# write.table(combined_loci, file=paste("rsid_list.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
loci_db141 <- merge(combined_loci, rsid37_build, by.x = 'rsid', by.y = 'V1', all.x = TRUE)
#`so `chr.pos in db142, and V2 is db141
rm(combined_loci, rsid37_build)
#################### from the untitled ##################

loci_db141 <- loci_db141[!duplicated(loci_db141$rsid), ]

# 1) Separate Chr:Pos to split the snps by chromosome number
loci_db141$chr <- unlist(lapply(strsplit(loci_db141$chr.pos, ":"), "[", 1))

# THIS STEP ONLY IF YOU NEE TO REDO THE LD IN SNPCLIP OR SO
# 2) Split into the list and write into files
chrSeparatedLoci <- split(loci_db141, loci_db141$chr)
lapply(names(chrSeparatedLoci),
       function(x, chrSeparatedLoci) write.table(chrSeparatedLoci[[x]], paste(x, ".txt", sep = ""),
                                                 col.names=FALSE, row.names=FALSE, sep="\t", 
                                                 quote=FALSE), chrSeparatedLoci)

# 3) Apply SNPClip and return the files here
chrSeparatedLoci[[1]] <- NULL
chrSeparatedLoci[[16]] <- NULL
allLociLD <- ldply(chrSeparatedLoci, rbind)
rm(chrSeparatedLoci, loci_db141)
# allLociLD <- data.frame()
# col.names <- c("rsid", "chr.pos", "phenotype", "locus", "V2", "chr", "Position", "Alleles", "Details", "LD")
# allLociLD <- rbind(allLociLD, c(rep(NA, length(col.names))))
# colnames(allLociLD) <- col.names
data <- list()
for (i in 1:23) {
  filename <- paste("LDclumped/", i, "_LD.txt", sep = '')
  print(filename)
  data[[i]] <- read.table(filename, sep = "\t", dec = ".", na.strings = "NA", stringsAsFactors = F, fill = T, header = T)
}

allLD <- ldply(data, rbind)
allLociLD <- addLD(allLD)
rm(allLD, data)

# to exclude the indels, delete all with symbol '-' present in V3
allLociLD$Alleles <- ifelse(grepl("-", allLociLD$Alleles), NA, allLociLD$Alleles)
allLociLD$Details <- ifelse(grepl("biallelic", allLociLD$Details), NA, allLociLD$Details)
allLociLD$Details <- ifelse(grepl("MAF", allLociLD$Details), NA, allLociLD$Details)
allLociLD$Details <- ifelse(grepl("1000G", allLociLD$Details), NA, allLociLD$Details)

allLociLD <- allLociLD[!is.na(allLociLD$Alleles) & !is.na(allLociLD$Details), ]
allLociLD <- allLociLD[!duplicated(allLociLD$rsid), ]


# OPTIONAL, A BRANCH OFF
dat1 <- data.frame(do.call(rbind, strsplit(as.vector(allLociLD$Alleles), split = ", ")))
dat11 <- data.frame(do.call(rbind, strsplit(as.vector(dat1$X1), split = "=")))
dat12 <- data.frame(do.call(rbind, strsplit(as.vector(dat1$X2), split = "=")))
alleles <- cbind(dat11, dat12)
allLociLDAlleles <- cbind(allLociLD, alleles)
allLociLD[, c('rsid')] <- allLociLD[, c('rsid')]
final_loci_list <- allLociLD[, c('rsid', 'chr.pos', 'phenotype', 'locus', 'LD')]

# DON'T DO SHADY THINGS LIKE THIS
# check before doing this
combined_loci <- combined_loci[ -646, ]
combined_loci <- combined_loci[ -c(3985, 3918, 3913), ]


# remove repeats
rsid_count <- count(combined_loci$rsid)
combined_loci <- merge(combined_loci, rsid_count, by.x = 'rsid', by.y = 'x', all.x = TRUE)
colnames(combined_loci)[5] <- 'original_freq'
combined_loci <- unique(combined_loci)
rsid_count <- count(combined_loci$rsid)
combined_loci <- merge(combined_loci, rsid_count, by.x = 'rsid', by.y = 'x', all.x = TRUE)
colnames(combined_loci)[6] <- 'new_freq'



loci_db141 <- loci_db141[!duplicated(loci_db141$rsid), ]
final_loci_list$chr.pos <- gsub('chr', '', final_loci_list$chr.pos)
final_loci_list$chr.pos <- gsub('X', '23', final_loci_list$chr.pos)

final_loci_list_37_unique <- unique(final_loci_list_37)

estLoci <- allLociLD
rm(allLociLD)
for (i in 1:nrow(estLoci)){
  if (!is.na(estLoci[i, c('V2')])) {
    estLoci[i, c('chr.pos')] <- estLoci[i, c('V2')]
  }
}

estLoci[, c(".id", "V2", "chr", "Position", "Alleles", "Details")] <- NULL
# change this to "delete all that != "rsid", "chr.pos", "phenotype", "locus", "LD"

# optional step - see what you can do by hands
count = 0
for (i in grep(':', loci_db141$rsid, fixed = TRUE)) {
  loci_db141[i, 1] = loci_db141[i, 2]
  count = count + 1
}

write.table(final_loci_list_37_unique, 'established_loci_combined.txt', sep = "\t", quote = FALSE)


