# Dissecting shared predisposition to cancer and diabetes using multi-dimensional analytical approach

This is a repository containing the code I used for my MEng individual project.
I have tried to put most of the documents into folders described below, but this isn't possible for pdf files.

**A full resolution version of the heatmap:**
at_euclcol_euclrow_not1d_truncated_uncorr_v5.pdf
Comment: the key for z-score coulours should be a continuous scale from blue to red through yellow. This must be a pdf upload bug.

**Additional findings in the [cancer loci x metabolome] plot:**
findings.xlsx - description of 6 loci with a significance overlap for several phenotypes

**Code for generating the table and clustering:**
make_locilist.R - combining loci 
functions_ztable.R - file containing functions used by other scripts
make_ztable.R
LD_data.txt
all_phenotype_table_v3.txt - the table used as an auxilary file for producing a table of z-scores
metabolites_phenotype_table.xlsx - the table used as an auxilary file for producing a table of z-scores for metabolites

**Code for preprocessing and postprocessing the files for imputation:**
ssimp-processing-workflow.pdf - the general workflow for imputation
generalised_script.R - the final generalised version of preprocessing script
postprocess_ssimp.sh
reinsert_snps.R
scan_post_imp.R
insert_originalGWAS_into_SSimp.R

**The tables of SNPs x Phenotypes including metabolites:**
v3_wtrepeats_wtmetabolites_wtannotations.txt

**Lists of loci:**
t2d_loci_mahajan.xlsx
Cardiometabolic_loci.xlsx
GWAS_catalog_Cancer_summary_14.02.19.xlsx - the compilation of outputs from GWAS catalog
cancer_loci.xlsx - a later curated version of cancer loci
combined_loci + metabolites loci.xlsx

**Supporting plots:**
Rplot.pdf - phenotype correlation
Rplot01.pdf - contributors to principal components
k3 with labels.pdf - k-means clustering, number of clusters = 3, with labels
k3.pdf - k-means clustering, number of clusters = 3, without labels




