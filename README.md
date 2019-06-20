# Dissecting shared predisposition to cancer and diabetes using multi-dimensional analytical approach <br />

This is a repository containing the code I used for my MEng individual project. <br />
I have sorted most of the documents into folders described below, but this isn't possible for pdf files. <br />

**A full resolution version of the heatmap:** <br />
at_euclcol_euclrow_not1d_truncated_uncorr_v5.pdf <br />
Comment: the key for z-score coulours should be a continuous scale from blue to red through yellow. This must be a pdf upload bug. <br />

**Additional findings in the [cancer loci x metabolome] plot:** <br />
findings.xlsx - description of 6 loci with a significance overlap for several phenotypes <br />

**Code for generating the table and clustering:** <br />
make_locilist.R - combining loci <br />
functions_ztable.R - file containing functions used by other scripts <br />
make_ztable.R <br />
LD_data.txt <br />
all_phenotype_table_v3.txt - the table used as an auxilary file for producing a table of z-scores <br />
metabolites_phenotype_table.xlsx - the table used as an auxilary file for producing a table of z-scores for metabolites <br />

**Code for preprocessing and postprocessing the files for imputation:** <br />
ssimp-processing-workflow.pdf - the general workflow for imputation <br />
generalised_script.R - the final generalised version of preprocessing script <br />
postprocess_ssimp.sh <br />
reinsert_snps.R <br />
scan_post_imp.R <br />
insert_originalGWAS_into_SSimp.R <br />

**The tables of SNPs x Phenotypes including metabolites:** <br />
v3_wtrepeats_wtmetabolites_wtannotations.txt <br />

**Lists of loci:** <br />
t2d_loci_mahajan.xlsx <br />
Cardiometabolic_loci.xlsx <br />
GWAS_catalog_Cancer_summary_14.02.19.xlsx - the compilation of outputs from GWAS catalog <br />
cancer_loci.xlsx - a later curated version of cancer loci <br />
combined_loci + metabolites loci.xlsx <br />

**Supporting plots:** <br />
Rplot.pdf - phenotype correlation <br />
Rplot01.pdf - contributors to principal components <br />
k3 with labels.pdf - k-means clustering, number of clusters = 3, with labels <br />
k3.pdf - k-means clustering, number of clusters = 3, without labels <br />




