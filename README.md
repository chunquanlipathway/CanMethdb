
# CanMethdb: a database for genome-wide DNA methylation annotation in cancers

## Table of contents AND pipline

  - [NearGenes](# Get candidate genes for each CpG)
     1. NearGenes.R
     2. run_NearGene.R
  - [Common samples](# Get DNA methylation and gene expression profiles of common samples)
	 1. profiles_com_sample.py
  - [Data pre-processing]
     1. gene_profile_pre-processing.py
	 2. methylation_profile_pre-processing.py
  - [ELMER](# Use ElMER to calculate target genes for each CpG in each cancer type)
     1. EMLMER.R
     2. run_ELMER.R
  - [Pearson](Use Pearson’s correlation coefficient method to calculate target genes for each CpG in each cancer type )
     1. pearson.py
  - [limma](# Using limma to find  differentially methylated CpGs)
     1. limma.R
 
***********
