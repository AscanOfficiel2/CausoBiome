# -----------------------------------------------
# Script: 03_subset_pathogenic_samples.py
# Authors: Suleiman & AbdulAziz
# Purpose: Extract pathogenic samples with their metadata ( Samples with in ARG/ VFS) before batch correction 
# -----------------------------------------------
#Important: First obtain a list of the ARGs samples and VF samples from the CARD and VFDB file before proceeding

import pandas as pd

# Paths
matrix_path       = "General_species_matrix.csv"
sample_list_path  = "vf_or_arg_samples"    #list of ARG/VF sample IDs from the VFDB and CARD files
meta_path         = "crc_meta.csv"
out_matrix_path   = "Species_subset_matrix.csv"   # subset matrix file of the samples

# 1) Load sample list
samples = pd.read_csv(sample_list_path, header=None, names=["Sample_ID"])
sample_ids = samples["Sample_ID"].astype(str).tolist()

# 2) Load and subset the matrix
mat = pd.read_csv(matrix_path, index_col=0)
common = [s for s in sample_ids if s in mat.columns]
mat_sub = mat[common]
mat_sub.to_csv(out_matrix_path)


