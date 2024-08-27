# sc
import scdrs
import pegasus as pg
import scanpy as sc
import anndata as ad

# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
import os
import csv
import glob
import re
import pynndescent
from scipy.stats import zscore

from synapseclient import Synapse
# Set up the Synapse client
syn = Synapse()
syn.login()  # Assuming you're already logged in or have set up your credentials

# Define the directory and create it if it doesn't exist
data_dir = "scDRS"
os.makedirs(data_dir, exist_ok=True)

# Download the entity to the specified directory
syn62147251 = syn.get(entity="syn62147251", downloadLocation=data_dir)

# Load the synapse intity
adata = sc.read_h5ad(syn62147251)
adata

# Preprocess ann data object via scdrs
scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)

# Load custom geneset
dict_gs = scdrs.util.load_gs('./custom_geneset.gs',
                            src_species="human",dst_species="human",to_intersect=adata.var_names)
dict_gs.keys()

# If you want specific gwas studies only
dict_you_want = dict_gs

# Assign connectivities slot
adata.obsp['connectivities'] = adata.obsp['W_pca_regressed_harmony']

# Calculate per cell/nuclei scDRS core
dict_df_score = dict()
for trait in dict_you_want:
    gene_list, gene_weights = dict_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False)

# Perform downstream analysis - calculate disease association z-score
for trait in dict_df_score:
    df_stats = scdrs.method.downstream_group_analysis(
        adata=adata,
        df_full_score=dict_df_score[trait],
        group_cols=["class"],
    )["class"]
    dict_df_score[trait].to_csv(f"./aging_class_{trait}.csv", sep="	", quoting=csv.QUOTE_NONE, escapechar=' ')
    df_stats.to_csv(f"./aging_class_{trait}_summary.csv", sep="   ", quoting=csv.QUOTE_NONE, escapechar=' ')

for trait in dict_df_score:
    df_stats = scdrs.method.downstream_group_analysis(
        adata=adata,
        df_full_score=dict_df_score[trait],
        group_cols=["subclass"],
    )["subclass"]
    dict_df_score[trait].to_csv(f"./aging_subclass_{trait}.csv", sep="    ", quoting=csv.QUOTE_NONE, escapechar=' ')
    df_stats.to_csv(f"./aging_subclass_{trait}_summary.csv", sep="        ", quoting=csv.QUOTE_NONE, escapechar=' ')
