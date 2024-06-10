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


adata = sc.read_h5ad("./Aging_freeze3.0_50perc_std1_000.h5ad")
adata

scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)

dict_gs = scdrs.util.load_gs('./all_ms_geneset_newOrd_agingsubset.gs',
                            src_species="human",dst_species="human",to_intersect=adata.var_names)
dict_gs.keys()

# If you want specific gwas studies only
dict_you_want = dict_gs

# assign connectivities slot
adata.obsp['connectivities'] = adata.obsp['W_pca_regressed_harmony']

dict_df_score = dict()
for trait in dict_you_want:
    gene_list, gene_weights = dict_gs[trait]
    dict_df_score[trait] = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=80,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )

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
