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
from joblib import dump

import multiprocessing


def score_trait(adata, trait, gene_list, gene_weights):
    score = scdrs.score_cell(
        data=adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=200,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
        verbose=False,
    )
    return trait, score


def scDRS_calculate(file_h5ad, file_gs):
    adata = sc.read_h5ad(file_h5ad)

    scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)

    dict_gs = scdrs.util.load_gs(file_gs,
                                 src_species="human", dst_species="human", to_intersect=adata.var_names)

    dict_df_score = dict()

    with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
        results = [
            pool.apply_async(score_trait, (adata, trait, dict_gs[trait][0], dict_gs[trait][1]))
            for trait in dict_gs
        ]

        for result in results:
            trait, score = result.get()
            dict_df_score[trait] = score

    return dict_df_score


def scDRS_calculate_simple(file_h5ad, file_gs):
    adata = sc.read_h5ad(file_h5ad)

    scdrs.preprocess(adata, n_mean_bin=20, n_var_bin=20, copy=False)

    dict_gs = scdrs.util.load_gs(file_gs,
                            src_species="human", dst_species="human", to_intersect=adata.var_names)


    dict_df_score = dict()
    for trait in dict_gs:
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
            verbose=False,
        )
    return(dict_df_score)