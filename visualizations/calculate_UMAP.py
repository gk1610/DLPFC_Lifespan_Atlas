# sc
import pegasus as pg
import scanpy as sc
import scanpy.external as sce
import anndata as ad
from anndata.tests.helpers import assert_equal
from anndata._core.sparse_dataset import SparseDataset
from anndata.experimental import read_elem, write_elem

# plotting
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns

# data
import numpy as np
import pandas as pd
from scipy import stats
from scipy import sparse
import h5py

# sys
import gc
from pathlib import Path

# harmony
from harmony import harmonize

# pge
import sys
import pge
pge.info()

## args
orig_h5ad = './Aging_freeze3.0_full.h5ad'
new_h5ad = './freeze3_AGING_clean_lcg_autosome_100perc.h5ad'
chunk_size = 500000

### load data
adata = pge.read_everything_but_X(orig_h5ad)

### get obs subset
subset_obs = (adata.obs_names!=None).tolist()

### get var subset
subset_var = ((adata.var.gene_type=='protein_coding')&(~adata.var.gene_chrom.isin(['MT', 'X', 'Y']))).tolist()

### run ondisk_subset
pge.ondisk_subset(orig_h5ad=orig_h5ad, new_h5ad=new_h5ad, subset_obs=subset_obs, subset_var=subset_var,
                  chunk_size=chunk_size, raw=True)

# check
with h5py.File(new_h5ad, 'r') as f:
    f["X"].visititems(print)


### args
args_input = './freeze3_AGING_clean_lcg_autosome_100perc.h5ad'
args_flavor = 'cell_ranger'
args_batch = 'Source'
args_n_top_genes = 6000 # was None

### HVF
hvg = pge.scanpy_hvf_h5ad(h5ad_file=args_input, flavor=args_flavor, batch_key=args_batch,
                          n_top_genes=args_n_top_genes, min_mean=0.0125, max_mean=3,
                          min_disp=0.5, protein_coding=True, autosome=True)
print(len(hvg))

# open file in write mode
with open('AGING_hvg6k.txt', 'w') as fp:
    for x in hvg:
        fp.write("%s\n" % x)
    print('Done')

with open('AGING_hvg6k.txt') as f:
    hvg = [line.rstrip() for line in f]
len(hvg)

### args
args_input = './freeze3_AGING_clean_lcg_autosome_100perc.h5ad'
args_output = './freeze3_AGING_clean_lcg_autosome_100perc_umap.h5ad'

args_batch = 'Source'
args_n_pcs = 30
args_n_neighbors = 100

### load data
data = pg.read_input(args_input, genome='GRCh38', modality='rna')
print(data)

### set scanpy hvg as hvf
data.var.highly_variable_features = False
data.var.loc[data.var.index.isin(hvg),'highly_variable_features'] = True

### final value counts
print(data.var.highly_variable_features.value_counts())
print(data.var[data.var.highly_variable_features==True].gene_chrom.value_counts())

### pca
pg.pca(data, n_components=args_n_pcs)
pg.elbowplot(data)
npc = min(data.uns["pca_ncomps"], args_n_pcs)
print('Using %i components for PCA' % npc)


from harmony import harmonize

### harmony/kNN Donghoon's setup
pg.regress_out(data, attrs=['n_counts','percent_mito','cycle_diff'])
pg.run_harmony(data, batch=args_batch, rep='pca_regressed', max_iter_harmony=20, n_comps=npc)
pg.neighbors(data, rep='pca_regressed_harmony', use_cache=False, dist='l2', K=100, n_comps=npc)

### umap
pg.umap(data, rep='pca_regressed_harmony', n_neighbors=args_n_neighbors, rep_ncomps=npc)

### save
pge.save(data, args_output)

