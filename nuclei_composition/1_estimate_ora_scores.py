import sys
import pandas as pd
import scanpy as sc
import decoupler as dc
import numpy as np

# Read Brain related go pathways and corresponding genes list
# get file from syn72106811 (split this file into chunks for efficiency)
# syn72106811 are GO pathways that contsins following brain_keywords = c(
#  "brain", "neuron", "neuro", "glia", "synapse", "axon", 
#  "dendrite", "astrocyte", "myelin", "cns", "nerve", 
#  "cortex", "prefrontal", "dlpfc", "pfc", "frontal", 
#  "microglia", "oligodendrocyte", "excitatory", "inhibitory")

go_bp_df = pd.read_csv(f'brain_related_GO_pathways.csv')

# Read AnnData from syn62147251
adata = sc.read_h5ad('aging_umap40pca_age_synapse.h5ad')

# Run ORA for chunk
dc.run_ora(
    mat=adata,
    net=go_bp_df,
    source='geneset',
    target='genesymbol',
    verbose=True,
    use_raw=False
)

# Save results
acts = dc.get_acts(adata, obsm_key='ora_estimate')
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e
pd.DataFrame(acts.obsm['ora_estimate']).to_csv(f'ORA_brain_related_GO_pathways.csv') 

