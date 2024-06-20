from optparse import OptionParser
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import rc_context
import seaborn as sns
import re
import scrublet as scr
import math
from scipy import stats
import pickle
import random
import scipy.sparse as sp
import fnmatch
from harmony import harmonize


sc.settings.verbosity = 1
sc.logging.print_header()
# set number of cores to use
sc.settings.n_jobs = 25
sc.settings.set_figure_params( dpi=300, fontsize=6)


##### UMAT #####
def umat_neighbors(stage_order, adata, pcaType='X_pca_harmony', n_neighbors=5, n_pcs=50, rand=123):
	# Use sparse matrix to save memory
    # if uns['neighbors'] doesnt exist run neighbors on whole lot
    if not('neighbors' in adata.uns):
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=rand, use_rep=pcaType)

    n_cells = adata.shape[0]

    # Initialize as sparse matrices
    conn = sp.lil_matrix((n_cells, n_cells))
    dist = sp.lil_matrix((n_cells, n_cells))
    inds = np.arange( n_cells)

    # neighbors iteratively for adjacent stages
    for itr1, itr2 in zip(stage_order[:-1], stage_order[1:]):
        itr_mk = np.in1d(adata.obs['stage_id'], [itr1, itr2])
        itr_adata = adata[itr_mk]

        sc.pp.neighbors(itr_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=rand, use_rep=pcaType)
        
        # No need to convert to dense matrix
        con_itr = itr_adata.obsp['connectivities']
        dis_itr = itr_adata.obsp['distances']
        
        for ii, (c_itr, d_itr) in zip( inds[itr_mk], zip( con_itr, dis_itr)):
            conn[ii,itr_mk] = c_itr
            dist[ii,itr_mk] = d_itr
        print( itr1, itr2, itr_mk.sum())

    adata.obsp['connectivities'] = conn.tocsr()
    adata.obsp['distances'] = dist.tocsr()

    return adata


def umat( stage_order, adata, pcaType='X_pca_harmony', n_neighbors=5, n_pcs=50, min_dist=0.5, n_comps=2, rand=123):
	adata = umat_neighbors( stage_order, adata, pcaType=pcaType, n_neighbors=n_neighbors, n_pcs=n_pcs, rand=rand)
	if ('X_umap' in adata.obsm):
		umap = adata.obsm['X_umap']
	# need to update function to include variables for umap
	sc.tl.umap( adata, random_state=rand, min_dist=min_dist, n_components=n_comps)
	adata.obsm['X_umat'] = adata.obsm['X_umap']
	if ('X_umap' in adata.obsm):
		adata.obsm['X_umap'] = umap
	else:
		del adata.obsm['X_umap']
	return( adata)


def umat_calculation(input_files, pca_type, n_pcs, n_neighbors, rand_state, out_label):
	adata = sc.read(input_files)
	high_adata = adata[:,adata.var.highly_variable.values]

	high_adata = umat( high_adata.uns['stage_order'], high_adata, pcaType=pca_type, n_pcs=int(n_pcs), n_neighbors=n_neighbors, rand=rand_state) #123, #12 #12345 #123456
	# set updated obsm to adata
	adata.obsm = high_adata.obsm
	adata.write(out_label + '_UMAT_adata.h5ad')
	print('UMAT calculation is finished!')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['Batch'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_Batch.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['SubID'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_SubID.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['poolID'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_poolID.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['prep'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_prep.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_stage_id.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['class'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_class.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subclass.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['cell_type_uni'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_cell_type_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass_uni'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subclass_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype_uni'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subtype_uni.png", format='png', bbox_inches='tight')


def main():
	parser = OptionParser()
	parser.add_option('-i', '--input_files', dest='input_files', help='Input files.')
	parser.add_option('-t', '--pca_type', dest = 'pca_type', help = 'Type of PCA.')
	parser.add_option('-p', '--n_pcs', dest='n_pcs', help='Number of PCs.')
	parser.add_option('-n', '--n_neighbors', dest='n_neighbors', help='Number of neighbors.')
	parser.add_option('-r', '--rand_state', dest='rand_state', help='Random state.')
	parser.add_option('-o', '--out_label', dest='out_label', help='Output file.')

	(options,args) = parser.parse_args()
	print(options)
	umat_calculation(options.input_files, options.pca_type, options.n_pcs, int(options.n_neighbors), int(options.rand_state), options.out_label)

main()