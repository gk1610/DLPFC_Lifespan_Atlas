from optparse import OptionParser
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.pyplot import rc_context
from harmony import harmonize
import scipy.sparse as sp


sc.settings.verbosity = 1
sc.logging.print_header()
sc.settings.set_figure_params( dpi=200, fontsize=6)
# set number of cores to use
sc.settings.n_jobs = 30

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

def lineage_analysis(lineageData, rand_state, hvgFile, n_pcs, n_neighbors, out_label, batch_key):
	# step0. Load data
	# read in lineage data
	adata = sc.read(lineageData) # 
	print('Data shape: ' + str(adata.shape))
	
	# step1_0. remove low count genes
	sc.pp.filter_genes( adata, min_cells=5)
	print('Data shape after removing low count genes: ' + str(adata.shape))
	adata.write(out_label + '_adata.h5ad')
	print('step1: removing low count genes is finished!')


	# step2. Feature Selection - select 4,000 highly variable genes
	if ((len(hvgFile.split('@')) > 1) &  (hvgFile.split('@')[0] == 'scanpy')):
		n_hvg = int(hvgFile.split('@')[1])
		adata.uns['log1p']["base"] = None
		adata_filtered = adata[:, (adata.var['gene_type-Aging'] == 'protein_coding') & (~adata.var['gene_chrom-Aging'].isin(['X', 'Y', 'MT']))].copy()
		sc.pp.highly_variable_genes( adata_filtered, n_top_genes=n_hvg, n_bins=20, flavor='seurat', inplace=True)
		adata.var['highly_variable'] = adata.var_names.isin(adata_filtered.var_names[adata_filtered.var['highly_variable']])
	else:
		hvgs = list(set(pd.read_csv(hvgFile, header=None)[0]) & set(adata.var_names))
		adata.var['highly_variable'] = False
		adata.var.loc[hvgs, 'highly_variable'] = True
	# ax = sc.pl.highly_variable_genes( adata, show=False)
	# plt.savefig( f"{out_label}_highly_variable_genes.png", format='png', bbox_inches='tight')

	adata.write(out_label + '_adata.h5ad')
	high_adata = adata[:,adata.var.highly_variable.values]
	print('Shape of high adata: ' + str(high_adata.shape))
	print("step2: feature selection is finished!")

	# step3. PCA
	if(n_pcs == 'half'):
		n_comps = min(500, high_adata.shape[1]-1)
		sc.pp.pca( high_adata, svd_solver='arpack', random_state=rand_state, n_comps=n_comps) #, use_highly_variable=True) # add random_stage at Oct-19-2023 11:30 AM
		n_pcs = ( np.cumsum( high_adata.uns['pca']['variance_ratio'])<0.50).sum()
		if( n_pcs==n_comps):
			n_comps = 1000
			sc.pp.pca( high_adata, svd_solver='arpack', random_state=rand_state, n_comps=n_comps) #, use_highly_variable=True) # add random_stage at Oct-19-2023 11:30 AM
			n_pcs = ( np.cumsum( high_adata.uns['pca']['variance_ratio'])<0.50).sum()
			if( n_pcs==n_comps):
				n_pcs = None
				print( "re-run PCA with more components")
	else:
		n_pcs = int(n_pcs)
		sc.pp.pca( high_adata, svd_solver='arpack', random_state=rand_state, n_comps=n_pcs) # add random_stage at Oct-19-2023 11:30 AM

	print('n_pcs used is: ' + str(n_pcs))
	ax = sc.pl.pca_variance_ratio( high_adata, show=False)
	plt.savefig(f"{out_label}_pca_variance_ratio.png", format='png', bbox_inches='tight')
	ax = sc.pl.pca( high_adata, color=['stage_id'], components=['1,2','3,2'], legend_fontsize=10, wspace=0.25, show=False)
	plt.savefig(f"{out_label}_pca_stage_id.png", format='png', bbox_inches='tight')
	high_adata.write(out_label + '_high_adata.h5ad')
	print('step3: PCA calculation is finished!')


	# step4. Remove batch effects
	if(batch_key != "no"):
		high_adata.obsm["X_pca_harmony"] = harmonize(high_adata.obsm["X_pca"], high_adata.obs, batch_key=batch_key, max_iter_harmony=20)
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step4: Batch effects removing is finished!')

		# step5. UMAT
		high_adata = umat( adata.uns['stage_order'], high_adata, pcaType='X_pca_harmony', n_neighbors=n_neighbors, n_pcs=n_pcs, rand=rand_state)
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step5: UMAT calculation is finished!')

		# step6. UMAP
		# step6_0. find nearest neighbors
		sc.pp.neighbors( high_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=rand_state, use_rep="X_pca_harmony")
		# step6_1. calc umap
		sc.tl.umap( high_adata, random_state=rand_state)
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step6: UMAP calculation is finished!')

	else:
		# step5. UMAT
		high_adata = umat( adata.uns['stage_order'], high_adata, pcaType='X_pca', n_neighbors=n_neighbors, n_pcs=n_pcs, rand=rand_state)
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step5: UMAT calculation is finished!')

		# step6. UMAP
		# step6_0. find nearest neighbors
		sc.pp.neighbors( high_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=rand_state, use_rep="X_pca")
		# step6_1. calc umap
		sc.tl.umap( high_adata, random_state=rand_state)
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step6: UMAP calculation is finished!')

	# Plot UMAT
	# swap order of x-axis to match umap, just makes visual comparison easier
	# high_adata.obsm['X_umat'][:,0] = high_adata.obsm['X_umat'][:,0] * -1
	ax = sc.pl.embedding( high_adata, basis='umat', color=['cell_type_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_cell_type_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subclass_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subtype_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['Batch'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_Batch.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['SubID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_SubID.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['poolID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_poolID.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['prep'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_prep.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_stage_id.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['class'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_class.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subclass.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subtype.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['cell_type'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_cell_type.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['Sex'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_Sex.png", format='png', bbox_inches='tight')

	# Plot UMAP
	ax = sc.pl.umap( high_adata, color=['cell_type_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_cell_type_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subclass_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subclass_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subtype_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subtype_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['Batch'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_batch.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['SubID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_SubID.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['poolID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_poolID.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['prep'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_prep.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_stage_id.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['class'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_class.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subclass'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subclass.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subtype'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subtype.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['cell_type'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_cell_type.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['Sex'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_Sex.png", format='png', bbox_inches='tight')

	# step7. Update full count matrix
	# set umap and pca to adata
	adata.obsm = high_adata.obsm
	adata.write(out_label + '_clustering.h5ad')

def main():
	parser = OptionParser()
	parser.add_option('-l', '--lineageData', dest='lineageData', help='Lineage data.')
	parser.add_option('-r', '--rand_state', dest='rand_state', help='Random state.')
	parser.add_option('-v', '--hvgFile', dest='hvgFile', help='Highly variable gene lists.')
	parser.add_option('-p', '--numPCs', dest='numPCs', help='Number of PCs.')
	parser.add_option('-n', '--numNeighbors', dest='numNeighbors', help='Number of neighbors.')
	parser.add_option('-o', '--out_label', dest='out_label', help='Output label.')
	parser.add_option('-k', '--batch_key', dest='batch_key', help='Batch or not?')

	(options,args) = parser.parse_args()
	lineage_analysis(options.lineageData, int(options.rand_state), options.hvgFile, options.numPCs, int(options.numNeighbors), options.out_label, options.batch_key)

main()
