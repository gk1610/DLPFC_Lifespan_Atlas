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
import scipy as sp
import fnmatch
from harmony import harmonize


sc.settings.verbosity = 1
sc.logging.print_header()
# set number of cores to use
sc.settings.n_jobs = 30
sc.settings.set_figure_params( dpi=200, fontsize=6)



###### Numerical age for fetal samples are negative based on 40 week gestation, with -280/365 being conception and 0 being birth ######
def calc_fetal( ga):
	return( ( ( ga * 7) - ( 40 * 7)) / 365)

# convert all age to a common numerical scale and add feature, i.e. years
def numerical_age( adata):
	# get numerical age
	str_age = list(np.unique(adata.obs['Age']))
	num_age = []
	for age_itr in str_age:
		age_itr = str(age_itr)
		digit = float(re.findall('\d+', age_itr)[0])
		if( "ga" in age_itr):
			num_age.append( calc_fetal( digit))
		elif( "d" in age_itr):
			num_age.append( digit / 365)
		elif( "yr" in age_itr):
			num_age.append( digit)
		else:
			num_age.append( digit)
	# add feature to anndata
	adata.obs['numerical_age'] = [num_age[str_age.index(ii)] for ii in adata.obs['Age']]
	return


###### stage of maturation, 
# ref https://www.nature.com/articles/nature10523/tables/1 ######
def stage_id( adata):
	numer_age = adata.obs['numerical_age'].values
	# array to hold stage ids
	stage_id  = np.zeros_like( adata.obs['Date-of-Collection'].values.astype(str))
	stage_id[:] = '!INVALID AGE!'
	
	# Fetal - up to ga38
	stage_id[numer_age<calc_fetal(38)] = 'Fetal'
	# Neonatal - ga38 up to 60d
	stage_id[(numer_age>=calc_fetal(38)) & (numer_age<60/365)] = 'Neonatal'
	# Infancy - 60d up to 1yr 
	stage_id[(numer_age>=60/365) & (numer_age<1)] = 'Infancy'
	# Childhood - 1yr up to 10yr
	stage_id[(numer_age>=1) & (numer_age<12)] = 'Childhood'
	# Adolescence - 10yr up to 20yr
	stage_id[(numer_age>=12) & (numer_age<20)] = 'Adolescen'
	# Adulthood - 20yr up to 40
	stage_id[(numer_age>=20) & (numer_age<40)] = 'Adult'
	# Late Adulthood - 40yr up to 60
	stage_id[(numer_age>=40) & (numer_age<60)] = 'Late_Adult'
	# Old - 60yr and up
	stage_id[numer_age>=60] = 'Old'

	adata.obs['stage_id'] = stage_id
	return


###### Manage features when intergration ######
def merge_data(input_files, batch_labels):
	adata_list = []
	for fl_itr in input_files.split('@'):
		adata_itr = sc.read(fl_itr)
		adata_list.append(adata_itr)
	adata = sc.AnnData.concatenate( *adata_list, join='inner', batch_key='Batch', batch_categories = batch_labels.split('@'))
	return(adata)

def manage_obs(adata):
	# cell_type_uni
	mc_en = ['L2-3_CUX2', 'L4_RORB', 'L5-6_TLE4', 'L5-6_THEMIS', 'PN_dev']
	mc_in = ['VIP', 'ID2', 'LAMP5_NOS1', 'CGE_dev', 'SST', 'PV', 'PV_SCUBE3', 'MGE_dev']
	mc_astro = ['Astro']
	mc_opc = ['OPC']
	mc_oligo = ['Oligo']
	mc_immune = ['Micro']
	mc_others = ['Vas', 'Poor-Quality']
	adata.obs.loc[adata.obs['major_clust'].isin(mc_en), 'cell_type_uni'] = 'EN'
	adata.obs.loc[adata.obs['major_clust'].isin(mc_in), 'cell_type_uni'] = 'IN'
	adata.obs.loc[adata.obs['major_clust'].isin(mc_astro), 'cell_type_uni'] = 'Astro'
	adata.obs.loc[adata.obs['major_clust'].isin(mc_opc), 'cell_type_uni'] = 'OPC'
	adata.obs.loc[adata.obs['major_clust'].isin(mc_oligo), 'cell_type_uni'] = 'Oligo'
	adata.obs.loc[adata.obs['major_clust'].isin(mc_immune), 'cell_type_uni'] = 'Immune'
	adata.obs.loc[adata.obs['major_clust'].isin(mc_others), 'cell_type_uni'] = 'Others'
	adata.obs.loc[adata.obs['Batch'] == 'Aging', 'cell_type_uni'] = adata[adata.obs['Batch'] == 'Aging', :].obs['class']

	# dev_mat_uni
	adata.obs.loc[adata.obs['Batch'] == 'Lister', 'dev_mat_uni'] = adata[adata.obs['Batch'] == 'Lister', :].obs['mat/dev']
	adata.obs.loc[adata.obs['Batch'] == 'Aging', 'dev_mat_uni'] = 'mat'

	# subclass_uni
	adata.obs.loc[adata.obs['Batch'] == 'Lister', 'subclass_uni'] = adata[adata.obs['Batch'] == 'Lister', :].obs['major_clust'].values.astype(str)
	adata.obs['subclass_uni'] = adata.obs['subclass_uni'].astype('str')
	adata.obs.loc[adata.obs['Batch'] == 'Aging', 'subclass_uni'] = adata[adata.obs['Batch'] == 'Aging', :].obs['subclass'].values.astype(str)
	adata.obs['subclass_uni'] = adata.obs['subclass_uni'].astype('category')

	# subtype_uni
	adata.obs.loc[adata.obs['Batch'] == 'Lister', 'subtype_uni'] = adata[adata.obs['Batch'] == 'Lister', :].obs['sub_clust'].values.astype(str)
	adata.obs['subtype_uni'] = adata.obs['subtype_uni'].astype('str')
	adata.obs.loc[adata.obs['Batch'] == 'Aging', 'subtype_uni'] = adata[adata.obs['Batch'] == 'Aging', :].obs['subtype'].values.astype(str)
	adata.obs['subtype_uni'] = adata.obs['subtype_uni'].astype('category')
	
	# numerical_age
	adata.obs['Age'] = adata.obs['Age'].astype(str)
	numerical_age( adata)

	# stage_id
	stage_id(adata)
	adata.obs['stage_id'] = adata.obs['stage_id'].astype('category')
	# set order and colors of stages
	stage_order_full = ['Fetal', 'Neonatal', 'Infancy', 'Childhood', 'Adolescen', 'Adult', 'Late_Adult', 'Old']
	adata.uns['stage_order'] =  [item for item in stage_order_full if item in np.unique(adata.obs['stage_id'])]
	# add colors and color dict to be used for plotting
	state_ordered_colors_full = ["#512568", "#4436F0", "#3D6B93", "#20988C", "#98CA43", "#F9E51B", "#B78600", "#CB0000"]
	stage_ordered_colors = [state_ordered_colors_full[pos] for pos in [stage_order_full.index(item) for item in adata.uns['stage_order'] if item in stage_order_full]]
	adata.uns['stage_colors_dict'] = dict( zip( adata.uns['stage_order'], stage_ordered_colors))
	adata.uns['stage_id_colors'] = [adata.uns['stage_colors_dict'][ii] for ii in adata.obs['stage_id'].cat.categories]
	
	# Sex
	adata.obs.loc[adata.obs['Sex'] == 'F', 'Sex'] = 'Female'
	adata.obs.loc[adata.obs['Sex'] == 'M', 'Sex'] = 'Male'
	
	# SubID
	adata.obs['SubID'] = adata.obs['SubID'].astype(str)
	adata.obs.loc[adata.obs['Batch'] == 'Lister', 'SubID'] = adata[adata.obs['Batch'] == 'Lister', :].obs['RL#'].values.astype(str)
	adata.obs['SubID'] = adata.obs['SubID'].astype('category')
	adata.uns['SubID_order'] = np.unique(adata.obs['SubID'])
	
	# poolID, prep, rep (not sure)
	# features should maintain
	adata.obs = adata.obs[['Batch', 'SubID', 'Age', 'Sex', 'PMI', 'cell_type_uni', 'subclass_uni', 'subtype_uni', 'dev_mat_uni', 'numerical_age', 'stage_id', 'Channel', 'Source', 'poolID_new', 'prep', 'rep', 'class', 'subclass', 'subtype', 'RL#', 'age', 'chem', 'concat_id', 'Race', 'Brain Regions*', 'Cause of Death', 'ICD-10 Code', 'ICD-10 category', 'Oxygen/No Oxygen', 'Date-of-Collection', 'Collection_year', 'Library Prep Date', 'Library Prep Lot', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'doublet_score', 'log10_gene_counts', 'log10_UMI_counts', 'percent_mito', 'percent_ribo', 'n_counts', 'leiden', 'mat/dev', 'cell_type', 'major_clust', 'sub_clust', 'combined-leiden']]
	adata.obs['Age'] = adata.obs['Age'].astype(str)
	return(adata)




###### standardize library sizes by removing nucs with low libs and downsampling the rest ######
def standardize_libs( adata, tar_lib_sz=1_000, seed=123):
	# make sure counts for each nuc are current
	sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True)
	# print number of nucs to be removed
	print( "Number of low lib nucs removed per run:")
	low_mask = ( adata.obs['total_counts']<=tar_lib_sz)
	print( pd.value_counts( adata.obs.SubID.values[low_mask]).to_string())
	
	# remove nucs with low library sizes
	#### comment out line below to keep low UMI count nuclei ####
	sc.pp.filter_cells( adata, min_counts=tar_lib_sz, copy=False)
	
	# downsample nuc total counts over target library size
	sc.pp.downsample_counts( adata, counts_per_cell=tar_lib_sz, copy=False, replace=False, random_state=seed)
	# due to scanpy bug in this version have to reset dtype
	adata.X = adata.X.astype( float)
	return



def whole_data_analysis(input_files, batch_labels, out_label, down_sample, hvgFile, n_pcs, n_neighbors, batch_key):
	# step0: Prepare data
	# step0_0: Whether need to concatenate data (Has completed)
	if (fnmatch.fnmatch(input_files, '*_adata.h5ad')):
		adata = sc.read(input_files)
	else:
		if (len(input_files.split('@')) > 1):
			adata = merge_data(input_files, batch_labels)
			adata = manage_obs(adata)
		else:
			adata = sc.read(input_files)
			adata.obs['Batch'] = batch_labels
		adata.write(out_label + '_adata.h5ad')

	# step2. Whether need to downsample
	if (down_sample != 'no'):
		# step2_0. Save raw current count matrix to raw
		adata.layers['raw-cts_pre-ds'] = adata.X.copy()
		adata = adata[:]

		# step2_1. only want to sample from genes that would pass preprocessing without downsampling
		sc.pp.filter_genes( adata, min_cells=5)
		print('Data shape before downsampling: ' + str(adata.shape))

		# step2_2. standardize the library
		standardize_libs( adata, tar_lib_sz = int(down_sample))
		print('Data shape after downsampling: ' + str(adata.shape))

		# step2_3. filter genes post downsampling
		sc.pp.filter_genes( adata, min_cells=5)
		print('Data shape after downsampling and filtering genes: ' + str(adata.shape))
		adata.write(out_label + '_adata.h5ad')

	# step3. Normalize (with standard library size, just scaling data) and scale to CPM
	# step3_0. re-calculate umi and gene counts per nuclei
	sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True)

	# step3_1. scale to CPM
	sc.pp.normalize_total( adata, target_sum=1e6, inplace=True)

	# step4. Variance Stabilizing Transformation - VST
	sc.pp.log1p( adata, copy=False)

	# step5. Feature selection
	if (hvgFile == 'scanpy'):
		n_hvg = 5_000
		adata.uns['log1p']["base"] = None
		sc.pp.highly_variable_genes( adata, n_top_genes=n_hvg, n_bins=20, flavor='seurat', inplace=True)
	else:
		hvgs = list(set(pd.read_csv(hvgFile, header=None)[0]) & set(adata.var_names))
		adata.var['highly_variable'] = False
		adata.var.loc[hvgs, 'highly_variable'] = True
	ax = sc.pl.highly_variable_genes( adata, show=False)
	plt.savefig( f"{out_label}_highly_variable_genes.png", format='png', bbox_inches='tight')

	high_adata = adata[:,adata.var.highly_variable.values]
	print('Shape of high_adata: ' + str(high_adata.shape))

	# step6. PCA
	if(n_pcs == 'half'):
		n_comps = 500
		sc.pp.pca( high_adata, n_comps=n_comps) #, use_highly_variable=True)
		n_pcs = ( np.cumsum( high_adata.uns['pca']['variance_ratio'])<0.50).sum()
		if( n_pcs==n_comps):
			n_comps = 1000
			sc.pp.pca( high_adata, svd_solver='arpack', n_comps=n_comps) #, use_highly_variable=True)
			n_pcs = ( np.cumsum( high_adata.uns['pca']['variance_ratio'])<0.50).sum()
			if( n_pcs==n_comps):
				n_pcs = None
				print( "re-run PCA with more components")
	else:
		n_pcs = int(n_pcs)
		sc.pp.pca( high_adata, svd_solver='arpack', n_comps=n_pcs)

	print('n_pcs used is: ' + str(n_pcs))
	ax = sc.pl.pca_variance_ratio( high_adata, show=False)
	plt.savefig(f"{out_label}_pca_variance_ratio.png", format='png', bbox_inches='tight')
	ax = sc.pl.pca( high_adata, color=['stage_id'], components=['1,2','3,2'], legend_fontsize=10, wspace=0.25, show=False)
	plt.savefig(f"{out_label}_pca.png", format='png', bbox_inches='tight')

	# step7. UMAP
	# step7_0. find nearest neighbors
	sc.pp.neighbors( high_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=6, use_rep="X_pca")
	# step7_a. calc umap
	sc.tl.umap( high_adata, random_state=6)
	ax = sc.pl.umap( high_adata, color=['cell_type_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_cell_type_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subclass_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subclass_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subtype_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subtype_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['dev_mat_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_dev_mat_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['Batch'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_batch.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['SubID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_SubID.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['poolID_new'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_poolID_new.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['prep'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_prep.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_stage_id.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['class'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_class.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subclass'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subclass.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subtype'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subtype.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['major_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_major_clust.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['sub_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_sub_clust.png", format='png', bbox_inches='tight')
	high_adata.write(out_label + '_high_adata.h5ad')
	print('UMAP calculation is finished!')

	# step8. removing batch effects
	if(batch_key != "no"):
		high_adata.obsm["X_pca_harmony"] = harmonize(high_adata.obsm["X_pca"], high_adata.obs, batch_key=batch_key, max_iter_harmony=20)
		high_adata.write(out_label + '_high_adata_rbe.h5ad')
		print('Batch effects removing is finished!')
		# step10_0. UMAP
		# step10_0_0. find nearest neighbors
		sc.pp.neighbors( high_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=6, use_rep="X_pca_harmony")
		# step10_0_1. calc umap
		sc.tl.umap( high_adata, random_state=6)
		ax = sc.pl.umap( high_adata, color=['cell_type_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_cell_type_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subclass_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subclass_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subtype_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subtype_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['dev_mat_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_dev_mat_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['Batch'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_batch.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['SubID'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_SubID.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['poolID_new'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_poolID_new.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['prep'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_prep.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_stage_id.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['class'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_class.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subclass'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subclass.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subtype'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subtype.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['major_clust'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_major_clust.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['sub_clust'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_sub_clust.png", format='png', bbox_inches='tight')
		high_adata.write(out_label + '_high_adata_rbe.h5ad')
		print('UMAP calculation after removing batch effects is finished!')

	# step9. clustering
	sc.tl.leiden( high_adata, resolution=3.0, random_state=0)
	sc.tl.rank_genes_groups( high_adata, 'leiden')
	# plot leiden clustering in black
	plt_adata = sc.pp.subsample( high_adata, fraction=1.0, random_state=0, copy=True)
	feat = 'leiden'
	ax = sc.pl.umap( plt_adata, color=feat, palette = 'tab20', size=2.5, sort_order=False, legend_loc='on data', show=False)
	plt.savefig(f"{out_label}_umap_leiden_black.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( plt_adata, color=feat, palette = 'tab20', size=2.5, sort_order=False, legend_loc='on data', show=False)
	# Adjust the text color of each text object in the Axes
	for text in ax.texts:
		text.set_color('white')
	
	plt.savefig(f"{out_label}_umap_leiden_white.png", format='png', bbox_inches='tight')

	# step10. update adata
	# set umap and pca to adata
	adata.obsm = high_adata.obsm
	# set connectivities and distances from umap to adata
	adata.obsp = high_adata.obsp
	# set clustering to adata
	adata.obs['leiden'] = high_adata.obs['leiden']
	# set downsampled count data to layer
	adata.layers['ds_norm_cts'] = adata.X.copy()
	adata.write(out_label + '_adata.h5ad')


def main():
	parser = OptionParser()
	parser.add_option('-i', '--input_files', dest='input_files', help='Input files.')
	parser.add_option('-b', '--batch_labels', dest='batch_labels', help='Batch labels.')
	parser.add_option('-o', '--out_label', dest='out_label', help='Output file.')
	parser.add_option('-d', '--down_sample', dest='down_sample', help='Down sample or not?')
	parser.add_option('-v', '--hvgFile', dest='hvgFile', help='Highly variable genes.')
	parser.add_option('-p', '--n_pcs', dest='n_pcs', help='Number of PCs.')
	parser.add_option('-n', '--n_neighbors', dest='n_neighbors', help='Number of PCs.')
	parser.add_option('-k', '--batch_key', dest='batch_key', help='Batch or not?')

	(options,args) = parser.parse_args()
	print(options)
	whole_data_analysis(options.input_files, options.batch_labels, options.out_label, options.down_sample, options.hvgFile, options.n_pcs, int(options.n_neighbors), options.batch_key)

main()
