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
sc.settings.n_jobs = 30
sc.settings.set_figure_params( dpi=200, fontsize=6)


###### Functions used to generate distinct list of hex colors ######
def get_random_color( pastel_factor=0.5):
	return [( x+pastel_factor)/( 1.0+pastel_factor) for x in [random.uniform( 0,1.0) for i in [1,2,3]]]

def color_distance( c1,c2):
	return sum( [abs( x[0]-x[1]) for x in zip( c1,c2)])

def generate_new_color( existing_colors, pastel_factor=0.5):
	max_distance = None
	best_color = None
	for i in range( 0,7_000):
		color = get_random_color( pastel_factor=pastel_factor)
		if not existing_colors:
			return color
		best_distance = min( [color_distance( color,c) for c in existing_colors])
		if not max_distance or best_distance > max_distance:
			max_distance = best_distance
			best_color = color
	return best_color

def random_hex_colors( n_colors, random_seed=0):	#Example:
	#To make your color choice reproducible, uncomment the following line:
	random.seed( random_seed) # is what we used for GABA fig 12345678
	colors = []
	for i in range( 0, n_colors):
		colors.append( generate_new_color( colors, pastel_factor=random.uniform( 0.0,0.10)))
	new_colors = np.array( [mpl.colors.to_hex(  ii) for ii in colors ])
	return( new_colors)


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
		digit = float(re.findall(r'\d+\.\d+|\d+', age_itr)[0])
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
	stage_id  = np.zeros_like( adata.obs['SubID'].values.astype(str), dtype = '<U20')
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
	stage_id[(numer_age>=12) & (numer_age<20)] = 'Adolescence'
	# Adulthood - 20yr up to 40
	stage_id[(numer_age>=20) & (numer_age<40)] = 'Young_Adulthood'
	# Late Adulthood - 40yr up to 60
	stage_id[(numer_age>=40) & (numer_age<60)] = 'Middle_Adulthood'
	# Late_Adulthood - 60yr and up
	stage_id[numer_age>=60] = 'Late_Adulthood'

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
	stage_order_full = ['Fetal', 'Neonatal', 'Infancy', 'Childhood', 'Adolescence', 'Young_Adulthood', 'Middle_Adulthood', 'Late_Adulthood']
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
	adata.obs = adata.obs[['Batch', 'SubID', 'Age', 'Sex', 'PMI', 'cell_type_uni', 'subclass_uni', 'subtype_uni', 'dev_mat_uni', 'numerical_age', 'stage_id', 'Channel', 'Source', 'poolID', 'prep', 'rep', 'class', 'subclass', 'subtype', 'RL#', 'age', 'chem', 'concat_id', 'Race', 'Brain Regions*', 'Cause of Death', 'ICD-10 Code', 'ICD-10 category', 'Oxygen/No Oxygen', 'Date-of-Collection', 'Collection_year', 'Library Prep Date', 'Library Prep Lot', 'n_genes_by_counts', 'log1p_n_genes_by_counts', 'total_counts', 'log1p_total_counts', 'doublet_score', 'log10_gene_counts', 'log10_UMI_counts', 'percent_mito', 'percent_ribo', 'n_counts', 'leiden', 'mat/dev', 'cell_type', 'major_clust', 'sub_clust']]
	adata.obs['Age'] = adata.obs['Age'].astype(str)
	return(adata)


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




###### Return adata with nuclei masked removed and print number removed from each SubID ######
def mask_cells( adata, mask):
	sample_cts = pd.value_counts( adata.obs.SubID[~mask])
	print( "Number of removed nuclei per sample:")
	print( sample_cts[sample_cts>0].to_string())
	print( "Total removed = ", sample_cts.values.sum())
	adata = adata[mask,:]
	return( adata)


###### plot log feature counts per nucs and return low nuc mask by median absolute deviation
# will return bool mask for nuc with less than 3 MADs of gene counts, but will not remove nuclei over 300 gene counts ######
def plot_mask_feature_per_nuc( adata, figure_label, feature, max_th, xlab, mad_scale=3.0):
	# use scanpy's built QCs metrics
	sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True);
	# max threshold in log space
	max_th = np.log10( max_th)
	low_feat_mask = np.array( [False]*adata.shape[0], dtype=bool)
	# seaborn plot settings
	sns.set()
	sns.set_context( "paper")
#	 palette = itertools.cycle( sns.color_palette( 'colorblind'))
	# subplot params
	n	= len( adata.uns['SubID_order'])
	cols = round( math.sqrt( n))
	rows = math.ceil( n / cols)
	size = cols * rows
	fig, axes = plt.subplots( cols, rows, figsize=( 5*cols, 4*rows), sharex=True)
	# loop through each sample
	for itr, (sub_itr, ax) in enumerate( zip( adata.uns['SubID_order'], axes.flatten())):
		sub_mk = adata.obs['SubID'].values==sub_itr
		# get feature counts as calculated by scanpy
		sub_feats = np.log10( adata.obs[feature][sub_mk] )
		# calc median absolute deviation
		med_th = min( max_th, np.median( sub_feats) - stats.median_abs_deviation( sub_feats)*mad_scale)
		low_feat_mask[sub_mk] = sub_feats<med_th
		# plot on new subplot for each histogram
		plt.sca(ax)
		sd = sns.distplot( sub_feats, color=adata.uns['SubID_colors_dict'][sub_itr])
		x, y = sd.get_lines()[0].get_data()
		plt.title( sub_itr)
		plt.vlines( med_th, 0, max(y), linestyles='dashed')
		plt.ylabel( "Nuclei Count Density", fontsize=8)
		plt.xlabel( xlab, fontsize=8)
	# clear rest of graphs
	for ii in range( itr+1, size):
		ax = axes.flatten()[ii]
		ax.axis('off')
	plt.savefig( f"{figure_label}", format='png', bbox_inches='tight')
	return( low_feat_mask)


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



def whole_data_analysis(input_files, rand_state, batch_labels, out_label, down_sample, hvgInfo, n_pcs, n_neighbors, batch_key):
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
		print("step0: data prepare is finished!")

		# step1. Nuclei and gene quality control
		# step1_0. QC nuclei
		# get scanpy quality metrics
		sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True)
		# remove genes with fewer than 5 total counts
		sc.pp.filter_genes( adata, min_counts=5, inplace=True)
		print('Data shape after removing genes with fewer than 5 total counts: ' + str(adata.shape))
	
		# step1_1. Remove dublets with scrublet
		# Remove doublets using scrublet on a per sample basis
		dub_mk = np.array( [False]*adata.shape[0], dtype=bool)
		dub_sc = np.array( [False]*adata.shape[0], dtype=float)
		for sub_itr in adata.uns['SubID_order']:
			sub_mk = adata.obs['SubID'].values==sub_itr
			adata_itr = adata[sub_mk,:]
			scrub = scr.Scrublet( adata_itr.X.tocsc(), expected_doublet_rate=0.1)
			doublet_scores, predicted_doublets = scrub.scrub_doublets( min_counts=2, 
				min_cells=3, 
				min_gene_variability_pctl=85, 
				n_prin_comps=10,
				verbose=False)
			dub_mk[sub_mk] = predicted_doublets
			dub_sc[sub_mk] = doublet_scores
			 
		print( "Number of potential dublet nucs found in each sample:")
		sample_cts = pd.value_counts( adata.obs.SubID[dub_mk])
		print( sample_cts.to_string())
		print( "Total doublets found = ", sample_cts.values.sum())

		# add scrublet doublet scores to adata
		adata.obs['doublet_score'] = dub_sc
		# remove predicted doublets
		adata = adata[~dub_mk]
		print('Data shape after removing doublets: ' + str(adata.shape))
	
		sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True)

		# step1_2. Library sizes and feature counts per sample
		# maximum number of colors in a colormap is tab20, since we have over 20 samples 
		# generate randomly distinct colors for each SubID-from functions in misc_code.py file
		adata.uns['SubID_colors'] = random_hex_colors( n_colors=len( adata.uns['SubID_order']), random_seed=rand_state)
		adata.uns['SubID_colors_dict'] = dict( zip( adata.obs['SubID'].cat.categories.tolist(), adata.uns['SubID_colors']))
	
		adata.obs['log10_gene_counts'] = np.log10( adata.obs['n_genes_by_counts'])
		adata.obs['log10_UMI_counts']  = np.log10( adata.obs['total_counts'])
	
		with rc_context({'figure.figsize': (25,10)}):
			ax = sc.pl.violin( adata, ['log10_gene_counts'], groupby='SubID', stripplot=True, 
				inner='box', order=adata.uns['SubID_order'], rotation=70, xlabel='SubID', show=False)
			plt.savefig( f"{out_label}_violin_feature-cts.png", format='png', bbox_inches='tight')
	
		with rc_context({'figure.figsize': (25,10)}):
			ax = sc.pl.violin( adata, ['log10_UMI_counts'], groupby='SubID', stripplot=True, 
				inner='box', order=adata.uns['SubID_order'], rotation=70, xlabel='SubID', show=False)
			plt.savefig( f"{out_label}_violin_library-sizes.png", format='png', bbox_inches='tight')
	
		# step1_3. Remove nuclei with low feature counts
		low_feat_mk = plot_mask_feature_per_nuc( adata, figure_label = out_label + '_nuclei_with_low_feature_counts.png', feature='n_genes_by_counts', max_th=300, xlab="log10 genes by counts")
		# remove cells with low features per cell
		adata = mask_cells( adata, ~low_feat_mk)
	
		# step1_4. Remove overly expressed genes before normalization
		# MALAT1 is known to be overexpressed and can affect normalization. Highly expressed genes also have a large variance, so normal flucuations can have an outsized effect on library normalization. 
		adata = adata[:,adata.var_names!='MALAT1']
		adata.write(out_label + '_adata.h5ad')
		print("step1: dublet removing is finished!")

	# step2. Whether need to downsample
	if (down_sample != 'no'):
		# step2_0. Save raw current count matrix to raw
		adata.layers['raw-cts_pre-ds'] = adata.X.copy()
		adata = adata[:]

		# step2_1. only want to sample from genes that would pass preprocessing without downsampling
		sc.pp.filter_genes( adata, min_cells=5)
		print('Data shape before downsampling: ' + str(adata.shape))

		# step2_2. standardize the library
		standardize_libs( adata, tar_lib_sz = int(down_sample), seed = rand_state)
		print('Data shape after downsampling: ' + str(adata.shape))

		# step2_3. filter genes post downsampling
		sc.pp.filter_genes( adata, min_cells=5)
		print('Data shape after downsampling and filtering genes: ' + str(adata.shape))
		adata.write(out_label + '_adata.h5ad')
		print("step2: downsampling is finished!")

	adata = sc.read(out_label + '_adata.h5ad')
	# step3. Normalize (with standard library size, just scaling data) and scale to CPM
	# step3_0. re-calculate umi and gene counts per nuclei
	sc.pp.calculate_qc_metrics( adata, percent_top=None, inplace=True)

	# step3_1. scale to CPM
	sc.pp.normalize_total( adata, target_sum=1e6, inplace=True)

	# step4. Variance Stabilizing Transformation - VST
	sc.pp.log1p( adata, copy=False)

	# step5. Feature selection
	if ((len(hvgInfo.split('@')) > 1) &  (hvgInfo.split('@')[0] == 'scanpy')):
		n_hvg = int(hvgInfo.split('@')[1])
		adata.uns['log1p']["base"] = None
		adata_filtered = adata[:, (adata.var['gene_type-Aging'] == 'protein_coding') & (~adata.var['gene_chrom-Aging'].isin(['X', 'Y', 'MT']))].copy()
		sc.pp.highly_variable_genes( adata_filtered, n_top_genes=n_hvg, n_bins=20, flavor='seurat', inplace=True)
		adata.var['highly_variable'] = adata.var_names.isin(adata_filtered.var_names[adata_filtered.var['highly_variable']])
	else:
		hvgs = list(set(pd.read_csv(hvgInfo, header=None)[0]) & set(adata.var_names))
		adata.var['highly_variable'] = False
		adata.var.loc[hvgs, 'highly_variable'] = True
	# ax = sc.pl.highly_variable_genes( adata, show=False)
	# plt.savefig( f"{out_label}_highly_variable_genes.png", format='png', bbox_inches='tight')

	adata.write(out_label + '_adata.h5ad')
	high_adata = adata[:,adata.var.highly_variable.values]
	print('Shape of high_adata: ' + str(high_adata.shape))
	print("step5: feature selection is finished!")

	# step6. PCA
	if(n_pcs == 'half'):
		n_comps = 600
		sc.pp.pca( high_adata, svd_solver='arpack', random_state=rand_state, n_comps=n_comps) #, use_highly_variable=True)
		n_pcs = ( np.cumsum( high_adata.uns['pca']['variance_ratio'])<0.50).sum()
		if( n_pcs==n_comps):
			n_comps = 1000
			sc.pp.pca( high_adata, svd_solver='arpack', random_state=rand_state, n_comps=n_comps) #, use_highly_variable=True)
			n_pcs = ( np.cumsum( high_adata.uns['pca']['variance_ratio'])<0.50).sum()
			if( n_pcs==n_comps):
				n_pcs = None
				print( "re-run PCA with more components")
	else:
		n_pcs = int(n_pcs)
		sc.pp.pca( high_adata, svd_solver='arpack', random_state=rand_state, n_comps=n_pcs)

	print('n_pcs used is: ' + str(n_pcs))
	ax = sc.pl.pca_variance_ratio( high_adata, show=False)
	plt.savefig(f"{out_label}_pca_variance_ratio.png", format='png', bbox_inches='tight')
	ax = sc.pl.pca( high_adata, color=['stage_id'], components=['1,2','3,2'], legend_fontsize=10, wspace=0.25, show=False)
	plt.savefig(f"{out_label}_pca.png", format='png', bbox_inches='tight')
	high_adata.write(out_label + '_high_adata.h5ad')
	print("step6: PCA is finished!")

	# step7. UMAP
	# step7_0. find nearest neighbors
	sc.pp.neighbors( high_adata, n_neighbors = n_neighbors, n_pcs = n_pcs, random_state = rand_state, use_rep="X_pca")
	# step7_a. calc umap
	sc.tl.umap( high_adata, random_state = rand_state)
	high_adata.write(out_label + '_high_adata.h5ad')
	print('step7: UMAP calculation is finished!')
	ax = sc.pl.umap( high_adata, color=['cell_type_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_cell_type_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subclass_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subclass_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subtype_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subtype_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['dev_mat_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_dev_mat_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['Batch'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_batch.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['SubID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_SubID.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['poolID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_poolID.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['prep'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_prep.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_stage_id.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['class'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_class.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subclass'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subclass.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['subtype'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_subtype.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['major_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_major_clust.png", format='png', bbox_inches='tight')
	ax = sc.pl.umap( high_adata, color=['sub_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_sub_clust.png", format='png', bbox_inches='tight')

	adata = sc.read(out_label + '_adata.h5ad')
	high_adata = sc.read(out_label + '_high_adata.h5ad')

	# step8. UMAT
	high_adata = umat( high_adata.uns['stage_order'], high_adata, pcaType='X_pca', n_pcs=int(n_pcs), n_neighbors=n_neighbors, rand=rand_state) #123, #12 #12345 #123456
	high_adata.write(out_label + '_high_adata.h5ad')
	print('step8: UMAT calculation is finished!')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['cell_type_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_cell_type_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subclass_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subtype_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['dev_mat_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_dev_mat_uni.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['Batch'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_batch.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['SubID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_SubID.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['poolID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_poolID.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['prep'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_prep.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_stage_id.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['class'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_class.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subclass.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_subtype.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['major_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_major_clust.png", format='png', bbox_inches='tight')
	ax = sc.pl.embedding( high_adata, basis='umat', color=['sub_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_sub_clust.png", format='png', bbox_inches='tight')


	# step9. removing batch effects
	if(batch_key != "no"):
		high_adata.obsm["X_pca_harmony"] = harmonize(high_adata.obsm["X_pca"], high_adata.obs, batch_key=batch_key, max_iter_harmony=20)
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step9: Batch effects removing is finished!')
		# step9_0. UMAP
		# step9_0_0. find nearest neighbors
		sc.pp.neighbors( high_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state = rand_state, use_rep="X_pca_harmony")
		# step9_0_1. calc umap
		sc.tl.umap( high_adata, random_state = rand_state)
		ax = sc.pl.umap( high_adata, color=['cell_type_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_cell_type_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subclass_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subclass_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subtype_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subtype_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['dev_mat_uni'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_dev_mat_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['Batch'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_batch.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['SubID'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_SubID.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['poolID'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_poolID.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['prep'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_prep.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_stage_id.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['class'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_class.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subclass'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subclass.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['subtype'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_subtype.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['major_clust'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_major_clust.png", format='png', bbox_inches='tight')
		ax = sc.pl.umap( high_adata, color=['sub_clust'], palette='tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umap_rbe_sub_clust.png", format='png', bbox_inches='tight')
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step9_0_UMAP calculation after removing batch effects is finished!')
		# step9_1. UMAT calculation
		high_adata = umat( high_adata.uns['stage_order'], high_adata, pcaType='X_pca_harmony', n_pcs=int(n_pcs), n_neighbors=n_neighbors, rand=rand_state) #123, #12 #12345 #123456
		high_adata.write(out_label + '_high_adata.h5ad')
		print('step9_1_UMAT calculation after removing batch effects is finished!')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['cell_type_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_cell_type_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_subclass_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_subtype_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['dev_mat_uni'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_dev_mat_uni.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['Batch'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_batch.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['SubID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_SubID.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['poolID'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_poolID.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['prep'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_prep.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['stage_id'], legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_stage_id.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['class'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_class.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['subclass'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_subclass.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['subtype'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_subtype.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['major_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_major_clust.png", format='png', bbox_inches='tight')
		ax = sc.pl.embedding( high_adata, basis='umat', color=['sub_clust'], palette = 'tab20', legend_fontsize=4, add_outline=False, size=5, legend_fontoutline=0.2, show=False);plt.savefig(f"{out_label}_umat_rbe_sub_clust.png", format='png', bbox_inches='tight')

	# step10. update adata
	# set umap and pca to adata
	adata.obsm = high_adata.obsm
	# set connectivities and distances from umap to adata
	adata.obsp = high_adata.obsp
	# set downsampled count data to layer
	if (down_sample != "no"):
		adata.layers['ds_norm_cts'] = adata.X.copy()
	adata.write(out_label + '_adata.h5ad')


def main():
	parser = OptionParser()
	parser.add_option('-i', '--input_files', dest='input_files', help='Input files.')
	parser.add_option('-r', '--rand_state', dest='rand_state', help='Random state.')
	parser.add_option('-b', '--batch_labels', dest='batch_labels', help='Batch labels.')
	parser.add_option('-o', '--out_label', dest='out_label', help='Output file.')
	parser.add_option('-d', '--down_sample', dest='down_sample', help='Down sample or not?')
	parser.add_option('-v', '--hvgInfo', dest='hvgInfo', help='Highly variable genes.')
	parser.add_option('-p', '--n_pcs', dest='n_pcs', help='Number of PCs.')
	parser.add_option('-n', '--n_neighbors', dest='n_neighbors', help='Number of PCs.')
	parser.add_option('-k', '--batch_key', dest='batch_key', help='Batch or not?')

	(options,args) = parser.parse_args()
	print(options)
	whole_data_analysis(options.input_files, int(options.rand_state), options.batch_labels, options.out_label, options.down_sample, options.hvgInfo, options.n_pcs, int(options.n_neighbors), options.batch_key)

main()