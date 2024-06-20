colors_batch = {
    'Aging': '#2277B2',
    'Lister': '#FF7F0E'
}

colors_class = {
    'EN': '#197EC0',
    'IN': '#1A9993',
    'Immune': '#C80813',
    'Mural': '#FD8CC1',
    'Endo': '#FD7446',
    'Astro': '#D2AF81',
    'Oligo': '#D5E4A2',
    'OPC': '#FED439',
    'Others': '#BFBFBF'
}

colors_subclass = {
    'EN_L2_3_IT': '#659EC7',
    'EN_L3_5_IT_1': '#95B9C7',
    'EN_L3_5_IT_2': '#6495ED',
    'EN_L3_5_IT_3': '#79BAEC',
    'EN_L5_6_NP': '#0020C2',
    'EN_L5_ET': '#4863A0',
    'EN_L6_CT': '#3BB9FF',
    'EN_L6_IT_1': '#B4CFEC',
    'EN_L6_IT_2': '#1589FF',
    'EN_L6B': '#488AC7',
    'PN_dev': '#0433FF',
    'IN_ADARB2': '#4E8975',
    'IN_LAMP5_LHX6': '#7BCCB5',
    'IN_LAMP5_RELN': '#3B9C9C',
    'IN_PVALB': '#B2C248',
    'IN_PVALB_CHC': '#008080',
    'IN_SST': '#728C00',
    'IN_VIP': '#89C35C',
    'Micro': '#F75D59',
    'PVM': '#C11B17',
    'Adaptive': '#6F4E37',
    'PC': '#E0B0FF',
    'SMC': '#4E387E',
    'VLMC': '#B93B8F',
    'Endo': '#FFA62F',
    'Astro': '#C19A6B',
    'Oligo': '#ECE5B6',
    'OPC': '#FFF380'
}

# midcolor
colors_stage_id = {
    'Fetal': '#F4ECD3', 
    'Neonatal': '#E9D8A6', 
    'Infancy': '#F1D4A4', 
    'Childhood': '#F9CFA1', 
    'Adolescence': '#FFCC66', 
    'Young_Adulthood': '#EE9B00', 
    'Middle_Adulthood': '#CA6702', 
    'Late_Adulthood': '#8C510A'
}
stage_id_order = ['Fetal', 'Neonatal', 'Infancy', 'Childhood', 'Adolescence', 'Young_Adulthood', 'Middle_Adulthood', 'Late_Adulthood']


colors_traj = {
    'EN_L2_3_IT': '#659EC7',
    'EN_L3_5_IT_1': '#95B9C7',
    'EN_L3_5_IT_2': '#6495ED',
    'EN_L3_5_IT_3': '#79BAEC',
    'EN_L6_CT': '#3BB9FF',
    'EN_L6_IT_1': '#B4CFEC',
    'EN_L6_IT_2': '#1589FF',
    'EN_L6B': '#488AC7',
    'EN_L5_6_NP': '#0020C2',
    'IN_ADARB2': '#4E8975',
    'IN_ADARB2_RAB37': '#58937F', 
    'IN_ADARB2_COL12A1': '#447F6B', 
    'IN_ADARB2_SYT10': '#4E9375', 
    'IN_ADARB2_SV2C': '#4E7F75',
    'IN_LAMP5_LHX6': '#7BCCB5',
    'IN_LAMP5_RELN': '#3B9C9C',
    'IN_PVALB': '#B2C248',
    'IN_PVALB_CHC': '#008080',
    'IN_SST': '#728C00',
    'IN_VIP': '#89C35C',
    'Micro': '#F75D59',
    'Astro_PLSCR1': '#9A7B56', 
    'Astro_ADAMTSL3' : '#F4CD9E',
    'Oligo': '#ECE5B6'
}


trait_info = {
    'sz3': 'SCZ_2022_Trubetskoy',
    'bip2': 'BD_2021_Mullins',
    'asd': 'ASD_2019_Grove',
    'adhd_ipsych': 'ADHD_2023_Demontis',
    'mdd_ipsych': 'MDD_2023_AlsBroad',
    'ocd': 'OCD_2018_IOCDF_GC',
    'insomn2': 'Insomnia_2019_Jansen',
    'alcohilism_2019': 'Alcoholism_2019_SanchezRoige',
    'tourette': 'Tourettes_2019_Yu',
    'intel': 'IQ_2018_Savage',
    'eduAttainment': 'Education_2018_Lee',

    
    'alzBellenguezNoApoe': 'AD_2022_Bellenguez',
    'ms': 'MS_2019_IMSGC',
    'pd_without_23andMe': 'PD_2019_Nalls',
    'migraines_2021': 'Migraines_2021_Donertas',
    'als2021': 'ALS_2021_vanRheenen',
    'stroke': 'Stroke_2018_Malik',
    'epilepsyFocal': 'Epilepsy_2018_ILAECCE',

                               
    'dm2': 'T2D_2012_Morris',
    'ibd': 'IBD_2015_Liu',
    'ra': 'RArthritis_2013_Okada',
    'lipidCholTotal': 'Cholesterol_2013_Wiler',
    'obesity': 'Obesity_2019_Watanabe',
    'uc': 'UC_2011_Anderson'
}

colors_plot = ["#2C7BB6", "white", "#D7191C"]
cmap_plot = mcolors.LinearSegmentedColormap.from_list("custom_cmap", colors_plot)

colors_niche = {
    '0': '#F4CD9E', 
    '1': '#9A7B56'}

from matplotlib.colors import ListedColormap
cmap_niche = ListedColormap([colors_niche[val] for val in colors_niche.keys()])
