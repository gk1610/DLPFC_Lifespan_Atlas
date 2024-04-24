library(ggplot2)
library(RColorBrewer)

# step0. Load MAGMA-associated scripts
.libPaths(c("/sc/arion/projects/roussp01a/jaro/programs/R_libs_4_2", .libPaths()))
source("/sc/arion/projects/roussp01a/jaro/atacseq_ad/NYGC_AD_R01/downstream_magma.R")

# step0_1. Traits used
Neurological <- c("alzBellenguezNoApoe", "ms", "pd_without_23andMe", "epilepsyFocal", "migraines_2021", "stroke", "als2021")
Psychiatric <- c("sz3", "bip2", "mdd_ipsych", "asd", "adhd_ipsych", "insomn2", "eduAttainment", "intel", "alcoholism_2019", "ocd", "tourette")
Others <- c("obesity", "dm2", "lipidCholTotal", "ra", "ibd", "uc")
scDRS_traits <- c(Neurological,Psychiatric,Others)

# step0_2. Plot LDsc / MAGMA heatmap (Modified based on script from Tereza)
ldscHeatmapPlotter <- function(df, plotCol="P", FDR_PER_TRAIT=F, myPadjustMethod="BH", markNominallySignificant=T,plotTextSize=11) {
    myPalette = colorRampPalette(c("#F7FCF5", "#E5F5E0", 
                                   "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45", 
                                   "#006D2C", "#00441B"), space = "Lab")
    df = data.frame(df)
    df$sumstatName <- recode(df$sumstatName, 
                             "alzBellenguezNoApoe"="AD_2022_Bellenguez",
                             "ms"="MS_2019_IMSGC",
                             "pd_without_23andMe"="PD_2019_Nalls",
                             "migraines_2021"="Migraines_2021_Donertas",
                             "als2021"="ALS_2021_vanRheenen",
                             "stroke"="Stroke_2018_Malik",
                             "epilepsyFocal"="Epilepsy_2018_ILAECCE",
                             
                             "sz3"="SCZ_2022_Trubetskoy",
                             "bip2"="BD_2021_Mullins",
                             "asd"="ASD_2019_Grove",
                             "adhd_ipsych"="ADHD_2023_Demontis",
                             "mdd_ipsych"="MDD_2023_AlsBroad",
                             "ocd"="OCD_2018_IOCDF_GC",
                             "insomn2"="Insomnia_2019_Jansen",
                             "alcohilism_2019"="Alcoholism_2019_SanchezRoige",
                             "tourette"="Tourettes_2019_Yu",
                             "intel"="IQ_2018_Savage",
                             "eduAttainment"="Education_2018_Lee",
                             # "anorexia"="Anorexia_2019_Watson",
                                                        
                             "dm2"="T2D_2012_Morris",
                             "ibd"="IBD_2015_Liu",
                             "ra"="RArthritis_2013_Okada",
                             "lipidCholTotal"="Cholesterol_2013_Wiler",
                             "obesity"="Obesity_2019_Watanabe",
                             "uc"="UC_2011_Anderson")
    
    df$plotLabel = ""
    if(markNominallySignificant) {   df$plotLabel[(10^-df[,plotCol]) < 0.05] = "·" }
    df$plotLabel[p.adjust((10^-df[,plotCol]),method=myPadjustMethod) < 0.05]="#"
    
    if(FDR_PER_TRAIT) {
      for(i in unique(df$annoID)) {
        df$adj.P.value[df$annoID==i] = fdrtool(df[,..plotCol][df$annoID==i], statistic=c("pvalue"), cutoff.method=c("fndr"), plot = F)$qval
      }
    } #else {
    #df$adj.P.value = p.adjust(df$P,method="BH")
    #}
    
    df$plotLabel[df$adj.P.value < 0.05] = "#"
    zz = ggplot(df, aes(sumstatName, annoName, fill = tempScoreCol)) + geom_tile() + scale_y_discrete(expand = c(0, 0)) + scale_x_discrete(expand = c(0, 0)) + xlab("Trait") + 
      ylab("Cluster") + theme_classic(base_size = plotTextSize) + theme(axis.text.x = element_text(colour = "black", angle = -45, hjust = 0)) + 
      coord_fixed() + theme(legend.title = element_text(size = 10, face = "bold"))
    
    print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradientn(colours = myPalette(100), name = "-logP") + geom_text(aes(label = plotLabel), size = plotTextSize * 0.55))
    print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradient2(high="#D7191C", mid="white", low="#2C7BB6", name = "-logP") + geom_text(aes(label = plotLabel), size = plotTextSize * 0.55))
    print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradientn(colours = myPalette(100), name = "-logP", trans = "sqrt") + geom_text(aes(label = plotLabel), size = plotTextSize * 0.55))
    print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradientn(colours = myPalette(100), name = "-logP"))
    print(zz + geom_tile(aes(fill = tempScoreCol)) + scale_fill_gradientn(colours = myPalette(100), name = "-logP", trans = "sqrt"))
}

# step0_3. Plot MAGMA results for specific trait (Modified based on script from Swadha)
plot_magma_trait <- function(magmasets, traits, annoName_order){
  magmasets_subset=magmasets[magmasets$gwasAcronym %in% traits, ]
  magmasets_subset$myCoefficient=magmasets_subset$BETA
  magmasets_subset$error_left = magmasets_subset$BETA - magmasets_subset$BETA_STD
  magmasets_subset$error_right = magmasets_subset$BETA + magmasets_subset$BETA_STD
  magmasets_subset$myLabel=""
  magmasets_subset$myLabel[magmasets_subset$P < .05] = "*"
  magmasets_subset$adj.P.value=p.adjust(magmasets_subset$P, method="BH")
  magmasets_subset$myLabel[magmasets_subset$adj.P.value<0.05]="#"
  magmasets_subset$VARIABLE[magmasets_subset$P<.05]
  magmasets_subset$gwasAcronym[magmasets_subset$P<.05]
  pallete1=brewer.pal(8,"Reds")
  myPalette = colorRampPalette(pallete1)

  magmasets_subset$annoName <- factor(magmasets_subset$VARIABLE, levels = annoName_order)
  g1=ggplot(data=magmasets_subset,aes(x = annoName,y = myCoefficient,color=-log10(P),label=myLabel)) + 
    geom_point(size=5)+geom_errorbar(aes(ymin = error_left,ymax = error_right),size=1)+ geom_text(color="black")+ylab(paste0("magma coeficient of ", traits)) + xlab("")#+facet_wrap(~Region,ncol=2)
  g1=g1+geom_hline(yintercept=0,linetype="dashed")+scale_y_continuous(expand = c(0,0))+scale_color_gradientn(colours = myPalette(10))
  
  good_plot=function(g1,sz,rot){
    g1=g1+theme_classic()+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz),axis.text.x=element_text(angle=rot,hjust=1,vjust=1))
    g1=g1+theme(legend.position="right",legend.title=element_blank())+theme(strip.text=element_text(face="bold",color="black"),strip.background=element_rect(fill="#eeeeee",color="#eeeeee"))
    g1+theme(axis.title= element_text(color="black", size=sz),axis.text = element_text(color="black", size=sz))+theme(aspect.ratio=1)
  }
  
  
  final_plot = good_plot(g1,12,90)
  print(final_plot)  
}

# step1. Prepare genesets (Official gene symbol to ensembl)
ASTSEMBL_ASTFO_hg38 = "/sc/arion/projects/roussp01b/resources/databases/gene-info/ensembl/my-downloads/muchEnsemblInfo_hg38.tsv.gz"

#make a table of gene-level ensembl info
ensemblInfo=read.delim(ASTSEMBL_ASTFO_hg38, stringsAsFactors=F)
ensemblGeneInfo=unique(ensemblInfo[,c("Ensembl.Gene.ID", "Strand", "Gene.End..bp.", "Gene.Start..bp.", "Associated.Gene.Name", "Gene.type", "Status..gene.", "Source..gene.", "Transcript.count", "Description", "Chromosome.Name","HGNC.symbol")])

cluster_1 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_1.csv", sep=",")$V1), "Ensembl.Gene.ID"])
cluster_2 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_2.csv", sep=",")$V1), "Ensembl.Gene.ID"])
cluster_3 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_3.csv", sep=",")$V1), "Ensembl.Gene.ID"])
cluster_4 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_4.csv", sep=",")$V1), "Ensembl.Gene.ID"])
cluster_5 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_5.csv", sep=",")$V1), "Ensembl.Gene.ID"])
cluster_6 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_6.csv", sep=",")$V1), "Ensembl.Gene.ID"])
cluster_7 <- unique(ensemblGeneInfo[ensemblGeneInfo$HGNC.symbol %in% as.character(read.table("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/AST_ordered_lineage_mean_peak_ks_7_cluster_7.csv", sep=",")$V1), "Ensembl.Gene.ID"])

genes <- list("Cluster_1" = cluster_1, "Cluster_2" = cluster_2, "Cluster_3" = cluster_3, "Cluster_4" = cluster_4, "Cluster_5" = cluster_5, "Cluster_6" = cluster_6, "Cluster_7" = cluster_7)


# step2. Run MAGMA
myAnalysis <- calcMagmaGeneSetPvals(outDir="/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/MAGMA_AST_ks_7", geneMetaSets=genes, pvalInfo=grabGenePvalInfo(selectedGeneTypeName="ensemblProtCodGenes35kb10kbAutosomesNoBmhc"), checkIfOutDirExists=F)

# step3. Plot based on MAGMA results
magmasets <- readr::read_tsv("/sc/arion/projects/CommonMind/aging/hui/files/freeze3/DEGs_moran_test/modified_10000/qvalue0.05/MAGMA_AST_ks_7/meta-files/geneSetResults.tsv.gz")

# step3_1. Plot LDsc / MAGMA heatmap
ngenes = lapply(unique(magmasets$VARIABLE), function(x) round(mean(magmasets[magmasets$VARIABLE==x,"NGENES"])))
names(ngenes) = unique(magmasets$VARIABLE)

magmasets = magmasets[(magmasets$gwasAcronym %in% scDRS_traits),]
magmasets = magmasets[order(magmasets$gwasAcronym),]
magmasets$sumstatName = magmasets$gwasAcronym
magmasets$sumstatName = ordered(magmasets$sumstatName, levels=scDRS_traits)
magmasets$annoName = factor(magmasets$VARIABLE, levels = rev(sort(unique(magmasets$VARIABLE))))

magmasets$adj.P.value = p.adjust(magmasets$P, method="BH") # alternative: fdrtool(magmasets$P, statistic=c("pvalue"), cutoff.method=c("fndr"), plot = F)$qval
magmasets$tempScoreCol = magmasets$LogP
ldscHeatmapPlotter(magmasets, plotCol="LogP", markNominallySignificant=T)

# step3_2. Plot MAGMA results for specific trait
plot_magma_trait(magmasets, "sz3", c("Cluster_1", "Cluster_2", "Cluster_3", "Cluster_4", "Cluster_5", "Cluster_6", "Cluster_7"))
