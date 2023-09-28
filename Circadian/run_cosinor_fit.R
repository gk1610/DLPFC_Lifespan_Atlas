#ENTER DATA AND SET UP FOR RUNNING ANALYSIS

#logcpm values format: column names are subject ids and row names are gene ids
data <- read.csv("~/Desktop/Panos Collab/v0_pseudoExpr/v0_Astro.csv", row.names = 1)
#TOD values and matching subject IDs
pheno <- read.csv ("~/Desktop/Panos Collab/pb_metadata.csv", row.names = 1)

#For the single cell pilot I did a bunch of filtering and subsetting - probably will not be necessary here but just in case:

#add in columns with "1" for each cell type - filter by cell type 
pheno_cell <- filter(pheno, pheno$Astro == "1")

#filter by group - NA = Non-AD, AD = Alzheimer's, Com = Small comparison group (no SZ)
pheno_cell_NA <- filter(pheno_cell, pheno_cell$Dx == "Control")
pheno_cell_AD <- filter(pheno_cell, pheno_cell$Dx == "AD")
pheno_cell_Com <- filter(pheno_cell, pheno_cell$COM == "1")

#THIS IS PROBABLY THE PART WE NEED 
#Make a data table- sub.data = data[,match(pheno$ID,colnames(data))]
cell.data = data[,match(pheno_cell$ID,colnames(data))]
colnames(data)

#Filter by group 
cell_Com.data = cell.data[,which((colnames(cell.data) %in% pheno_cell_Com$ID) == TRUE)]

#Way to check your names - all(pheno_Astro$ID == colnames(data))
all(pheno_cell_Com$ID == colnames(cell_Com.data))
SubID <- pheno_Astro_Com$SubID

#END DATA SET UP 


#ACTUAL RUNNING OF THE COSINOR CODE
out.list = mclapply(1:nrow(cell_Com.data), function(i){ 
  if(i%%1000==0) print(i)
  res = one_cosinor_OLS(tod = as.numeric(pheno_cell_Com$TOD), y = as.numeric(cell_Com.data[i,]))
  res.onerow = as.data.frame(list(offset = res$M$est, offset.ll = res$M$CI[1], offset.ul = res$M$CI[2], 
                                  A = res$A$est, A.sd = res$A$sd, 
                                  phase = res$phase$est, phase.sd = res$phase$sd, 
                                  pvalue = res$test$pval, R2 = res$test$R2))
  return(res.onerow)
}, mc.cores = 1) 

out.tab = do.call(rbind.data.frame, out.list)
out.tab$qvalue = p.adjust(out.tab$pvalue, "BH")
out.tab$peak = 24 - out.tab$phase*24/(2*pi)
row.names(out.tab) = row.names(cell.data)
setwd("~/Desktop/Panos Collab/v0_rhythmicity/Astro//")
write.csv(out.tab,paste0("observed_Astro_Com.csv"))
