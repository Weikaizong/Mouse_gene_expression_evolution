library("edgeR") 
library('DESeq2')
setwd('D:/R/New_House_mouse_Transcriptome/Sex_biased_gene_analysis')
# Winsorization function
winsor.fun <- function(Y, X, quan) {
  N <- estimateSizeFactors(DESeqDataSetFromMatrix(Y, DataFrame(X), ~X))$sizeFactor
  P <- t(t(Y) / N)
  cut <- apply(P, 1, quantile, quan)
  Cut <- matrix(rep(cut, ncol(Y)), nrow(Y))
  ind <- P > Cut
  P[ind] <- Cut[ind]
  Y <- round(t(t(P) * N))
  return(Y)
}

Coding_Gene_List <- read.table("GFF110_Protein_Coding_Gene.txt")
Coding_Gene_List <- Coding_Gene_List$V1

##################################### Get SB genes with edgeR ###############################################

#Different populations use their corresponding Tissue IDs
Tissue_ID <- c("DOM_Bra", "DOM_Hea", "DOM_Kid", "DOM_Liv", "DOM_Tes-Ova", "DOM_Mam", "DOM_Vas-Ovi", "DOM_Epi-Ute")
for (i in 1:length(Tissue_ID)) {
  Target_Tissue <- Tissue_ID[i]
  Tissue_Count <- read.table(paste0("D:/R/New_House_mouse_Transcriptome/Sex_biased_gene_analysis/Input_count_file/Protein coding/DOM_Coding_FeatureCounts/", 
                                    Target_Tissue, "_FeatureCounts_Fragment.count", sep = ""), sep = "\t", header = TRUE)
  # rownames(Tissue_Count) <- Tissue_Count$Geneid
  rownames(Tissue_Count) <- Tissue_Count[, 1] 
  Tissue_Count <- Tissue_Count[, 7:ncol(Tissue_Count)]
  Tissue_Count <- as.matrix(Tissue_Count)
  Tissue_Count <- Tissue_Count[Coding_Gene_List, ]
  
  ### FR
  Tissue_Annotation <- read.table(paste0("Sex_biased_Analysis/FR-Male_vs_Female/", Target_Tissue, "_Sample_Annotation.tsv", sep = ""), header = TRUE, row.names = 1)
  Tissue_Count_Subset <- Tissue_Count[, rownames(Tissue_Annotation)]
  Tissue_Group <- Tissue_Annotation$condition 
  
  # Apply Winsorization
  Tissue_Count_Subset <- winsor.fun(Y = Tissue_Count_Subset, X = Tissue_Group, quan = 0.95)
  
  Tissue_DGEList <- DGEList(counts = Tissue_Count_Subset, group = Tissue_Group, remove.zeros = TRUE)
  
  ##use print.function(filterByExpr.default) to see more details on filterByExpr function
  # Here we set min.count = 10,large.n = 10 (>the small group size of 8); 
  # This means that we only keep genes with at least 10 counts/fragment (~1 CPM) within the samples of at least the size of small group (usually 8)
  Tissue_DGEList_Keep <- filterByExpr(Tissue_DGEList, min.count = 10, large.n = 10)
  Tissue_DGEList_Final <- Tissue_DGEList[Tissue_DGEList_Keep, , keep.lib.sizes = FALSE]
  Tissue_DGEList_Final <- calcNormFactors(Tissue_DGEList_Final)  # Normalize libary sizes with TMM
  dim(Tissue_DGEList_Final$counts)  # To see the numbers of remaining genes and samples
  boxplot(Tissue_DGEList_Final$samples$norm.factors)  # To show the range of scaling factors for Tissue libraries
  plotMDS(Tissue_DGEList_Final)  # Plot the sample distance with leading biological coefficient of variation (BCV)
  
  Tissue_DGEList_Final <- estimateDisp(Tissue_DGEList_Final)  # Estimate the dispersion
  Tissue_DGEList_Final$common.dispersion  # The common dispersion
  plotBCV(Tissue_DGEList_Final)
  
  Tissue_fit <- glmFit(Tissue_DGEList_Final)  # Test for DE genes
  Tissue_lrt <- glmLRT(Tissue_fit)
  FDR_Result <- topTags(Tissue_lrt, n = length(Tissue_lrt$table$PValue), adjust.method = "fdr", p.value = 0.05)
  write.table(FDR_Result$table, file = paste0("Sex_biased_Analysis/FR-Male_vs_Female/", Target_Tissue, "_SB_Gene.txt", sep = ""), quote = FALSE, sep = "\t")
}