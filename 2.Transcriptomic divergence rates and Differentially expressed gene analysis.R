###Transcriptomic divergence rate analysis using 1−Spearman’s ρ.###

library(readr)
setwd('D:/R/New_House_mouse_Transcriptome/Spearman_correlations_8_Tissues/Protein_coding')
data <- read_delim('All_Female_8_tissues_FPKM_means.txt', delim = "\t")
prefixes <- c("Bra", "Hea", "Kid", "Liv", "Ova", 'Mam','Ute', 'Ovi')
results <- data.frame(Comparison = character(), Spearman_Correlation = numeric(), stringsAsFactors = FALSE)
for (prefix in prefixes) {
  dom_mean <- data[[paste0(prefix, "_DOM_mean")]]
  mus_mean <- data[[paste0(prefix, "_MUS_mean")]]
  dm_mean <- data[[paste0(prefix, "_DOM_MUS_mean")]]
  spr_mean <- data[[paste0(prefix, "_SPR_mean")]]

  spearman_dom_mus <- cor(dom_mean, mus_mean, method = "spearman", use = "complete.obs")
  spearman_dm_spr <- cor(dm_mean, spr_mean, method = "spearman", use = "complete.obs")
  results <- rbind(results, 
                   data.frame(Comparison = paste0(prefix, "_DOM_mean vs ", prefix, "_MUS_mean"), 
                              Spearman_Correlation = spearman_dom_mus),
                   data.frame(Comparison = paste0(prefix, "_DOM_MUS_mean vs ", prefix, "_SPR_mean"), 
                              Spearman_Correlation = spearman_dm_spr))
}
output_file_path <- 'Female_spearman_correlations_8_Tissues.txt'
write.table(results, file = output_file_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Results saved to", output_file_path, "\n")



###Differentially expressed (DE) gene analysis using edgeR.###
library("edgeR") 
library('DESeq2')
setwd('D:/R/New_House_mouse_Transcriptome')
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
##################################### Get DE genes with edgeR ###############################################
Tissue_ID <- c("Bra", "Hea", "Kid", "Liv", "Tes", "Mam", "Vas", "Epi")
for (i in 1:length(Tissue_ID)) {
  Target_Tissue <- Tissue_ID[i]
  Tissue_Count <- read.table(paste0("Coding_Tissue_Female_count/", Target_Tissue, "_FeatureCounts_Fragment.count", sep = ""), sep = "\t", header = TRUE)
  # rownames(Tissue_Count) <- Tissue_Count$Geneid
  rownames(Tissue_Count) <- Tissue_Count[, 1]  
  Tissue_Count <- Tissue_Count[, 7:ncol(Tissue_Count)]
  Tissue_Count <- as.matrix(Tissue_Count)
  Tissue_Count <- Tissue_Count[Coding_Gene_List, ]
  
  # 1) FR_VS_KA
  Tissue_Annotation <- read.table(paste0("DE_gene_Analysis/FR_VS_KA/", Target_Tissue, "_Sample_Annotation.tsv", sep = ""), header = TRUE, row.names = 1)
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
  write.table(FDR_Result$table, file = paste0("DE_gene_Analysis/FR_VS_KA/", Target_Tissue, "_DE_Gene.txt", sep = ""), quote = FALSE, sep = "\t")
  # 2) FR-KA_VS_SP
  Tissue_Annotation <- read.table(paste0("DE_gene_Analysis/FR-KA_VS_SP/", Target_Tissue, "_Sample_Annotation.tsv", sep = ""), header = TRUE, row.names = 1)
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
  write.table(FDR_Result$table, file = paste0("DE_gene_Analysis/FR-KA_VS_SP/", Target_Tissue, "_DE_Gene.txt", sep = ""), quote = FALSE, sep = "\t")
}
