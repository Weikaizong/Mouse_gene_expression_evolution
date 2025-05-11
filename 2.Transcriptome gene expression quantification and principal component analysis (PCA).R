###Gene expression quantification and principal component analysis (PCA) of the global transcriptome.###

setwd('D:/R/New_House_mouse_Transcriptome')
All_Count<-read.table("Coding-lncRNA_gene_Male_FeatureCounts.count",sep="\t",header = T)
rownames(All_Count)<-All_Count$Geneid
All_Count<-All_Count[,6:ncol(All_Count)]
All_Count<-as.matrix(All_Count)
Gene_Len<-All_Count[,1]
All_Count<-All_Count[,-1]

### 1. Expression quantification of combined protein-coding and long non-coding RNA (lncRNA) genes using edgeR, incorporating low counts filtering and normalization with the TMM method.

library("edgeR") #v3.30.3
All_DGEList<-DGEList(counts = All_Count,genes = Gene_Len, remove.zeros = F) #remove rows that have 0 total count.

All_DGEList_Keep <- filterByExpr(All_DGEList,min.count = 0, min.total.count=0)
All_DGEList_Final<-All_DGEList[All_DGEList_Keep, , keep.lib.sizes=FALSE]
All_DGEList_Final<-calcNormFactors(All_DGEList_Final) ##Normalize libary sizes with TMM
boxplot(All_DGEList_Final$samples$norm.factors) #To show the range of scaling factors for all the libraries
dim(All_DGEList_Final$counts)
##average count to be added to each observation to avoid taking log of zero.
All_FPKM_Matrix<-rpkm(All_DGEList_Final, gene.length = All_DGEList_Final$genes$genes, normalized.lib.sizes = TRUE,log = TRUE) 
All_FPKM_Matrix<-t(All_FPKM_Matrix)
write.table(All_FPKM_Matrix,file = "Coding-lncRNA_Male_FeatureCounts_8_Tissues_FPKM.txt", quote = F, sep = "\t")




### 2. Principal component analysis (PCA) of combined protein-coding and long non-coding RNA (lncRNA) gene expression profiles in the House mouse, revealing global clustering patterns across all organs.

library(ggfortify)# ‘0.4.16’
#Add the "Sample_ID" element in the header before the following code
All_Gene_Expression_Matrix <- read.table("Coding-lncRNA_Male_FeatureCounts_8_Tissues_FPKM.txt", header = T)
sample_Annotation <- read.table("All_Male_Sample_Annotation_FR_KA_SP_8_Tissues.tsv", header = T)
All_Gene_Expression_Data <- merge.data.frame(sample_Annotation, All_Gene_Expression_Matrix, by.x = "Sample_ID", by.y = "Sample_ID")
rownames(All_Gene_Expression_Data) <- All_Gene_Expression_Data$Sample_ID
All_Gene_Expression_Value <- All_Gene_Expression_Data[, 5:ncol(All_Gene_Expression_Data)]
alpha_value <- 1.3
size_point <- 5
size_label <- 10
tissues <- c("Brain", "Epididymus","Heart", "Kidney", "Liver", "Mammary", "Testis",   "Vas_Deferens") 
tissue_colors <- c("#E1CBA0", "#3182bd", "#FA5959", "#C4D7B7", "#DC7F61", "#756bb1", "#B6CBDF", "#43a2ca")
populations <- c("M.m.dom_FR", "M.m.mus_KA", "M.spretus_SP" )
shapes <- c( 1, 2, 3)
#------------------------ Plot with the FPKM data with log2 transformation ---------------------------------------------------------
######Panel A (all tissues) 
autoplot(prcomp(All_Gene_Expression_Value), data = All_Gene_Expression_Data, colour = 'Tissue', shape = 'Population', size = size_point, alpha = alpha_value) +
  theme(text = element_text(size = size_label,face = "bold"), 
        axis.text = element_text(size = 13), 
        axis.title = element_text(size = 13), 
        legend.key = element_rect(fill = "transparent", colour = "transparent"), 
        legend.position = 'none',
        legend.box = "vertical",
        #legend.margin = margin(3, 7, 3, 7),  
        legend.text = element_text(size =9),
        legend.title = element_blank(),
        #legend.background = element_rect(fill = "gray97", colour = "gray", size = 0.1, linetype = "solid"),
        panel.background = element_rect(fill = "gray100", colour = "black", size = 0.5, linetype = "solid")) +
  scale_color_manual(values = tissue_colors) +  
  scale_shape_manual(values = shapes) +  
  guides(title = NA, col = guide_legend(nrow = 8)) 

ggsave("Coding-lncRNA_gene_Male_FR-KA-SP_8_Tissues_pca.pdf", width = 12, height = 11, units = "cm", dpi = 500)
ggsave("Coding-lncRNA_gene_Male_FR-KA-SP_8_Tissues_pca.png", width = 12, height = 11, units = "cm", dpi = 500) 


###Gene expression quantification and principal component analysis (PCA) of the organ-specific transcriptome.###

folder_path <- "Coding_Tissue_Male_count"  
conditions <- c("Bra", "Hea", "Kid", "Liv", "Tes", "Epi", "Mam", "Vas")

for (condition in conditions) {
  count_file <- paste0(condition, "_FeatureCounts_Fragment.count")
  count_file_path <- file.path(folder_path, count_file)
  All_Count <- read.table(count_file_path, sep = "\t", header = TRUE)
  rownames(All_Count) <- All_Count$Geneid
  All_Count <- All_Count[, 6:ncol(All_Count)]
  All_Count <- as.matrix(All_Count)
  Gene_Len <- All_Count[, 1]
  All_Count <- All_Count[, -1]
  All_DGEList <- DGEList(counts = All_Count, genes = Gene_Len, remove.zeros = FALSE)
  All_DGEList_Keep <- filterByExpr(All_DGEList, min.count = 0, min.total.count = 0)
  All_DGEList_Final <- All_DGEList[All_DGEList_Keep, , keep.lib.sizes = FALSE]
  All_DGEList_Final <- calcNormFactors(All_DGEList_Final)
  All_FPKM_Matrix <- rpkm(All_DGEList_Final, gene.length = All_DGEList_Final$genes$genes, normalized.lib.sizes = TRUE, log = TRUE) 
  All_FPKM_Matrix <- t(All_FPKM_Matrix)
  output_file <- paste0(condition, "_coding_FPKM_Matrix_STAR.txt")
  write.table(All_FPKM_Matrix, file = output_file, quote = FALSE, sep = "\t")
}
file_names <- c("Bra_coding_FPKM_Matrix_STAR.txt",
                "Hea_coding_FPKM_Matrix_STAR.txt",
                "Kid_coding_FPKM_Matrix_STAR.txt",
                "Liv_coding_FPKM_Matrix_STAR.txt",
                "Tes_coding_FPKM_Matrix_STAR.txt",
                "Epi_coding_FPKM_Matrix_STAR.txt",
                "Mam_coding_FPKM_Matrix_STAR.txt",
                "Vas_coding_FPKM_Matrix_STAR.txt")
annotation_files <- c("Brain_annotation.tsv",
                      "Heart_annotation.tsv",
                      "Kidney_annotation.tsv",
                      "Liver_annotation.tsv",
                      "Testis_annotation.tsv",
                      "Epididymus_annotation.tsv",
                      "Mammary_annotation.tsv",
                      "Vas_Deferens_annotation.tsv")
tissue_colors <- c(
  "Brain" = "#E1CBA0",
  "Epididymus" = "#3182bd",
  "Heart" = "#FA5959",
  "Kidney" = "#C4D7B7",
  "Liver" = "#DC7F61",
  "Mammary" = "#756bb1",
  "Testis" = "#B6CBDF",
  "Vas_Deferens" = "#43a2ca"
)
populations <- c("M.m.dom_FR", "M.m.mus_KA", "M.spretus_SP")
shapes <- c(1, 2, 3)
for (i in seq_along(file_names)) {
  All_Gene_Expression_Matrix <- read.table(file_names[i], header = TRUE)
  sample_Annotation <- read.table(annotation_files[i], header = TRUE)
  All_Gene_Expression_Data <- merge.data.frame(sample_Annotation, All_Gene_Expression_Matrix, by.x = "Sample_ID", by.y = "Sample_ID")
  rownames(All_Gene_Expression_Data) <- All_Gene_Expression_Data$Sample_ID
  All_Gene_Expression_Value <- All_Gene_Expression_Data[, 5:ncol(All_Gene_Expression_Data)]
  alpha_value <- 1
  size_point <- 13
  size_label <- 19
  shapes <- c(1, 2, 3)
  plot <- autoplot(prcomp(All_Gene_Expression_Value), data = All_Gene_Expression_Data, 
                   colour = 'Tissue', shape = 'Population', size = size_point, alpha = alpha_value) +
    geom_point(aes(shape = Population, colour = Tissue), size = size_point, stroke = 3) +
    theme(plot.margin = margin(t = 11, r = 11, b = 11, l = 11))+ 
    theme(
      text = element_text(size = size_label,face = "bold"), 
      legend.key = element_rect(fill = "transparent", colour = "transparent"), 
      legend.position = "none", 
      legend.box = "vertical",
      legend.margin = margin(),
      legend.text = element_text(size = 30), 
      legend.title = element_blank(),
      panel.background = element_rect(fill = "white", colour = "black", size = 1.7, linetype = "solid"),
      axis.text.x = element_text(size = 30),  
      axis.text.y = element_text(size = 30),  
      axis.title.x = element_text(size = 34), 
      axis.title.y = element_text(size = 34) ) + 
    scale_color_manual(values = tissue_colors) +
    scale_shape_manual(values = shapes) +
    guides(title = NA, col = guide_legend(nrow = 5))
  plot_file_name <- paste(substr(file_names[i], 1, 3), "_PCA_Plot.png", sep = "")
  pdf_file_name <- paste(substr(file_names[i], 1, 3), "_PCA_Plot.pdf", sep = "")
  
  ggsave(plot_file_name, plot, width = 9, height = 8.5, units = "in", dpi = 1500)
  ggsave(pdf_file_name, plot, width = 9, height = 8.5, units = "in", dpi = 1500)
}