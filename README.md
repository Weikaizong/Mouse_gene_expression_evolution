This repository contains a pipeline for analyzing transcriptomic data in different sexes and populations of the house mouse (Mus musculus). The pipeline includes preprocessing of raw data, gene expression quantification, principal component analysis (PCA), transcriptomic divergence rate analysis, differentially expressed gene analysis, and sex-biased gene analysis.

File Structure

1.Preprocessing of raw data.sh: Preprocesses raw sequencing data, including quality control, alignment, and gene counting.
2.Transcriptome gene expression quantification and principal component analysis (PCA).R : Performs gene expression quantification of combined protein-coding and long non-coding RNA(lncRNA) genes using edgeR with low counts filtering and normalization via TMM method, and conducts PCA of gene expression profiles in house mouse organs to reveal clustering patterns.
3.Transcriptomic divergence rates and Differentially expressed gene analysis.R : Analyzes transcriptomic divergence rates using 1−Spearman’s ρ and carries out differentially expressed gene analysis using edgeR.
4.Sex Biased genes analysis in Different Populations of House Mouse.R : Conducts sex biased genes analysis in different populations of house mouse using edgeR.

Dependencies

The following tools are required for data preprocessing:fastp、STAR、featureCounts
The following R packages are required to run the R scripts: edgeR、ggfortify、readr、DESeq2
You can install these packages using the following commands: 
	install.packages("edgeR")
	install.packages("ggfortify")
	install.packages("readr")
	install.packages("DESeq2")

Usage

Preprocessing of Raw Data
Run the shell script to perform quality control, alignment, and gene counting: bash 1.Preprocessing of raw data.sh

Transcriptome Gene Expression Quantification and PCA
1.Set the working directory to the location of your count data file (Coding-lncRNA_gene_Male_FeatureCounts.count) and sample annotation file (All_Male_Sample_Annotation_FR_KA_SP_8_Tissues.tsv).
2.Run the script 1.Transcriptome gene expression quantification and principal component analysis (PCA).R.

Transcriptomic Divergence Rates and Differentially Expressed Gene Analysis
1.Set the working directory to the location of your FPKM means data file (All_Male_8_tissues_FPKM_means.txt) and count data files for different tissues.
2.Run the script 2.Transcriptomic divergence rates and Differentially expressed gene analysis.R.

Sex Biased Genes Analysis
1.Set the working directory to the location of your protein coding gene list file (GFF110_Protein_Coding_Gene.txt) and count data files for different tissues and populations.
2.Run the script 3.Sex Biased genes analysis in Different Populations of House Mouse.R.

Results
The scripts will generate the following results:
-Quality-Controlled FASTQ Files: Cleaned FASTQ files after quality control using fastp.
-Aligned BAM Files: BAM files generated after alignment using STAR.
-Gene Count Files: Count files generated after gene expression quantification using featureCounts.
-FPKM matrices : Log-transformed FPKM values of genes in different tissues and populations.
-PCA plots : Visualizations of the principal components of gene expression data, showing the clustering patterns of different tissues and populations.
-Spearman correlation results : Transcriptomic divergence rates calculated using 1−Spearman’s ρ.
-Differentially expressed gene lists : Lists of genes that are differentially expressed between different groups.
-Sex biased gene lists : Lists of genes that show sex-biased expression in different tissues and populations.