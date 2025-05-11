#!/bin/bash

# Set environment variables and tool paths
export LD_LIBRARY_PATH=/public/home/others/zhangwenyu_lab/wangpengfei/software/mummer/lib

# fastp configuration
export fastp=/public/home/others/zhangwenyu_lab/Common_Files/01.Software/fastp-0.20.0/fastp
export raw_data=/public/home/others/zhangwenyu_lab/Common_Files/02.Data/House_Mouse/sex_bias/PRJEB50011_DE_female
export cleandata=/public/home/others/zhangwenyu_lab/Common_Files/02.Data/House_Mouse/sex_bias/DE_female_align/01.clean_data

# STAR configuration
export STAR=/public/home/others/zhangwenyu_lab/Common_Files/01.Software/STAR-2.7.10a_alpha_220601/source/STAR
export GenomeDir=/public/home/others/zhangwenyu_lab/Common_Files/02.Data/House_Mouse/Genome/GRCm39_genomeDir_STAR_Long150
export bam=/public/home/others/zhangwenyu_lab/Common_Files/02.Data/House_Mouse/sex_bias/DE_female_align/02.bam

# featureCounts configuration
export featureCounts=/public/home/others/zhangwenyu_lab/Common_Files/01.Software/subread-1.6.3-Linux-x86_64/bin/featureCounts

# 1. Perform quality control with fastp
$fastp \
  -i $raw_data/ERR9112877_1.fastq.gz \
  -o $cleandata/ERR9112877_QC_1.fastq \
  -I $raw_data/ERR9112877_2.fastq.gz \
  -O $cleandata/ERR9112877_QC_2.fastq \
  --cut_tail --average_qual 20 --length_required 50 \
  -j ERR9112877.json -h ERR9112877.html -w 5
echo "fastp is Done!"

# 2. Perform alignment with STAR
$STAR \
  --runThreadN 5 \
  --genomeDir $GenomeDir \
  --readFilesIn $cleandata/ERR9112877_QC_1.fastq $cleandata/ERR9112877_QC_2.fastq \
  --outFileNamePrefix $bam/ERR9112877 \
  --outFilterMismatchNmax 30 --scoreDelOpen -1 --scoreDelBase -1 \
  --scoreInsOpen -1 --scoreInsBase -1 --seedSearchStartLmax 25 \
  --winAnchorMultimapNmax 100 \
  --outSAMtype BAM SortedByCoordinate \
  --twopassMode Basic --twopass1readsN -1
echo "Align is Done!"

# 3. Calculate gene expression using featureCounts
$featureCounts \
  -p -T 2 -s 2 -t exon -g gene_id \
  -a /public/home/others/zhangwenyu_lab/wangpengfei/Genome/Mus_musculus.GRCm39.110.gtf \
  -o /public/home/others/zhangwenyu_lab/Common_Files/02.Data/House_Mouse/sex_bias/DE_female_align/03.count/ERR9112877_FeatureCounts.count \
  $bam/ERR9112877Aligned.sortedByCoord.out.bam
echo "Count is Done!"
