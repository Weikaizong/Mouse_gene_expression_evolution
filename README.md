The repository contains related datasets and major analysis methods for the manuscript "Full-length transcript sequencing traces brain isoform diversity in house mouse natural populations".


- **Datasets**
1. PacBio Iso-Seq data of one tested mouse individual, for the optimization of Iso-Seq library preparation protocols. The raw subread data can be found in European Nucleotide Archive (ENA; https://www.ebi.ac.uk/ena) under accession number: PRJEB54001. The primer sequences of Clontech and Lexogen are available at the ./Data folder of this repository.
2. PacBio Iso-Seq data of the 48 mice individuals reared under laboratory environment (main experiment). The raw subread data can be found in European Nucleotide Archive (ENA; https://www.ebi.ac.uk/ena) under accession number: PRJEB54000.
3. Illumina RNA-Seq data of the 48 mice individuals reared under laboratory environment (main experiment). The raw fastq data can be found in European Nucleotide Archive (ENA; https://www.ebi.ac.uk/ena) under accession number: PRJEB53988.
4. Public available Illumina RNA-Seq data on 13 tissues of C57BL/6 mouse inbred line. The dataset was retrieved from Lin S et al 2014 (PMID: 25413365).
5. Alignment bam files, GTF track data, and SNP VCF files are stored at ftp site of MPI: https://wwwuser.gwdg.de/~evolbio/evolgen/wildmouse/mouse_population_isoform/.


- **Methods**
1. Iso-Seq_Library_Protocol_Selection: Analysis code and description concerning the test on the performance of three PacBio Iso-Seq library preparation protocols. 
2. LongRead_Sequence_Processing: Analysis code and description concerning the processing of PacBio Iso-Seq dataset. 
3. ShortRead_Sequence_Processing: Analysis code and description concerning the processing of Illumina RNA-Seq dataset. 
