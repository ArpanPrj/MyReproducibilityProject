# RNA-Seq Data Analysis Project
This repository contains the complete pipeline and associated scripts for RNA-seq data analysis performed as part of a research project investigating transcriptional changes in *Fusarium* species response to specific biological conditions. The workflow spans from raw sequencing data preprocessing to differential gene expression analysis and downstream functional annotation.

**Key steps in this project include**

- [Data retrival and quality assessment](0_1_download_QualityCheck.sh)

- [Trimming of raw reads and post trim quality check](2_cleanTrimmomatic_QualityFastQC.sh)

- [Alignment to reference genome and transcript quantification](3_mapper_hisat2.sh)

- Differential expression analysis

- Functional enrichment and visualization

The repository aims to ensure reproducibility and transparency in analysis. All code is version-controlled, documented, and organized to allow others to replicate or adapt the workflow for their own studies.


## Data dowloading and quality checking
This step involves downloading RNA-seq data directly from the NCBI Sequence Read Archive (SRA). The data comes from a specific ........... experiment that includes samples ............. The SRA toolkit is used to retrieve the raw sequencing data using unique run accession IDs. Each downloaded sample is automatically split into forward and reverse reads, reflecting the paired-end nature of the sequencing setup. This ensures the data is correctly formatted for subsequent alignment and analysis steps.

Once the sequencing data is downloaded, the next step is to evaluate its quality using FastQC. FastQC generates a detailed report for each sample, highlighting important metrics such as per-base sequence quality, GC content, sequence duplication levels, and potential adapter contamination. These metrics help identify any issues in the raw data that might affect the reliability of downstream analyses. The resulting reports are packaged into a compressed archive for easy transfer and local inspection.

## Trimming of raw reads and post trim quality check
This part of the RNA-seq pipeline focuses on improving the quality of raw sequencing reads by removing unwanted adapter sequences and filtering out low-quality bases. The script uses Trimmomatic to trim adapters, crop low-quality regions, and discard short reads, resulting in cleaned paired and unpaired FASTQ files.

Following trimming, FastQC is used to assess the quality of the cleaned reads. This ensures that the trimming process was effective and that the data is suitable for accurate downstream analysis. The FastQC results are compiled into an archive file, making it easy to transfer and review quality metrics such as per-base quality scores, adapter contamination, and sequence length distribution.

## Alignment to reference genome and transcript quantification
This stage of the RNA-seq analysis pipeline involves aligning cleaned sequencing reads to a reference genome and quantifying gene and transcript expression levels using a combination of established bioinformatics tools.

The process begins with the preparation of the reference genome and its annotation using gffread, which converts the annotation file into a format compatible with downstream tools. The reference genome is then indexed with HISAT2, a fast and sensitive splice-aware aligner designed specifically for RNA-seq data. Cleaned paired-end reads are mapped to this indexed genome using HISAT2, allowing for accurate alignment across exon-exon junctions.

Post-alignment, SAMtools is used to convert, sort, and index the resulting alignment files, as well as to generate basic mapping statistics. These alignments are then passed to StringTie, a transcript assembler and quantification tool, which calculates expression levels using a reference-guided approach. StringTie outputs are structured in a way that is compatible with downstream visualization and statistical tools.

Finally, prepDE.py, a utility script provided by the StringTie suite, is used to compile all individual expression estimates into two comprehensive count matrices: one at the gene level and one at the transcript level. These matrices are essential inputs for differential expression analysis in tools such as DESeq2, edgeR, or Ballgown.
