# RNA-Seq Data Analysis Project
This repository contains the complete pipeline and associated scripts for RNA-seq data analysis performed as part of a research project investigating transcriptional changes in *Fusarium* species response to specific biological conditions. The workflow spans from raw sequencing data preprocessing to differential gene expression analysis and downstream functional annotation.

**Key steps in this project include**
_ Data retrival and quality assessment

- Trimming of raw reads and post trim quality check

- Alignment to reference genome

- Transcript quantification

- Differential expression analysis

- Functional enrichment and visualization

The repository aims to ensure reproducibility and transparency in analysis. All code is version-controlled, documented, and organized to allow others to replicate or adapt the workflow for their own studies.


## Data dowloading and quality checking
### Data Acquisition
This step involves downloading RNA-seq data directly from the NCBI Sequence Read Archive (SRA). The data comes from a specific ........... experiment that includes samples ............. The SRA toolkit is used to retrieve the raw sequencing data using unique run accession IDs. Each downloaded sample is automatically split into forward and reverse reads, reflecting the paired-end nature of the sequencing setup. This ensures the data is correctly formatted for subsequent alignment and analysis steps.

### Quality Assessment
Once the sequencing data is downloaded, the next step is to evaluate its quality using FastQC. FastQC generates a detailed report for each sample, highlighting important metrics such as per-base sequence quality, GC content, sequence duplication levels, and potential adapter contamination. These metrics help identify any issues in the raw data that might affect the reliability of downstream analyses. The resulting reports are packaged into a compressed archive for easy transfer and local inspection.
