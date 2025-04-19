# RNA-Seq Data Analysis Project

This repository contains the complete pipeline and associated scripts for RNA-seq data analysis performed as part of a research project investigating transcriptional changes in *Fusarium* species in response to host-specific biological extracts. The workflow spans from raw sequencing data preprocessing to differential gene expression analysis and downstream functional annotation.

The project is divided into two independent branches of analysis:

- *Fusarium cucurbiticola* treated with zucchini extract  
- *Fusarium bataticola* treated with sweet potato extract  

All code and outputs are organized in two subdirectories: `Fusarium_cucurbiticola` and `Fusarium_bataticola`, each containing a complete pipeline.

---

## Key Steps in the Workflow

- Data retrieval and quality assessment  
- Trimming of raw reads and post-trim quality check  
- Alignment to reference genome and transcript quantification  
- Differential expression analysis  
- Functional enrichment and visualization  

Each step is documented through version-controlled scripts to ensure reproducibility and transparency.

---

## Folder Structure and Script Links

### Fusarium cucurbiticola

Directory: [`Fusarium_cucurbiticola`](./Fusarium_cucurbiticola)

- [0_1_download_QualityCheck.sh](./Fusarium_cucurbiticola/0_1_download_QualityCheck.sh): Download RNA-seq data and perform initial quality check using FastQC  
- [2_cleanTrimmomatic_QualityFastQC.sh](./Fusarium_cucurbiticola/2_cleanTrimmomatic_QualityFastQC.sh): Trim raw reads with Trimmomatic and perform post-trim QC  
- [3_mapper_hisat2.sh](./Fusarium_cucurbiticola/3_mapper_hisat2.sh): Align reads to reference genome using HISAT2 and quantify transcripts with StringTie  

---

### Fusarium bataticola

Directory: [`Fusarium_bataticola`](./Fusarium_bataticola)

- [0_1_download_QualityCheck.sh](./Fusarium_bataticola/0_1_download_QualityCheck.sh): Download RNA-seq data and perform initial quality check using FastQC  
- [2_cleanTrimmomatic_QualityFastQC.sh](./Fusarium_bataticola/2_cleanTrimmomatic_QualityFastQC.sh): Trim raw reads with Trimmomatic and perform post-trim QC  
- [3_mapper_hisat2.sh](./Fusarium_bataticola/3_mapper_hisat2.sh): Align reads to reference genome using HISAT2 and quantify transcripts with StringTie  

---

## Software and Dependencies

The analysis pipeline was developed and tested using the following tools and versions:

| Tool         | Version  | Description |
|--------------|----------|-------------|
| SRA Toolkit  | 3.0.0    | Downloading raw RNA-seq data from NCBI SRA |
| FastQC       | 0.11.9   | Quality check of raw and trimmed reads |
| Trimmomatic  | 0.39     | Adapter trimming and quality filtering |
| HISAT2       | 2.2.1    | Splice-aware alignment of reads to the reference genome |
| SAMtools     | 1.15.1   | Conversion, sorting, indexing, and stats of alignment files |
| StringTie    | 2.2.1    | Transcript assembly and expression quantification |
| prepDE.py    | (from StringTie 2.2.1) | Generation of count matrices for differential expression |
| gffread      | 0.12.7   | GFF/GTF conversion for StringTie compatibility |
| DESeq2       | 1.38.3   | (R package) Differential gene expression analysis |
| R            | 4.2.2    | Statistical computing and plotting |
| Python       | 3.9      | Scripting and workflow automation |

> All tools were run on a Unix-based high-performance computing environment.

---

## Notes

- The `prepDE.py` script from the StringTie suite is used to generate count matrices for DE analysis  
- All tools can be installed via `conda`, `apt`, or compiled from source  
- The project adheres to FAIR data principles for bioinformatics workflows  

---

## Citation

Please cite this repository if you use any part of the pipeline or code in your own work. For questions or collaboration, contact the repository maintainer.
