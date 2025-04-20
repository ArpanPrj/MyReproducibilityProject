#  RNA-Seq Data Analysis Project

Welcome to the RNA-seq analysis pipeline designed to unravel how *Fusarium* species tweak their gene expression when they encounter host-derived biological extracts. This study dives into two specific fungal pathogens:

- *Fusarium cucurbiticola*, treated with zucchini extract  
- *Fusarium bataticola*, treated with sweet potato extract

Each fungus gets its own little sandbox — the entire analysis is neatly compartmentalized into two separate folders: `Fusarium_cucurbiticola/` and `Fusarium_bataticola/`, so there's no cross-contamination or confusion.

---

## Workflow Overview

This RNA-seq pipeline is built as a stepwise and reproducible workflow — from raw read downloads all the way to functional interpretations of differentially expressed genes.

---

### Data Retrieval and Quality Check

We begin by getting our hands on some raw data. Using the **SRA Toolkit (v3.0.0)**, we pull down paired-end RNA-seq reads directly from the NCBI SRA. The reads are automatically split into forward and reverse files for easy handling.

Next up, we run **FastQC (v0.11.9)**, which acts like a preliminary health check for your reads. It tells us about:
- Base-wise quality scores
- Presence of adapter sequences
- Sequence duplication
- GC content quirks

All FastQC reports are stored and can be viewed later to see how good (or bad) your raw data looks.

Scripts:
- [`Fusarium_cucurbiticola/0_1_download_QualityCheck.sh`](./Fusarium_cucurbiticola/0_1_download_QualityCheck.sh)  
- [`Fusarium_bataticola/0_1_download_QualityCheck.sh`](./Fusarium_bataticola/0_1_download_QualityCheck.sh)

---

### Adapter Trimming and Post-Trim QC

Once we know what's wrong with the raw reads, we clean them up using **Trimmomatic (v0.39)**. It clips off the pesky Illumina adapter sequences, shaves off low-quality ends, and kicks out reads that are too short to be useful.

We then rerun **FastQC** on these trimmed reads to confirm that we've improved things. If the reports look greener, we move on with confidence.

Scripts:
- [`Fusarium_cucurbiticola/2_cleanTrimmomatic_QualityFastQC.sh`](./Fusarium_cucurbiticola/2_cleanTrimmomatic_QualityFastQC.sh)  
- [`Fusarium_bataticola/2_cleanTrimmomatic_QualityFastQC.sh`](./Fusarium_bataticola/2_cleanTrimmomatic_QualityFastQC.sh)

---

### Read Mapping and Expression Quantification

Now the fun begins — we align the clean reads to their respective reference genomes. Here's how this step unfolds:

- We first process the GFF/GTF annotation files using **gffread (v0.12.7)** to make sure they play nicely with downstream tools.
- The reference genome is indexed using **HISAT2 (v2.2.1)** — a spliced aligner that smartly handles reads spanning exon-exon junctions.
- Reads are aligned to the genome, producing SAM files, which are promptly converted to sorted and indexed BAMs using **SAMtools (v1.15.1)**.
- Then comes **StringTie (v2.2.1)**, which assembles transcripts and estimates expression levels. It works in a reference-guided mode, keeping things tidy and biologically accurate.
- Finally, we use the **prepDE.py** script bundled with StringTie to extract gene- and transcript-level count matrices — all ready for DE analysis.

Scripts:
- [`Fusarium_cucurbiticola/3_mapper_hisat2.sh`](./Fusarium_cucurbiticola/3_mapper_hisat2.sh)  
- [`Fusarium_bataticola/3_mapper_hisat2.sh`](./Fusarium_bataticola/3_mapper_hisat2.sh)

---

### Differential Expression Analysis

With count data in hand, we jump into R to perform statistical sleuthing using **DESeq2 (v1.38.3)**. Here's what happens:

- Raw counts are normalized to account for differences in library size and sequencing depth.
- We run statistical tests to detect genes that are significantly up- or down-regulated between treatment conditions.
- MA plots and volcano plots are generated to help visualize the expression landscape and spot key players.
- DEGs are filtered based on log2 fold-change and adjusted p-values (because false discoveries are a party we don't want to crash).

Scripts: *Coming soon*

---

### Functional Enrichment and Biological Insights

We wrap up the analysis by asking the all-important question: **What do these DEGs actually do?**

Using R packages like `clusterProfiler`, we run GO and KEGG enrichment analysis to identify overrepresented functions and pathways.

These results are visualized using dot plots, barplots, and curated pathway diagrams, giving us a systems-level understanding of how the fungi respond to their host extracts.

Scripts: *Coming soon*

---

## Software and Dependencies

Everything in this pipeline was run on a Unix-based HPC setup. Here’s the toolkit used and their respective versions:

| Tool         | Version  | Description |
|--------------|----------|-------------|
| SRA Toolkit  | 3.0.0    | Fetch raw sequencing data from NCBI SRA |
| FastQC       | 0.11.9   | Assess sequence quality before and after trimming |
| Trimmomatic  | 0.39     | Remove adapters and low-quality reads |
| HISAT2       | 2.2.1    | Splice-aware aligner for RNA-seq reads |
| SAMtools     | 1.15.1   | BAM/SAM manipulation toolkit |
| StringTie    | 2.2.1    | Transcript assembler and quantifier |
| prepDE.py    | bundled with StringTie | Generates count matrices for DESeq2 |
| gffread      | 0.12.7   | Annotation fixer and converter |
| DESeq2       | 1.38.3   | Differential gene expression analysis |
| R            | 4.2.2    | Statistical computing and visualization |
| Python       | 3.9      | Used for scripting various steps |



---

## Notes

- The project is built with reproducibility in mind, adhering to FAIR principles.
- It’s modular, version-controlled, and easy to extend.
- Each script is prefixed and named according to the step number, so they’re always in logical order.

---

## Citation and Contact

If this pipeline helps your research, please cite this repository in your work. For questions, bugs, or collaborations, feel free to reach out to the maintainer.
