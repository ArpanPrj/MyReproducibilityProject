#! /bin/bash

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories       

  ## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scratch/aubaxp004/Reproducibility/Fusarium_bataticola/RawData 
WD=/scratch/aubaxp004/Reproducibility/Fusarium_bataticola				
RDQ=RawDataQuality
 
##  make the directories in SCRATCH for holding the raw data 
## -p tells it to make any upper level directories that are not there. Notice how this will also make the WD.
mkdir -p ${DD}
## move to the Data Directory
cd ${DD}

##########  Download data files from NCBI: SRA using the Run IDs
  ### from SRA use the SRA tool kit - see NCBI website https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/
	## this downloads the SRA file and converts to fastq
	## -F 	Defline contains only original sequence name.
	## -I 	Append read id after spot id as 'accession.spot.readid' on defline.
	## splits the files into R1 and R2 (forward reads, reverse reads)

## These samples are transcriptome of Fusarium bataticola submitted to NCBI. This experiment has  3 samples on zucchini extract treated, 5 samples on sweet potato extract treated Fusarium RNA.

vdb-config --interactive

fastq-dump -F --split-files SRR10229809
fastq-dump -F --split-files SRR10229810
fastq-dump -F --split-files SRR10229811 
fastq-dump -F --split-files SRR10229812
fastq-dump -F --split-files SRR10229813
fastq-dump -F --split-files SRR10229814
fastq-dump -F --split-files SRR10229815
fastq-dump -F --split-files SRR10229816

##### Extra ####
## If you are downloading data from a sequencing company instead of NCBI, using wget for example, then calculate the md5sum values of all the files in the folder (./*), and read into a text file.
## then you can compare the values in this file with the ones provided by the company.
#md5sum ./* > md5sum.txt

############## FASTQC to assess quality of the sequence data
## FastQC: run on each of the data files that have 'All' to check the quality of the data
## The output from this analysis is a folder of results and a zipped file of results and a .html file for each sample
mkdir ${WD}/${RDQ}
fastqc *.fastq --outdir=${WD}/${RDQ}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
cd ${WD}/${RDQ}
tar cvzf ${RDQ}.tar.gz  ${WD}/${RDQ}/*
## when finished use scp or rsync to bring the tarballed .gz results file to your computer and open the .html file to evaluate the quality of your raw data.
