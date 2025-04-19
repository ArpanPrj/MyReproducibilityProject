#! /bin/bash

########## Load Modules
source /apps/profiles/modules_asax.sh.dyn
module load sra
module load fastqc/0.10.1

##########  Define variables and make directories
## Replace the numbers in the brackets with Your specific information
  ## make variable for your ASC ID so the directories are automatically made in YOUR directory
MyID=[1]          ## Example: MyID=aubtss

  ## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=[2]   			## Example: DD=/scratch/${MyID}/PracticeRNAseq/RawData
WD=[3]				## Example: WD=/scratch/${MyID}/PracticeRNAseq
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

## These samples are from Bioproject PRJNA437447. An experiment on Daphnia pulex, 5 samples on ad lib feed, 5 samples on caloric restriction diet
## https://www.ncbi.nlm.nih.gov/bioproject?LinkName=sra_bioproject&from_uid=5206312
## For class only do the 6 that you are assigned, delete the other 4 from this list

vdb-config --interactive

fastq-dump -F --split-files SRR6819014
fastq-dump -F --split-files SRR6819015
fastq-dump -F --split-files SRR6819016
fastq-dump -F --split-files SRR6819017
fastq-dump -F --split-files SRR6819018
fastq-dump -F --split-files SRR6819019
fastq-dump -F --split-files SRR6819020
fastq-dump -F --split-files SRR6819021
fastq-dump -F --split-files SRR6819022
fastq-dump -F --split-files SRR6819023

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
