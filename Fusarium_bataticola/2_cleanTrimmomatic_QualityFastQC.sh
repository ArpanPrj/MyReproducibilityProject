#! /bin/bash

## Purpose: The purpose of this script is to trim sequencing adapters and low quality regions from the read data.
## Input Data: Raw R1 & R2 reads (FASTQ); Adapter sequences to remove (FASTA)
## Output Data: Trimmed R1 & R2 paired and unpaired reads (FASTQ)
## More Information: http://www.usadellab.org/cms/?page=trimmomatic

# Modules
	#  load the module
source /apps/profiles/modules_asax.sh.dyn
module load trimmomatic/0.39
module load fastqc/0.10.1

##########  Define variables and make directories
  ## Make variable that represents YOUR working directory(WD) in scratch, your Raw data directory (DD) and the pre or postcleaned status (CS).
DD=/scratch/aubaxp004/Reproducibility/Fusarium_bataticola/RawData 
WD=/scratch/aubaxp004/Reproducibility/Fusarium_bataticola                              
CD=/scratch/aubaxp004/Reproducibility/Fusarium_bataticola/CleanData
PCQ=PostCleanQuality
adapters=AdaptersToTrim_All.fa  ## This is a fasta file that has a list of adapters commonly used in NGS sequencing. 
				## In the future, for your data, you will likely need to edit this for other projects based on how your libraries 
				## were made to search for the correct adapters for your project
				
## make the directories to hold the Cleaned Data files, and the directory to hold the results for assessing quality of the cleaned data.
mkdir ${CD}
mkdir ${WD}/${PCQ}

################ Trimmomatic ###################################
## Move to Raw Data Directory
cd ${DD}

### Make list of file names to Trim
        ## this line is a set of piped (|) commands
        ## ls means make a list, 
        ## grep means grab all the file names that end in ".fastq", 
        ## cut that name into elements every where you see "_" and keep the first element (-f 1)
        ## sort the list and keep only the unique names and put it into a file named "list"
ls | grep ".fastq" |cut -d "_" -f 1 | sort | uniq > list

### Copy over the list of Sequencing Adapters that we want Trimmomatic to look for (along with its default adapters)
cp /scratch/aubaxp004/Illumina/AdaptersToTrim_All.fa . 

### Run a while loop to process through the names in the list and Trim them with the Trimmomatic Code
while read i
do

        ### Run Trimmomatic in paired end (PE) mode with 6 threads using phred 33 quality score format. 
        ## STOP & DISCUSS: Check out the trimmomatic documentation to understand the parameters in line 77
	#From Brittany
                 #/apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar 
        java -jar /apps/x86-64/apps/spack_0.19.1/spack/opt/spack/linux-rocky8-zen3/gcc-11.3.0/trimmomatic-0.39-iu723m2xenra563gozbob6ansjnxmnfp/bin/trimmomatic-0.39.jar   \
	PE -threads 6 -phred33 \
        "$i"_1.fastq "$i"_2.fastq  \
        ${CD}/"$i"_1_paired.fastq ${CD}/"$i"_1_unpaired.fastq  ${CD}/"$i"_2_paired.fastq ${CD}/"$i"_2_unpaired.fastq \
        ILLUMINACLIP:AdaptersToTrim_All.fa:2:35:10 HEADCROP:10 LEADING:30 TRAILING:30 SLIDINGWINDOW:6:30 MINLEN:36
        
                ## Trim read for quality when quality drops below Q30 and remove sequences shorter than 36 bp
                ## PE for paired end phred-score-type  R1-Infile   R2-Infile  R1-Paired-outfile R1-unpaired-outfile R-Paired-outfile R2-unpaired-outfile  Trimming paramenter
                ## MINLEN:<length> #length: Specifies the minimum length of reads to be kept.
                ## SLIDINGWINDOW:<windowSize>:<requiredQuality>  #windowSize: specifies the number of bases to average across  
                ## requiredQuality: specifies the average quality required.

	############## FASTQC to assess quality of the Cleaned sequence data
	## FastQC: run on each of the data files that have 'All' to check the quality of the data
	## The output from this analysis is a folder of results and a zipped file of results

fastqc ${CD}/"$i"_1_paired.fastq --outdir=${WD}/${PCQ}
fastqc ${CD}/"$i"_2_paired.fastq --outdir=${WD}/${PCQ}

done<list			# This is the end of the loop

#########################  Now compress your results files from the Quality Assessment by FastQC 
## move to the directory with the cleaned data
cd ${WD}/${PCQ}

#######  Tarball the directory containing the FASTQC results so we can easily bring it back to our computer to evaluate.
tar cvzf ${PCQ}.tar.gz ${WD}/${PCQ}/*

## when finished use scp or rsync to bring the .gz file to your computer and open the .html file to evaluate the quality of the data.
