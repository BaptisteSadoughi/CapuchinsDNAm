#!/bin/bash

# Note: It may be more convenient to first run trimgalore, and once completed run the remaining steps in one job.

# Load Required Software
module load mamba/latest
module load bowtie2-2.4.2-gcc-11.2.0
module load samtools-1.13-gcc-11.2.0
module load py-cutadapt-2.10-gcc-11.2.0

# Activate python environment
export PATH=$PATH:/..../FastQC:/..../bin:/..../TrimGalore-0.6.6:/..../bismark2/Bismark-0.24.0

# path to fastq files
fastq_path=/path/to/fastq

# Assign samples by referring to a manifest file listing the libraries.
sampleID=`sed -n ${SLURM_ARRAY_TASK_ID}p /path/to/fastq/manifest.txt`

# set working directory
cd /path/to/workingdirectory

genome_path=/path/to/cimit_wholegenome/cimit
out_path=./mapped
cov_out=./cov
trim_path=/path/to/trimmed_maps/trimmed
finalcov_out=./final_cov


## trim adaptors
trim_galore --paired -j 4 -o ${trim_path} \
 --basename ${sampleID} --gzip ${fastq_path}/${sampleID}_S1_L005_R1_001.fastq.gz ${fastq_path}/${sampleID}_S1_L005_R2_001.fastq.gz  

# map the reads
bismark --genome ${genome_path} \
    -o ${out_path} \
    --score_min L,0,-0.6 -R 10 \
    -p 12 \
    -1 ${trim_path}/${sampleID}*val_1.fq.gz \
    -2 ${trim_path}/${sampleID}*val_2.fq.gz

# extract methylation data
bismark_methylation_extractor -s -o ${cov_out} \
	--gzip \
	--bedGraph --comprehensive --parallel 12 \
	--merge_non_CpG --genome_folder ${genome_path} \
	${out_path}/${sampleID}*.bam

coverage2cytosine --merge_CpG --gzip --genome_folder $genome_path \
	-o ${sampleID} \
	--dir ${finalcov_out}/ \
	${cov_out}/${sampleID}*.cov.gz

#submitted with: -t 2-00:00:00 --cpus-per-task=4 --mem-per-cpu=30G -p general --array=1-TotalNumberOfLibraries trim_map_capuch_feces.sh
