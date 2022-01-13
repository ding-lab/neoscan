# neoscan pipeline v1.3 for wustl computer1 cluster#

Pipeline for detecting neoantigen from snvs and indels

## Install the third-party software ##

Install netMHC4 by following http://www.cbs.dtu.dk/services/doc/netMHC-4.0.readme

## Usage ##

Step 0: 
set environment for LSF job on compute1 by adding the following to ~/.bashrc file:

export PATH=/storage1/fs1/songcao/Active/Software/anaconda3/bin:$PATH

export STORAGE2=/storage1/fs1/dinglab/Active export SCRATCH2=/storage1/fs1/dinglab/

export STORAGE1=/storage1/fs1/songcao/Active export SCRATCH1=/storage1/fs1/songcao/

export LSF_DOCKER_VOLUMES="$STORAGE1:$STORAGE1 $STORAGE2:$STORAGE2"

step 2: Enter the directory where you downloaded neoscan pipeline, and run it step by step (total 7 steps)

perl neoscan.pl --rdir --log --bamfq --bed --step --rna --refdir --q --groupname --users

 <rdir> = full path of the folder holding files for this sequence run

 <log> = full path of the folder saving log files 

 <bam> = 1, input is bam; 0, input is fastq: default 1

 <rna> =1, input data is rna, otherwise is dna

 <bed> = bed file for annotation: refseq: /storage1/fs1/songcao/Active/Database/hg38_database/refseq/refseq_hg38_june29/proteome.bed 
 
 <refdir> = ref directory: /storage1/fs1/songcao/Active/Database/hg38_database/refseq/refseq_hg38_june29

 <groupname> = job group name

 <users> = user name for job group

 <q> = which queue for submitting job; research-hpc, ding-lab, long (default)
 
<step>  =run this pipeline step by step. (running the whole pipeline if step number is 0)

# files required in the running directory ##
 - vcf file format for snvs with columns: chromosome, start position, ref allele, alt allele, gene hugo symbol, HGSV short, is it somatic or germline mutation. Filename: <id>.snp.vcf

        1       113854971       C       G       PTPN22  p.E207Q Somatic

        1       113900168       C       A       AP4B1   p.A284S Somatic

        1       117623561       C       G       FAM46C  p.I231M Somatic
 
 - vcf file for indels with the same columns. Filename <id>.indel.vcf

        3       161086280       T       -       B3GALNT1        p.T159Pfs*8     Somatic
 
 - RNA-Seq or exome bam  or fastq file for HLA type

All three input files should be in one folder. One set of files per sample


## Contact ##

Song Cao, scao@wustl.edu 

