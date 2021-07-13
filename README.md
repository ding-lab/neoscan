# neoscan pipeline v1.3 #

Pipeline for detecting neoantigen from snvs and indels

## Install the third-party software ##

conda install -c bioconda optitype

Change the path for OptiPathPipeline.py in script neoscan.pl to where you install optitype

## Usage ##

perl neoscan.pl --rdir <rdir> --log <log> --bamfq <bamfq> --bed <bed> --rna <rna> --refdir <refdir> --step <step_number>

        <rdir> = full path of the folder holding files for this sequence run

        <log> = full path of the folder saving log files

        <bamfq> = 1, input is bam; 0, input is fastq: default 1

        <rna> =1, input data is rna, otherwise is dna. For HLA genotype

        <bed> = bed file for annotation: ensembl: /gscmnt/gc2518/dinglab/scao/db/ensembl38.85/proteome-first.bed

         refseq: /gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29/proteome.bed

        <refdir> = ref directory: /gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29

        <step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

# files required in the running directory ##
 - vcf file format for snvs with columns: chromosome, start position, ref allele, alt allele, gene hugo symbol, HGSV short, is it somatic or germline mutation. Filename: <id>.snp.vcf

        1       113854971       C       G       PTPN22  p.E207Q Somatic

        1       113900168       C       A       AP4B1   p.A284S Somatic

        1       117623561       C       G       FAM46C  p.I231M Somatic
 
 - vcf file for indels with the same columns. Filename <id>.indel.vcf

        3       161086280       T       -       B3GALNT1        p.T159Pfs*8     Somatic
 
 - RNA-Seq or exome bam  or fastq file for HLA type

All three input files should be in one folder. One set of files per sample

# Some hints for running neoscan pipleine in WU internal MGI cluster # 
    

1. Copy the tool to a folder. This command will create neoscan/ folder in your current folder:

2. git clone https://github.com/ding-lab/neoscan.git

This is required to be able to do git checkout, Needed just once:

LSF_DOCKER_PRESERVE_ENVIRONMENT=true bsub -Is -R "select[mem>15000] rusage[mem=15000]" -M 32000000 -q docker-interactive -a "docker(scao/dailybox)" /bin/bash

3. install optitype 

conda install -c bioconda optitype

Change the path for OptiPathPipeline.py in script neoscan.pl to where you install optitype

To prepare vcf and bam input files, you can follow the example at /gscmnt/gc2524/dinglab/akarpova/cptac3/CCRCC_neoscan_test.  
 
perl /gscmnt/gc2524/dinglab/akarpova/software/neoscan/neoscan.pl --rdir /gscmnt/gc2524/dinglab/akarpova/cptac3/CCRCC_neoscan_test --log /gscmnt/gc2524/dinglab/akarpova/cptac3 --bamfq 1 --bed /gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29/proteome.bed --rna 1 --refdir /gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29 --step 1

Then change --step 2/3/4/5

After finishing running step 5, you can get the final result in the followint two files:

*SAMPLE*.neo.snv.summary

*SAMPLE*.neo.indel.summary

Then you may want to filter out peptides found in human cells in general. Just grep every single peptide in this database 
/gscmnt/gc2518/dinglab/scao/db/ensembl38.85/Homo_sapiens.GRCh38.pep.all.fa.cleaned.fa

## Contact ##

Song Cao, scao@wustl.edu 

