#########Song Cao###########
# last updated date: 12/21/2016 #
### updated 1/3/2018 ###

### updated 08/11/2018 ##

### use docker image for submitting jobs to research-hpc queue ##

#!/usr/bin/perl
use strict;
use warnings;
#use POSIX;
my $version = 1.0;
#color code
my $red = "\e[31m";
my $gray = "\e[37m";
my $yellow = "\e[33m";
my $green = "\e[32m";
my $purple = "\e[35m";
my $cyan = "\e[36m";
my $normal = "\e[0m";
#usage information
(my $usage = <<OUT) =~ s/\t+//g;

This script will predict neoantigen for somatic variants in cancer sample 
Pipeline version: $version

$yellow     Usage: perl $0 <run_folder> <step_number> $normal

<run_folder> = full path of the folder holding files for this sequence run

<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)
$red     	[1] Generate fasta snv and indel
		 	[2] Generate peptide for snv and indel

$green      [3]  Run BWA 
			[4] parse bwa

$purple		[5]  Run netMHC 
		    [6] parse netMHC result
$yellow 	[7] get ref read count
			[8] Generate final report
$normal
OUT

die $usage unless @ARGV == 2;
my ( $run_dir, $step_number ) = @ARGV;
if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

die $usage unless ($step_number >=0)&&(($step_number <= 10) || ($step_number>=12));
my $email = "scao\@wustl\.edu";
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];

my $HOME1="/gscmnt/gc2524/dinglab";
#store job files here
if (! -d $HOME1."/tmpNeo") {
    `mkdir $HOME1"/tmpNeo"`;
}

my $job_files_dir = $HOME1."/tmpNeo";
#store SGE output and error files here

if (! -d $HOME1."/LSF_DIR_Neo") {
    `mkdir $HOME1"/LSF_DIR_Neo"`;
}

my $lsf_file_dir = $HOME1."/LSF_DIR_Neo";

## hlaminer for genotype, netMHC for neoantigen prediction ##

my $db_hla_abc_cds="/gscmnt/gc2523/dinglab/neoantigen/human_DB/HLA_ABC_CDS.fasta";
my $hlaminer="/gscmnt/gc2523/dinglab/neoantigen/HLAminer_v1.3.1/bin/HLAminer_song.pl";
my $hla_nom_p="/gscmnt/gc2523/dinglab/neoantigen/HLAminer_v1.3.1/database/hla_nom_p.txt";
my $f_allele="/gscmnt/gc2523/dinglab/neoantigen/netMHC-4.0/Linux_x86_64/data/allelelist";
my $netMHC="/gscmnt/gc2523/dinglab/neoantigen/netMHC-4.0/netMHC";
my $db_ref_bed="/gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29/proteome.bed";
my $h38_fa="/gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29";
my $db_cdna="/gscmnt/gc3027/dinglab/medseq/fasta/human/Homo_sapiens.GRCh37.70.cdna.all.fa";
my $db_sanger_qs="/gscmnt/gc2523/dinglab/neoantigen/neoantigen-scan/quality_sanger.table.2col.tsv";
#my $f_h37="/gscmnt/gc2522/dinglab/ref_tophat/bowtie2-index/Homo_sapiens.GRCh37.75";
#my $f_h37_gtf="/gscmnt/gc2522/dinglab/ref_tophat/bowtie2-index/Homo_sapiens.GRCh37.75.gtf";
#my $f_g="/gscmnt/gc2522/dinglab/ref_tophat/Bowtie2Index/genome";
# obtain script path
my $run_script_path = `dirname $0`;
chomp $run_script_path;
my $run_script_path_perl = "/usr/bin/perl ".$run_script_path."/";
my $run_script_path_python = "/usr/bin/python ".$run_script_path."/";
#my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

if ($step_number < 8 || $step_number>=12) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                $current_job_file="";
                if($step_number == 0 || $step_number>=12)
                { 
					if($step_number==0)
					{
						&bsub_fa();}
					if($step_number<=12)
					{
						&bsub_pep(); }
					if($step_number<=13)
					{
						&bsub_bwa();}
                    if($step_number<=14)
                    {
						&bsub_parse_bwa();}	
					if($step_number<=15)
					{
						&bsub_netmhc();}
					if($step_number<=16)
					{
						&bsub_parsemhc();}
					if($step_number<=17)
                    {
                        &bsub_get_ref_count();}
				}
				elsif ($step_number == 1) {
                    &bsub_fa();
                } 
                elsif ($step_number == 2) {
                    &bsub_pep(1);
                } elsif ($step_number == 3) {
                    &bsub_bwa(1);
				} elsif ($step_number == 4) {
                    &bsub_parse_bwa(1);
                } elsif ($step_number == 5) {
                    &bsub_netmhc(1);
                }elsif ($step_number == 6) {
                    &bsub_parsemhc(1);
                } elsif ($step_number == 7) {
                    &bsub_get_ref_count(1);
                }
				
            }
        }
    }
}

if($step_number==0 || $step_number>=12) 
	{
	print "generate the final report\n"; 
	&bsub_final_report(); 
	}

if($step_number==8) 
	{
	print "generate the final report\n";
    &bsub_final_report(1);
    }
#######################################################################
# send email to notify the finish of the analysis

if (($step_number == 0) || ($step_number == 9)) {
    print $yellow, "Submitting the job for sending an email when the run finishes ",$sample_name, "...",$normal, "\n";
    $hold_job_file = $current_job_file;
    $current_job_file = "Email_run_".$$.".sh";
    open(EMAIL, ">$job_files_dir/$current_job_file") or die $!;
    print EMAIL "#!/bin/bash\n";
    print EMAIL "#BSUB -n 1\n";
    print EMAIL "#BSUB -o $lsf_file_dir","\n";
    print EMAIL "#BSUB -e $lsf_file_dir","\n";
    print EMAIL "#BSUB -J $current_job_file\n";
    print EMAIL "#BSUB -w \"$hold_job_file\"","\n";
	print EMAIL "#BSUB -q ding-lab\n";
    print EMAIL $run_script_path."send_email.pl ".$run_dir." ".$email."\n";
    close EMAIL;
    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #$bsub_com = "qsub -V -hold_jid $hold_job_file -e $lsf_file_dir -o $lsf_file_dir $job_files_dir/$current_job_file\n";
    system ($bsub_com);
}
#######################################################################
if ($step_number == 0) {
    print $green, "All jobs are submitted! You will get email notification when this run is completed.\n",$normal;
}

### step 1: get fasta sequence for transcript with indel and snv mutation ##

sub bsub_fa{

    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j1_fa_".$sample_name.".sh";
    #my $IN_sam = $sample_full_path."/".$sample_name.".exome.sam"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
 
   open(FA, ">$job_files_dir/$current_job_file") or die $!;
    print FA "#!/bin/bash\n";
    #print FA "#BSUB -n 1\n";
    #print FA "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print FA "#BSUB -M 30000000\n";
   # print FA "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
   # print FA "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
   # print FA "#BSUB -q ding-lab\n";
	print FA "#BSUB -J $current_job_file\n";
   # print FA "#BSUB -w \"$hold_job_file\"","\n";
    print FA "f_vcf_indel=".$sample_full_path."/$sample_name.indel.vcf\n";
    print FA "f_vcf_snp=".$sample_full_path."/$sample_name.snp.vcf\n";
    print FA "f_fa_indel_wt=".$sample_full_path."/$sample_name.indel.vcf.fa-indel-wt.fasta\n";
    print FA "f_fa_indel_mut=".$sample_full_path."/$sample_name.indel.vcf.fa-indel-mut.fasta\n";
    print FA "f_fa_snp_wt=".$sample_full_path."/$sample_name.snp.vcf.fa-snv-wt.fasta\n";
    print FA "f_fa_snp_mut=".$sample_full_path."/$sample_name.snp.vcf.fa-snv-mut.fasta\n";
    print FA "f_fa_all=".$sample_full_path."/$sample_name.transcript.fa\n";
    print FA 'if [ -s $f_vcf_snp ]',"\n";
    print FA "then\n";
    print FA " ".$run_script_path_perl."fasta_seq_for_snv_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_snp}"."\n";
    print FA "cat \${f_fa_snp_mut} > \${f_fa_all}","\n";
    print FA "fi\n\n";
    print FA 'if [ -s $f_vcf_indel ]',"\n";
    print FA "then\n";
    print FA " ".$run_script_path_perl."fasta_seq_for_indel_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_indel}"."\n";
    print FA "cat \${f_fa_indel_mut} >> \${f_fa_all}","\n";
    print FA "fi\n\n";
    print FA "";
    close FA;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 

    system ( $bsub_com );
}

##step 2: get peptide sequences for transcript with indel and snv mutations ##

sub bsub_pep{
	
	my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }
	else{
        $hold_job_file = $current_job_file;
    }
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
  
	$current_job_file = "j2_pep_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;


    #my $IN_sam = $sample_full_path."/".$sample_name.".exome.sam"; 
    open(PEP, ">$job_files_dir/$current_job_file") or die $!;
    print PEP "#!/bin/bash\n";
#    print PEP "#BSUB -n 1\n";
#    print PEP "#BSUB -R \"rusage[mem=30000]\"","\n";
#    print PEP "#BSUB -M 30000000\n";
#	print PEP "#BSUB -q ding-lab\n";
#    print PEP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
#    print PEP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
#    print PEP "#BSUB -J $current_job_file\n";
#    print PEP "#BSUB -w \"$hold_job_file\"","\n";
    print PEP "f_vcf_indel=".$sample_full_path."/$sample_name.indel.vcf\n";
    print PEP "f_vcf_snp=".$sample_full_path."/$sample_name.snp.vcf\n";
    print PEP "f_pep_indel_wt=".$sample_full_path."/$sample_name.indel.vcf.proteome-indel-wt.fasta\n";
    print PEP "f_pep_indel_mut=".$sample_full_path."/$sample_name.indel.vcf.proteome-indel-mut.fasta\n";
    print PEP "f_pep_snp_wt=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-wt.fasta\n";
    print PEP "f_pep_snp_mut=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-mut.fasta\n";
    print PEP "f_pep_all=".$sample_full_path."/$sample_name.pep.fa\n";
    print PEP 'if [ -s $f_vcf_snp ]',"\n";
    print PEP "then\n";
    print PEP " ".$run_script_path_perl."protein_seq_for_snv_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_snp}"."\n";
    print PEP "cat \${f_pep_snp_mut} > \${f_pep_all}","\n";
    print PEP "fi\n\n";
    print PEP 'if [ -s $f_vcf_indel ]',"\n";
    print PEP "then\n";
    print PEP " ".$run_script_path_perl."protein_seq_for_indel_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_indel}"."\n";
    print PEP "cat \${f_pep_indel_mut} >> \${f_pep_all}","\n";
    print PEP "fi\n\n";
    print PEP "";
    close PEP;
    my $sh_file=$job_files_dir."/".$current_job_file;

   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}


##step 3: use bwa to align the hla reference to get HLA genotypet ## 

sub bsub_bwa{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j3_bwa_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    my $IN_bam = $sample_full_path."/".$sample_name.".rnaseq.bam";
    my $f_fq_1=$sample_full_path."/".$sample_name.".1.fq";
    my $f_fq_2=$sample_full_path."/".$sample_name.".2.fq";
#	my $RNASEQ_MT_sam=$sample_full_path."/".$sample_name.".rnaseq.mapped.sam";

    if ((! -e $IN_bam) && ((!-e $f_fq_1) || (! -e $f_fq_2))) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: Died because there is no input bam or fq file for bwa:\n";
        print "File $IN_bam does not exist!\n";
        die "Please check command line argument!", $normal, "\n\n";

    }
    #if (! -s $IN_bam) {#make sure input fasta file is not empty
    #    print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
     #   die "Warning: Died because $IN_bam is empty!", $normal, "\n\n";
    #}
    open(RNASEQ, ">$job_files_dir/$current_job_file") or die $!;
  # my $f_fq_1=$sample_full_path."/".$sample_name.".1.fq";
   # my $f_fq_2=$sample_full_path."/".$sample_name.".2.fq";

    print RNASEQ "#!/bin/bash\n";
    #print RNASEQ "#BSUB -n 1\n";
    #print RNASEQ "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print RNASEQ "#BSUB -M 30000000\n";
    #print RNASEQ "#BSUB -q ding-lab\n";
	#print RNASEQ "#BSUB -w \"$hold_job_file\"","\n";
	#print RNASEQ "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print RNASEQ "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print RNASEQ "#BSUB -J $current_job_file\n";
    print RNASEQ "RNASEQ_IN=".$sample_full_path."/".$sample_name.".rnaseq.bam\n";
    print RNASEQ "RNASEQ_sorted=".$sample_full_path."/".$sample_name.".rnaseq.sorted\n";
    print RNASEQ "RNASEQ_sorted_bam=".$sample_full_path."/".$sample_name.".rnaseq.sorted.bam\n";
	print RNASEQ "UNMAPPED_IN=".$sample_full_path."/".$sample_name.".rnaseq.unaligned_reads.bam\n";
    print RNASEQ "UNMAPPED_sorted=".$sample_full_path."/".$sample_name.".rnaseq.unaligned_reads.sorted\n";
    print RNASEQ "UNMAPPED_sorted_bam=".$sample_full_path."/".$sample_name.".rnaseq.unaligned_reads.sorted.bam\n";
	print RNASEQ "RNASEQ_fq_1=".$sample_full_path."/".$sample_name.".1.fq\n";
    print RNASEQ "RNASEQ_fq_2=".$sample_full_path."/".$sample_name.".2.fq\n";
	print RNASEQ "RNASEQ_sai_1=".$sample_full_path."/".$sample_name.".1.sai\n";
	print RNASEQ "RNASEQ_sai_2=".$sample_full_path."/".$sample_name.".2.sai\n";
	print RNASEQ "RNASEQ_sam=".$sample_full_path."/".$sample_name.".rnaseq.sam\n";	
	print RNASEQ "INDEX_T=".$sample_full_path."/".$sample_name.".transcript.fa\n";
	print RNASEQ "RNASEQ_sai_T1=".$sample_full_path."/".$sample_name.".T1.sai\n";
    print RNASEQ "RNASEQ_sai_T2=".$sample_full_path."/".$sample_name.".T2.sai\n";
    print RNASEQ "RNASEQ_T_sam=".$sample_full_path."/".$sample_name.".T.rnaseq.sam\n";
	print RNASEQ "HLA_log=".$sample_full_path."/HLAminer_HPRA.log\n";
    print RNASEQ "HLA_csv=".$sample_full_path."/HLAminer_HPRA.csv\n";
    print RNASEQ "HLA_tsv=".$sample_full_path."/HLAminer_alleles.tsv\n";
    #print HLA "HLA_bam=".$sample_full_path."/$sample_name.hla.bam\n";
	##mapped reads##
	print RNASEQ 'if [ -f $IN_bam ]',"\n"; # input file exist
	print RNASEQ "then\n";
    print RNASEQ "samtools sort -n \${RNASEQ_IN} \${RNASEQ_sorted}","\n";
	print RNASEQ "samtools view \${RNASEQ_sorted_bam} | perl -ne \'\$l=\$_; \$f_q1=\"$f_fq_1\"; \$f_q2=\"$f_fq_2\"; if(\$first==0) { open(OUT1,\">\$f_q1\"); open(OUT2,\">\$f_q2\");  \$first=1;}  \@ss=split(\"\\t\",\$l); \$flag=\$ss[1]; \$cigar=\$ss[5]; if((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) { next; } \$id=\$ss[0]; \$seq=\$ss[9]; \$q=\$ss[10];  if(\$id=~/\\/1\$/ || (\$flag & 0x40) ) { \$r1=\$id; \$r1=~s/\\/1\$//g; \$seq1=\$seq; \$q1=\$q; } if(\$id=~/\\/2\$/ || (\$flag & 0x80)) { \$r2=\$id; \$r2=~s/\\/2\$//g; \$seq2=\$seq; \$q2=\$q; } if((\$r1 eq \$r2)) { print OUT1 \"\@\",\$r1,\"/1\",\"\\n\"; print OUT1 \$seq1,\"\\n\"; print OUT1 \"+\",\"\\n\"; print OUT1 \$q1,\"\\n\"; print OUT2 \"\@\",\$r1,\"/2\",\"\\n\"; print OUT2 \$seq2,\"\\n\"; print OUT2 \"+\",\"\\n\"; print OUT2 \$q2,\"\\n\";}\'","\n";
	print RNASEQ "  fi\n\n";
	### unmapped reads ###
	print RNASEQ 'if [ -f $UNMAPPED_IN ]',"\n"; # input file exist
    print RNASEQ "then\n";
    print RNASEQ "samtools sort -n \${UNMAPPED_IN} \${UNMAPPED_sorted}","\n";    
    print RNASEQ "samtools view \${UNMAPPED_sorted_bam} | perl -ne \'\$l=\$_; \$f_q1=\"$f_fq_1\"; \$f_q2=\"$f_fq_2\"; if(\$first2==0) { open(OUT1,\">>\$f_q1\"); open(OUT2,\">>\$f_q2\");  \$first2=1;}  \@ss=split(\"\\t\",\$l); \$flag=\$ss[1]; \$cigar=\$ss[5]; if((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) { next; } \$id=\$ss[0]; \$seq=\$ss[9]; \$q=\$ss[10];  if(\$id=~/\\/1\$/ || (\$flag & 0x40) ) { \$r1=\$id; \$r1=~s/\\/1\$//g; \$seq1=\$seq; \$q1=\$q; } if(\$id=~/\\/2\$/ || (\$flag & 0x80)) { \$r2=\$id; \$r2=~s/\\/2\$//g; \$seq2=\$seq; \$q2=\$q; } if((\$r1 eq \$r2)) { print OUT1 \"\@\",\$r1,\"/1\",\"\\n\"; print OUT1 \$seq1,\"\\n\"; print OUT1 \"+\",\"\\n\"; print OUT1 \$q1,\"\\n\"; print OUT2 \"\@\",\$r1,\"/2\",\"\\n\"; print OUT2 \$seq2,\"\\n\"; print OUT2 \"+\",\"\\n\"; print OUT2 \$q2,\"\\n\";}\'","\n";
    #print RNASEQ "bwa0.7.12 mem -t 8 $db_hla_abc_cds \${RNASEQ_fq_1} \${RNASEQ_fq_2}  > \${RNASEQ_sam}"."\n"; 
    print RNASEQ "  fi\n\n";
	### bwa alignement ###
	#print RNASEQ "bwa0.7.12 index \${INDEX_T}"."\n";
    print RNASEQ "rm \${RNASEQ_sorted_bam}","\n";
    print RNASEQ 'if [ -f $UNMAPPED_sorted_bam ]',"\n"; # input file exist
    print RNASEQ "then\n";
    print RNASEQ "rm \${UNMAPPED_sorted_bam}","\n";
    print RNASEQ "fi\n";
	print RNASEQ "bwa0.7.12 aln -e 0 -o 0 -t 4 -f \${RNASEQ_sai_1} $db_hla_abc_cds \${RNASEQ_fq_1}"."\n";
	print RNASEQ "bwa0.7.12 aln -e 0 -o 0 -t 4 -f \${RNASEQ_sai_2} $db_hla_abc_cds \${RNASEQ_fq_2}"."\n";
	print RNASEQ "bwa0.7.12 sampe -f \${RNASEQ_sam} $db_hla_abc_cds \${RNASEQ_sai_1} \${RNASEQ_sai_2} \${RNASEQ_fq_1} \${RNASEQ_fq_2}"."\n";
    print RNASEQ "perl $hlaminer -h $db_hla_abc_cds -p $hla_nom_p -a \${RNASEQ_sam} -o $sample_full_path"."\n";
    print RNASEQ  " ".$run_script_path_perl."parseHLAresult.pl \${HLA_csv} \${HLA_tsv}"."\n";
	print RNASEQ "rm \${RNASEQ_sam}","\n";
	print RNASEQ "bwa0.7.12 index \${INDEX_T}"."\n";
    print RNASEQ "bwa0.7.12 aln -e 0 -o 0 -t 4 -f \${RNASEQ_sai_T1} \${INDEX_T} \${RNASEQ_fq_1}"."\n";
    print RNASEQ "bwa0.7.12 aln -e 0 -o 0 -t 4 -f \${RNASEQ_sai_T2} \${INDEX_T} \${RNASEQ_fq_2}"."\n";
    print RNASEQ "bwa0.7.12 sampe -f \${RNASEQ_T_sam} \${INDEX_T} \${RNASEQ_sai_T1} \${RNASEQ_sai_T2} \${RNASEQ_fq_1} \${RNASEQ_fq_2}"."\n";
	print RNASEQ "rm \${RNASEQ_fq_1}","\n";
 	print RNASEQ "rm \${RNASEQ_fq_2}","\n";
    print RNASEQ "rm \${RNASEQ_sai_1}","\n";
	print RNASEQ "rm \${RNASEQ_sai_2}","\n";
    print RNASEQ "rm \${RNASEQ_sai_T1}","\n";
    print RNASEQ "rm \${RNASEQ_sai_T2}","\n";
    close RNASEQ;
    my $sh_file=$job_files_dir."/".$current_job_file;

   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );
}

sub bsub_parse_bwa{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

	#my $RNASEQ_MT_sam=$sample_full_path."/".$sample_name.".rnaseq.mapped.sam";

    $current_job_file = "j4_parse_bwa_".$sample_name.".sh";

    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    my $RNASEQ_MT_sam=$sample_full_path."/".$sample_name.".rnaseq.T.mapped.sam";

    open(BWAP, ">$job_files_dir/$current_job_file") or die $!;

    print BWAP "#!/bin/bash\n";
    #print BWAP "#BSUB -n 1\n";
    #print BWAP "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print BWAP "#BSUB -q ding-lab\n";
  	#print BWAP "#BSUB -M 30000000\n";
    #print BWAP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print BWAP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
    #print BWAP "#BSUB -J $current_job_file\n";
    #print BWAP "#BSUB -w \"$hold_job_file\"","\n";
	print BWAP "RNASEQ_T_sam=".$sample_full_path."/".$sample_name.".T.rnaseq.sam\n";	
    print BWAP "CDNA_sai=".$sample_full_path."/".$sample_name.".cdna.sai\n";
    print BWAP "CDNA_sam=".$sample_full_path."/".$sample_name.".cdna.sam\n";
    print BWAP "RNASEQ_T_fq=".$sample_full_path."/".$sample_name.".T.rnaseq.fq\n";
    print BWAP "RNASEQ_TF_fq=".$sample_full_path."/".$sample_name.".T.rnaseq.filtered.fq\n";
    print BWAP "f_fa_indel_wt=".$sample_full_path."/$sample_name.indel.vcf.fa-indel-wt.fasta\n";
    print BWAP "f_fa_indel_mut=".$sample_full_path."/$sample_name.indel.vcf.fa-indel-mut.fasta\n";
    print BWAP "f_fa_snp_wt=".$sample_full_path."/$sample_name.snp.vcf.fa-snv-wt.fasta\n";
    print BWAP "f_fa_snp_mut=".$sample_full_path."/$sample_name.snp.vcf.fa-snv-mut.fasta\n";
	print BWAP "f_novel_seq=".$sample_full_path."/$sample_name.novelseq.fa\n";
	print BWAP "f_novel_q20_seq=".$sample_full_path."/$sample_name.novelseq.q20.fa\n";
	print BWAP "f_indel_pos=".$sample_full_path."/$sample_name.indel.fa.pos.tsv\n";
	print BWAP "f_snv_pos=".$sample_full_path."/$sample_name.snv.fa.pos.tsv\n";
	print BWAP "f_mut_pos=".$sample_full_path."/$sample_name.mut.fa.pos.tsv\n";
    print BWAP "f_fa_all=".$sample_full_path."/$sample_name.transcript.fa\n";
	print BWAP "f_sup=".$sample_full_path."/$sample_name.transcript.novel.support.rna.tsv\n";
	print BWAP  " ".$run_script_path_perl."get_perferct_mapped_reads.pl \${RNASEQ_T_sam} \${RNASEQ_T_fq}"."\n";	
	print BWAP "bwa0.7.12 aln $db_cdna \${RNASEQ_T_fq} > \${CDNA_sai}"."\n";
 	print BWAP "bwa0.7.12 samse $db_cdna \${CDNA_sai} \${RNASEQ_T_fq} > \${CDNA_sam}"."\n";
	print BWAP " ".$run_script_path_perl."filter_reads_in_cdna.pl \${CDNA_sam} \${RNASEQ_T_fq} \${RNASEQ_TF_fq}"."\n";
	print BWAP  " ".$run_script_path_perl."extract_novel_sequence.pl \${f_fa_all} \${f_fa_snp_wt} \${f_fa_indel_wt} \${RNASEQ_TF_fq} \${f_novel_seq}"."\n";
	print BWAP " ".$run_script_path_perl."q20_filter.pl $db_sanger_qs \${f_novel_seq} \${RNASEQ_TF_fq} \${f_novel_q20_seq}"."\n";
	print BWAP 'if [ -s $f_fa_indel_wt ]',"\n";
    print BWAP "then\n";
	print BWAP " ".$run_script_path_perl."find_first_diff_pos_mut_wt_indel.pl \${f_fa_indel_wt} \${f_fa_indel_mut} \${f_indel_pos}"."\n";	
	print BWAP "fi\n\n";
	print BWAP " ".$run_script_path_perl."find_first_diff_pos_mut_wt_snv.pl \${f_fa_snp_wt} \${f_fa_snp_mut} \${f_snv_pos}"."\n";
	print BWAP "cat \${f_indel_pos} \${f_snv_pos} > \${f_mut_pos}","\n"; 
	print BWAP " ".$run_script_path_perl."extract_support_reads_from_novel_sequence.pl \${f_fa_all} \${f_novel_q20_seq} \${f_mut_pos} \${f_sup}"."\n";
	print BWAP "rm \${RNASEQ_T_sam}","\n";
	print BWAP "rm \${CDNA_sam}","\n";
	close BWAP;
    my $sh_file=$job_files_dir."/".$current_job_file;

   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );

    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );	
}

sub bsub_netmhc{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j5_bind_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    #my $IN_sam = $sample_full_path."/".$sample_name.".exome.sam"; 
	open(MHC, ">$job_files_dir/$current_job_file") or die $!;
    print MHC "#!/bin/bash\n";
    #print MHC "#BSUB -n 1\n";
    #print MHC "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print MHC "#BSUB -M 30000000\n";
	#print MHC "#BSUB -q ding-lab\n";
   # print MHC "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
   # print MHC "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
   # print MHC "#BSUB -J $current_job_file\n";
   # print MHC "#BSUB -w \"$hold_job_file\"","\n";
	print MHC "f_pep=".$sample_full_path."/$sample_name.pep.fa\n";
    print MHC "HLA_tsv=".$sample_full_path."/HLAminer_alleles.tsv\n";
	print MHC "f_netMHC_result=".$sample_full_path."/netMHC4.0.out.append.txt\n";
	print MHC "f_out=".$sample_full_path."/result_neoantigen\n";
    print MHC  " ".$run_script_path_python."runNetMHC4.py -a \${HLA_tsv} -f \${f_pep} -p 8,9,10,11 -o $sample_full_path -n $netMHC -v $f_allele"."\n";
    close MHC;

    my $sh_file=$job_files_dir."/".$current_job_file;

  #  $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -o $lsf_out -e $lsf_err sh $sh_file\n";

    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n";

    system ( $bsub_com );
	
}

sub bsub_parsemhc{
    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }   
 #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j6_parsebind_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    #my $IN_sam = $sample_full_path."/".$sample_name.".exome.sam"; 
    open(PMHC, ">$job_files_dir/$current_job_file") or die $!;
    print PMHC "#!/bin/bash\n";
   # print PMHC "#BSUB -n 1\n";
   # print PMHC "#BSUB -R \"rusage[mem=30000]\"","\n";
   # print PMHC "#BSUB -M 30000000\n";
#	print PMHC "#BSUB -q ding-lab\n";
#    print PMHC "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
 #   print PMHC "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
  #  print PMHC "#BSUB -J $current_job_file\n";
#	print PMHC "#BSUB -w \"$hold_job_file\"","\n";
    print PMHC "f_indel=".$sample_full_path."/$sample_name.indel.vcf\n";
    print PMHC "f_snv=".$sample_full_path."/$sample_name.snp.vcf\n";	
    print PMHC "f_indel_wt_fa=".$sample_full_path."/$sample_name.indel.vcf.proteome-indel-wt.fasta\n";
	print PMHC "f_snv_wt_fa=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-wt.fasta\n";
    print PMHC "f_indel_mut_fa=".$sample_full_path."/$sample_name.indel.vcf.proteome-indel-mut.fasta\n";
    print PMHC "f_snv_mut_fa=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-mut.fasta\n";
    print PMHC "f_netMHC_result=".$sample_full_path."/netMHC4.0.out.append.txt\n";
    print PMHC "f_out=".$sample_full_path."/$sample_name.neoantigen.tsv\n";
	print PMHC "f_sum=".$sample_full_path."/$sample_name.neo.summary\n";
	print PMHC "f_min=".$sample_full_path."/$sample_name.neo.summary.min\n";
    print PMHC  " ".$run_script_path_perl."parseNetMHC4result.pl \${f_netMHC_result} \${f_indel_wt_fa} \${f_snv_wt_fa} \${f_out}"."\n";
    print PMHC	" ".$run_script_path_perl."reportSummary.pl \${f_out} \${f_indel} \${f_snv} \${f_indel_mut_fa} \${f_snv_mut_fa} \${f_sum}"."\n";
	print PMHC  " ".$run_script_path_perl."get_min_result.pl \${f_sum} \${f_min}"."\n";
	print PMHC	"rm \${f_netMHC_result}"."\n";
	close PMHC;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );
    my $sh_file=$job_files_dir."/".$current_job_file;
   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}

sub bsub_get_ref_count()
 {

	my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }
 #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j7_getrefcount_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    #my $IN_sam = $sample_full_path."/".$sample_name.".exome.sam"; 
	my $IN_bam = $sample_full_path."/".$sample_name.".rnaseq.bam";
    open(REFC, ">$job_files_dir/$current_job_file") or die $!;
    print REFC "#!/bin/bash\n";
  #  print REFC "#BSUB -n 1\n";
  #  print REFC "#BSUB -R \"rusage[mem=30000]\"","\n";
  #  print REFC "#BSUB -M 30000000\n";
  #  print REFC "#BSUB -q ding-lab\n";
   # print REFC "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
   # print REFC "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
   # print REFC "#BSUB -J $current_job_file\n";
   # print REFC "#BSUB -w \"$hold_job_file\"","\n";
	print REFC "f_fa_all=".$sample_full_path."/$sample_name.transcript.fa\n";
	print REFC " ".$run_script_path_perl."extract_reads_supporting_ref.pl $IN_bam  $db_ref_bed \${f_fa_all}"."\n";
	close REFC;
#    $bsub_com = "bsub < $job_files_dir/$current_job_file\n";
#    system ( $bsub_com );
    my $sh_file=$job_files_dir."/".$current_job_file;
   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
	}

sub bsub_final_report()
	{   

	my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

 	#my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    $current_job_file = "j8_final_report_".$working_name.".sh";
    #my $IN_sam = $sample_full_path."/".$sample_name.".exome.sam"; 
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;

    open(REP, ">$job_files_dir/$current_job_file") or die $!;
    print REP "#!/bin/bash\n";
    #print REP "#BSUB -n 1\n";
    #print REP "#BSUB -R \"rusage[mem=30000]\"","\n";
    #print REP "#BSUB -M 30000000\n";
    #print REP "#BSUB -o $lsf_file_dir","/","$current_job_file.out\n";
    #print REP "#BSUB -e $lsf_file_dir","/","$current_job_file.err\n";
	#print REP "#BSUB -q ding-lab\n";
    #print REP "#BSUB -J $current_job_file\n";
    #print REP "#BSUB -w \"$hold_job_file\"","\n";
	print REP " ".$run_script_path_perl."generate_report_summary.pl $run_dir $version"."\n";
	close REP;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );

    my $sh_file=$job_files_dir."/".$current_job_file;
   # $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -q research-hpc -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
    
}




