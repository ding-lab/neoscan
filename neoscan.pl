#########Song Cao###########
# last updated date: 12/21/2016 #
### updated 1/3/2018 ###

### updated 08/11/2018 ##

### use docker image for submitting jobs to dinglab queue ##

#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#use POSIX;
my $version = 1.2;
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

$yellow     Usage: perl $0 --rdir --log --bamfq --bed --step --rna --refdir --q --groupname --users $normal

<rdir> = full path of the folder holding files for this sequence run

<log> = full path of the folder saving log files 

<bam> = 1, input is bam; 0, input is fastq: default 1

<rna> =1, input data is rna, otherwise is dna

<bed> = bed file for annotation: refseq: /storage1/fs1/songcao/Active/Database/hg38_database/refseq/refseq_hg38_june29/proteome.bed
 
<refdir> = ref directory:  /storage1/fs1/songcao/Active/Database/hg38_database/refseq/refseq_hg38_june29

<groupname> = job group name

<users> = user name for job group

<q> which queue for submitting job; research-hpc, ding-lab, long (default)
 
<step_number> run this pipeline step by step. (running the whole pipeline if step number is 0)

$red     	[1] Generate fasta snv and indel
		 	[2] Generate peptide for snv and indel
$green      [3]  Run HLA type
$purple		[4]  Run netMHC 
		    [5] parse netMHC result
$cyan 		[6] generate report summary
$normal
OUT

my $step_number = -1;
my $help = 0;
my $run_dir="";
my $log_dir="";
my $hla = 1; 
my $s_bam_fq =1; 
my $s_rna=1; 
my $compute_username="";
my $group_name="";
my $q_name="";
#my $db_ref_bed="/gscmnt/gc2518/dinglab/scao/db/ensembl38.85/proteome-first.bed";
#my $db_ref_bed="/gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29/proteome.bed";
#my $h38_fa="/gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29";

my $db_ref_bed;
my $h38_fa; 
my $status = &GetOptions (
      "step=i" => \$step_number,
	  "rdir=s" => \$run_dir,
	  "bed=s" => \$db_ref_bed,
	  "refdir=s" => \$h38_fa,
              "groupname=s" => \$group_name,
      "users=s" => \$compute_username,
    "q=s" => \$q_name,
 	  "bam=i" => \$s_bam_fq,
	  "rna=i" => \$s_rna,		
      "log=s"  => \$log_dir,
	  "help" => \$help
    );

if ($help || $run_dir eq "" || $log_dir eq "" || $group_name eq "" || $compute_username eq "" || $step_number<0) {
      print $usage;
      exit;
   }
		
#die $usage unless @ARGV == 2;
#my ( $run_dir, $step_number ) = @ARGV;

if ($run_dir =~/(.+)\/$/) {
    $run_dir = $1;
}

#die $usage unless ($step_number >=0)&&(($step_number <= 10) || ($step_number>=12));
#my $email = "scao\@wustl\.edu";
my $HOME = $ENV{HOME};
my $working_name= (split(/\//,$run_dir))[-1];

my $HOME1=$log_dir;;

if (! -d $HOME1)
{
`mkdir $HOME1`;
}

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

#my $db_hla_abc_cds="/gscmnt/gc2523/dinglab/neoantigen/human_DB/HLA_ABC_CDS.fasta";
my $optitype="OptiTypePipeline.py"; 
my $f_allele="/storage1/fs1/songcao/Active/Software/netMHC-4.0/Linux_x86_64/data/allelelist";
my $netMHC="/storage1/fs1/songcao/Active/Software/netMHC-4.0/netMHC";
my $samtools="/storage1/fs1/songcao/Active/Software/anaconda3/bin/samtools";
#my $db_ref_bed="/gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29/proteome.bed";
#my $h38_fa="/gscmnt/gc2518/dinglab/scao/db/refseq_hg38_june29";
#my $f_opti_config = "/gscmnt/gc2518/dinglab/scao/home/git/neoscan/config.ini";
 
#my $db_cdna="/gscmnt/gc3027/dinglab/medseq/fasta/human/Homo_sapiens.GRCh37.70.cdna.all.fa";
#my $db_sanger_qs="/gscmnt/gc2523/dinglab/neoantigen/neoantigen-scan/quality_sanger.table.2col.tsv";

my $run_script_path = `dirname $0`;
chomp $run_script_path;
my $run_script_path_perl = "/usr/bin/perl ".$run_script_path."/";
#my $run_script_path_python = "/storage1/fs1/songcao/Active/Software/anaconda3/bin/python ".$run_script_path."/";
my $run_script_path_python = "/usr/bin/python ".$run_script_path."/";

#/usr/bin/python
#my $hold_RM_job = "norm";
my $current_job_file = "";#cannot be empty
my $hold_job_file = "";
my $bsub_com = "";
my $sample_full_path = "";
my $sample_name = "";
opendir(DH, $run_dir) or die "Cannot open dir $run_dir: $!\n";
my @sample_dir_list = readdir DH;
close DH;

if ($step_number < 6) {
    #begin to process each sample
    for (my $i=0;$i<@sample_dir_list;$i++) {#use the for loop instead. the foreach loop has some problem to pass the global variable $sample_name to the sub functions
        $sample_name = $sample_dir_list[$i];
        if (!($sample_name =~ /\./ || $sample_name=~/worklog/)) {
            $sample_full_path = $run_dir."/".$sample_name;
            if (-d $sample_full_path) { # is a full path directory containing a sample
                print $yellow, "\nSubmitting jobs for the sample ",$sample_name, "...",$normal, "\n";
                $current_job_file="";
				if ($step_number == 1) {
                    &bsub_fa();
                } 
                elsif ($step_number == 2) {
                    &bsub_pep(1);
                } elsif ($step_number == 3) {
                    &bsub_hla(1);
				} elsif ($step_number == 4) {
                    &bsub_netmhc(1);
                }elsif ($step_number == 5) {
                    &bsub_parsemhc(1);
                } 				
            }
        }
    }
}


if($step_number==6) 
	{
	print "generate the final report\n"; 
	&bsub_final_report(); 
	}


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
    #print FA "#BSUB -J $current_job_file\n";
   # print FA "#BSUB -w \"$hold_job_file\"","\n";
    print FA "f_vcf_indel=".$sample_full_path."/$sample_name.indel.vcf\n";
    print FA "f_vcf_snp=".$sample_full_path."/$sample_name.snp.vcf\n";
    print FA "f_fa_indel_wt=".$sample_full_path."/$sample_name.indel.vcf.fa-indel-wt.fasta\n";
    print FA "f_fa_indel_mut=".$sample_full_path."/$sample_name.indel.vcf.fa-indel-mut.fasta\n";
    print FA "f_fa_snv_wt=".$sample_full_path."/$sample_name.snp.vcf.fa-snv-wt.fasta\n";
    print FA "f_fa_snv_mut=".$sample_full_path."/$sample_name.snp.vcf.fa-snv-mut.fasta\n";
    print FA "f_fa_all=".$sample_full_path."/$sample_name.transcript.fa\n";
    print FA 'if [ -s $f_vcf_snp ]',"\n";
    print FA "then\n";
    print FA " ".$run_script_path_perl."fasta_seq_for_snv_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_snp}"."\n";
    print FA "cat \${f_fa_snv_mut} > \${f_fa_all}","\n";
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

    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>20000] rusage[mem=20000]\" -M 20000000 -a \'docker(ubuntu)\' -o $lsf_out -e $lsf_err \'sh $sh_file\'\n";

    #$bsub_com = "bsub -q dinglab -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -o $lsf_out -e $lsf_err sh $sh_file\n"; 

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
    print PEP "f_pep_snv_wt=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-wt.fasta\n";
    print PEP "f_pep_snv_mut=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-mut.fasta\n";
    #print PEP "f_pep_all=".$sample_full_path."/$sample_name.pep.fa\n";
    print PEP 'if [ -s $f_vcf_snp ]',"\n";
    print PEP "then\n";
    print PEP " ".$run_script_path_perl."protein_seq_for_snv_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_snp}"."\n";
    #print PEP "cat \${f_pep_snv_mut} > \${f_pep_all}","\n";
    print PEP "fi\n\n";
    print PEP 'if [ -s $f_vcf_indel ]',"\n";
    print PEP "then\n";
    print PEP " ".$run_script_path_perl."protein_seq_for_indel_using_refseq_bed.pl $db_ref_bed $h38_fa \${f_vcf_indel}"."\n";
   # print PEP "cat \${f_pep_indel_mut} >> \${f_pep_all}","\n";
    print PEP "fi\n\n";
    print PEP "";
    close PEP;
    my $sh_file=$job_files_dir."/".$current_job_file;

   # $bsub_com = "bsub -q dinglab -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>20000] rusage[mem=20000]\" -M 20000000 -a \'docker(ubuntu)\' -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
}


##step 3: use bwa to align the hla reference to get HLA genotypet ## 

sub bsub_hla{
    #my $cdhitReport = $sample_full_path."/".$sample_name.".fa.cdhitReport";
    my ($step_by_step) = @_;

    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }

    $current_job_file = "j3_hla_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    my $IN_bam = $sample_full_path."/".$sample_name.".bam";
    my $f_fq_1=$sample_full_path."/".$sample_name.".1.fq";
    my $f_fq_2=$sample_full_path."/".$sample_name.".2.fq";
#	my $HLA_MT_sam=$sample_full_path."/".$sample_name.".rnaseq.mapped.sam";
	my $dir_hla=$sample_full_path."/hla";

	#if(-d $dir_hla) { `rm $dir_hla`; }

	my $f_hla_type=`find $dir_hla -name \'*.tsv\'`;
	chomp($f_hla_type);
	print $f_hla_type,"\n"; 
	
#	if(-e $f_hla_type) { print "existing \n"; }

    if ((! -e $IN_bam) && ((!-e $f_fq_1) || (! -e $f_fq_2)) && (!-e $f_hla_type)) {#make sure there is a input fasta file 
        print $red,  "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
        print "Warning: there is no input bam or fq file for bwa:\n";
        print "File $IN_bam does not exist!\n";
		#last;
        #die "Please check command line argument!", $normal, "\n\n";
    }

    #if (! -s $IN_bam) {#make sure input fasta file is not empty
    #    print $red, "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&\n";
     #   die "Warning: Died because $IN_bam is empty!", $normal, "\n\n";
    #}
    open(HLA, ">$job_files_dir/$current_job_file") or die $!;
    print HLA "#!/bin/bash\n";
    print HLA "HLA_IN=".$sample_full_path."/".$sample_name.".bam\n";
    #print HLA "HLA_sorted=".$sample_full_path."/".$sample_name.".sorted\n";
    print HLA "HLA_sorted_bam=".$sample_full_path."/".$sample_name.".sorted.bam\n";
    print HLA "HLA_tsv=".$sample_full_path."/HLA_alleles.tsv\n";
    print HLA "f_optitype_hla=".$sample_full_path."/hla/notexisting.tsv\n";
    print HLA "if [ $s_bam_fq -eq 1 ]\n";
    print HLA "then\n";
    print HLA 'if [ -f $IN_bam ]',"\n"; # input file exist
    print HLA "then\n";
    print HLA "$samtools sort -n \${HLA_IN} -o \${HLA_sorted_bam}","\n";
    print HLA "$samtools view \${HLA_sorted_bam} | perl -ne \'\$l=\$_; \$f_q1=\"$f_fq_1\"; \$f_q2=\"$f_fq_2\"; if(\$first==0) { open(OUT1,\">\$f_q1\"); open(OUT2,\">\$f_q2\");  \$first=1;}  \@ss=split(\"\\t\",\$l); \$flag=\$ss[1]; \$cigar=\$ss[5]; if((\$flag & 0x100) || (\$flag & 0x800) || (\$cigar=~/H/)) { next; } \$id=\$ss[0]; \$seq=\$ss[9]; \$q=\$ss[10];  if(\$id=~/\\/1\$/ || (\$flag & 0x40) ) { \$r1=\$id; \$r1=~s/\\/1\$//g; \$seq1=\$seq; \$q1=\$q; } if(\$id=~/\\/2\$/ || (\$flag & 0x80)) { \$r2=\$id; \$r2=~s/\\/2\$//g; \$seq2=\$seq; \$q2=\$q; } if((\$r1 eq \$r2)) { print OUT1 \"\@\",\$r1,\"/1\",\"\\n\"; print OUT1 \$seq1,\"\\n\"; print OUT1 \"+\",\"\\n\"; print OUT1 \$q1,\"\\n\"; print OUT2 \"\@\",\$r1,\"/2\",\"\\n\"; print OUT2 \$seq2,\"\\n\"; print OUT2 \"+\",\"\\n\"; print OUT2 \$q2,\"\\n\";}\'","\n";
    print HLA "  fi\n";
    print HLA "fi\n"; 
    print HLA "if [ -f $f_fq_1 ] && [ -f $f_fq_2 ]","\n"; # input file exist
    print HLA "then\n";
    print HLA "if [ $s_rna -eq 1 ]\n";
    print HLA "then\n";	
    print HLA "if [ -d $dir_hla ]","\n";
    print HLA "then\n";		 
    print HLA "f_optitype_hla=`find $dir_hla -name '*tsv'`","\n";
    print HLA "fi\n";
    print HLA "if [ -z \${f_optitype_hla} ] || [ ! -f \${f_optitype_hla} ]\n";
    print HLA "then\n";
    print HLA "if [ -d $dir_hla ]","\n";
    print HLA "then\n";
    print HLA "rm -rf $dir_hla\n";   
    print HLA "fi\n";
    print HLA "$optitype -i $f_fq_1 $f_fq_2 -r -o $dir_hla"."\n";
    print HLA "fi\n";
    print HLA "else\n";
    print HLA "if [ ! -f \${f_optitype_hla} ]\n";
    print HLA "then\n";
    print HLA "if [ -d $dir_hla ]","\n";
    print HLA "then\n";
	print HLA "rm -rf $dir_hla\n";	
    print HLA "fi\n";  
        print HLA "$optitype -i $f_fq_1 $f_fq_2 -d -o $dir_hla"."\n";	
  	#print HLA "$optitype -i $f_fq_1 $f_fq_2 -c $f_opti_config --dna -v -o $dir_hla"."\n";
	print HLA "  fi\n";	
    print HLA "  fi\n";	
 	print HLA "if [ -d $dir_hla ]","\n";
    print HLA "then\n";
    print HLA "f_optitype_hla=`find $dir_hla -name '*tsv'`","\n";
    print HLA "fi\n";
	print HLA "if [ ! -z \${f_optitype_hla} ] && [  -f \${f_optitype_hla} ]\n";
	print HLA "then\n";
	print HLA "rm \${HLA_sorted_bam}","\n";
	print HLA "rm $f_fq_1","\n";
 	print HLA "rm $f_fq_2","\n"; 
	print HLA "  fi\n";
	print HLA "  fi\n";	
    print HLA "if [ -d $dir_hla ]","\n";
    print HLA "then\n";
    print HLA "f_optitype_hla=`find $dir_hla -name '*tsv'`","\n";
    print HLA "fi\n";
    print HLA "if [ ! -z \${f_optitype_hla} ] && [  -f \${f_optitype_hla} ]\n";	
    print HLA "then\n"; 
    print HLA  " ".$run_script_path_perl."parseHLAresult.pl \${f_optitype_hla} \${HLA_tsv}"."\n";
	print HLA "fi\n";
   	close HLA;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "LSF_DOCKER_ENTRYPOINT=/bin/bash LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>200000] rusage[mem=200000]\" -M 200000000 -a \'docker(fred2/optitype)\' -o $lsf_out -e $lsf_err \'sh $sh_file\'\n";
    system ( $bsub_com );

	}

### calculating binding affinity between epitopes and MHC ##
sub bsub_netmhc{

    my ($step_by_step) = @_;
    if ($step_by_step) {
        $hold_job_file = "";
    }else{
        $hold_job_file = $current_job_file;
    }
    $current_job_file = "j4_bind_".$sample_name.".sh";
    my $lsf_out=$lsf_file_dir."/".$current_job_file.".out";
    my $lsf_err=$lsf_file_dir."/".$current_job_file.".err";
    `rm $lsf_out`;
    `rm $lsf_err`;
    open(MHC, ">$job_files_dir/$current_job_file") or die $!;
    print MHC "#!/bin/bash\n";
    print MHC "f_pep_indel_wt=".$sample_full_path."/$sample_name.indel.vcf.proteome-indel-wt.fasta\n";
    print MHC "f_pep_indel_mut=".$sample_full_path."/$sample_name.indel.vcf.proteome-indel-mut.fasta\n";
    print MHC "f_pep_snv_wt=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-wt.fasta\n";
    print MHC "f_pep_snv_mut=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-mut.fasta\n";
    print MHC "f_pep_snv_mut_v1=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-mut.v1.fasta\n";		
    print MHC "f_pep_snv_mut_v2=".$sample_full_path."/$sample_name.snp.vcf.proteome-snv-mut.v2.fasta\n";		
    print MHC "HLA_tsv=".$sample_full_path."/HLA_alleles.tsv\n";
    print MHC "f_netMHC_result=".$sample_full_path."/netMHC4.0.out.append.txt\n";
    print MHC "f_netMHC_result_snv=".$sample_full_path."/netMHC4.0.out.append.snv.txt\n";
    print MHC "f_netMHC_result_indel=".$sample_full_path."/netMHC4.0.out.append.indel.txt\n";
    print MHC "f_out=".$sample_full_path."/result_neoantigen\n";
    print MHC " ".$run_script_path_perl."generate_mut_peptide_snv.pl $db_ref_bed \${f_pep_snv_wt} \${f_pep_snv_mut} \${f_pep_snv_mut_v1}\n";
    print MHC " ".$run_script_path_perl."remove_duplicate_mut_peptide_snv.pl \${f_pep_snv_mut_v1} \${f_pep_snv_mut_v2}\n";
    print MHC  " ".$run_script_path_python."runNetMHC4.py -a \${HLA_tsv} -f \${f_pep_snv_mut_v2} -p 8,9,10,11 -o $sample_full_path -n $netMHC -v $f_allele"."\n";
    print MHC "mv \${f_netMHC_result} \${f_netMHC_result_snv}\n"; 
    print MHC  " ".$run_script_path_python."runNetMHC4.py -a \${HLA_tsv} -f \${f_pep_indel_mut} -p 8,9,10,11 -o $sample_full_path -n $netMHC -v $f_allele"."\n";
    print MHC "mv \${f_netMHC_result} \${f_netMHC_result_indel}\n";

### check if netMHC prediction is successfully finished, indel
    print MHC 'if [ -s $f_pep_indel_mut ] ',"\n"; # file exist  
    print MHC "then\n";
    print MHC 'if [ -f $f_netMHC_result_indel ] ',"\n"; # file exist
    print MHC "then\n";
    print MHC ' grep "Error" ${f_netMHC_result_indel}',"\n";
    print MHC ' CHECK=$?',"\n";
    print MHC ' while [ ${CHECK} -eq 0 ] ',"\n"; # grep success, file not finish
    print MHC " do\n";
    print MHC  " ".$run_script_path_python."runNetMHC4.py -a \${HLA_tsv} -f \${f_pep_indel_mut} -p 8,9,10,11 -o $sample_full_path -n $netMHC -v $f_allele"."\n";
    print MHC "mv \${f_netMHC_result} \${f_netMHC_result_indel}\n";
    print MHC '     grep "Error" ${f_netMHC_result_indel}',"\n";
    print MHC '     CHECK=$?',"\n";
    print MHC " done\n";
    print MHC "     fi\n";
    print MHC "     fi\n";
### check if netMHC prediction is successfully finished, snv  
   
    print MHC 'if [ -s $f_pep_snv_mut_v2 ] ',"\n"; # file exist  
    print MHC "then\n"; 
    print MHC 'if [ -f $f_netMHC_result_snv ] ',"\n"; # file exist
    print MHC "then\n";
    print MHC ' grep "Error" ${f_netMHC_result_snv}',"\n";
    print MHC ' CHECK=$?',"\n";
    print MHC ' while [ ${CHECK} -eq 0 ] ',"\n"; # grep success, file not finish
    print MHC " do\n";
    print MHC  " ".$run_script_path_python."runNetMHC4.py -a \${HLA_tsv} -f \${f_pep_snv_mut_v2} -p 8,9,10,11 -o $sample_full_path -n $netMHC -v $f_allele"."\n";
    print MHC "mv \${f_netMHC_result} \${f_netMHC_result_snv}\n";
    print MHC '     grep "Error" ${f_netMHC_result_snv}',"\n";
    print MHC '     CHECK=$?',"\n";
    print MHC " done\n";
    print MHC "     fi\n";
    print MHC "     fi\n";
    close MHC;

    my $sh_file=$job_files_dir."/".$current_job_file;

    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>20000] rusage[mem=20000]\" -M 20000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
   #    system ( $bsub_com );

#    $bsub_com = "LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q dinglab -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(scao/dailybox)\' -o $lsf_out -e $lsf_err sh $sh_file\n";

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
    $current_job_file = "j5_parsebind_".$sample_name.".sh";
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
    print PMHC "f_netMHC_result_snv=".$sample_full_path."/netMHC4.0.out.append.snv.txt\n";
	print PMHC "f_netMHC_result_indel=".$sample_full_path."/netMHC4.0.out.append.indel.txt\n";
    print PMHC "f_snv_out=".$sample_full_path."/$sample_name.neoantigen.snv.tsv\n";
	print PMHC "f_indel_out=".$sample_full_path."/$sample_name.neoantigen.indel.tsv\n";
	print PMHC "f_snv_sum=".$sample_full_path."/$sample_name.neo.snv.summary\n";
	print PMHC "f_indel_sum=".$sample_full_path."/$sample_name.neo.indel.summary\n";
#	print PMHC "f_min=".$sample_full_path."/$sample_name.neo.summary.min\n";
  	print PMHC  " ".$run_script_path_perl."parseNetMHC4result.pl \${f_netMHC_result_snv} \${f_indel_wt_fa} \${f_snv_wt_fa} \${f_snv_out}"."\n";
    print PMHC	" ".$run_script_path_perl."reportSummary.pl \${f_snv_out} \${f_indel} \${f_snv} \${f_indel_mut_fa} \${f_snv_mut_fa} \${f_snv_sum}"."\n";
	print PMHC  " ".$run_script_path_perl."parseNetMHC4result.pl \${f_netMHC_result_indel} \${f_indel_wt_fa} \${f_snv_wt_fa} \${f_indel_out}"."\n";
    print PMHC  " ".$run_script_path_perl."reportSummary.pl \${f_indel_out} \${f_indel} \${f_snv} \${f_indel_mut_fa} \${f_snv_mut_fa} \${f_indel_sum}"."\n";
	close PMHC;
    my $sh_file=$job_files_dir."/".$current_job_file;
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>20000] rusage[mem=20000]\" -M 20000000 -a \'docker(ubuntu)\' -o $lsf_out -e $lsf_err sh $sh_file\n";

    #$bsub_com = "bsub -q dinglab -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
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
    $current_job_file = "j6_final_report_".$working_name.".sh";
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
	print REP " ".$run_script_path_perl."generate_report_summary_2.pl $run_dir"."\n";
	close REP;
    #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    #system ( $bsub_com );

    my $sh_file=$job_files_dir."/".$current_job_file;
   # $bsub_com = "bsub -q dinglab -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w $hold_job_file -o $lsf_out -e $lsf_err sh $sh_file\n"; 
    $bsub_com = "bsub -g /$compute_username/$group_name -q $q_name -n 1 -R \"select[mem>20000] rusage[mem=20000]\" -M 20000000 -a \'docker(ubuntu)\' -o $lsf_out -e $lsf_err sh $sh_file\n";
    #$bsub_com = "bsub -q dinglab -n 1 -R \"select[mem>30000] rusage[mem=30000]\" -M 30000000 -a \'docker(registry.gsc.wustl.edu/genome/genome_perl_environment)\' -w \"$hold_job_file\" -o $lsf_out -e $lsf_err sh $sh_file\n";
   #$bsub_com = "bsub < $job_files_dir/$current_job_file\n";
    system ( $bsub_com );
    
}




