#!/usr/bin/perl

use strict;

my $usage = "
This script will read corresponding files for each sample and generate a report which 
perl $0 <run folder>
<run folder> = full path of the folder holding files for this sequence run
";
die $usage unless scalar @ARGV == 1;
my ($dir) = @ARGV;

my @temp = split("/", $dir);
my $run_name = pop @temp;
my $outFile = $dir."/Neoantigen_Report_".$run_name;
my %neo=();

open (OUT, ">$outFile") or die "can not open file $outFile!\n";

#&generate_table($dir);

print OUT "Sample","\t","Chr","\t","Pos","\t","WT","\t","Var","\t","Gene","\t","Type","\t","Refseq","\t","HLA","\t","Neoantigen","\t","Binding affinity","\n";

foreach my $s (`ls $dir`) { 

my $str=$s; 
chomp($str); 
my $dirs=$dir."/".$str;
my $f_neo_snv=$dirs."/".$str.".neo.snv.summary"; 
my $f_neo_ind=$dirs."/".$str.".neo.indel.summary"; 
if(-f $f_neo_snv) 
{ 
foreach my $l (`cat $f_neo_snv`) 
{ 
my $ltr=$l; 
chomp($ltr); 
print OUT $str,"\t",$ltr,"\n"; 
} 
}  

if(-f $f_neo_ind) 
{ 
foreach my $l (`cat $f_neo_ind`) 
{ 
my $ltr=$l; 
chomp($ltr); print OUT $str,"\t",$ltr,"\n"; 
} 
}
}
close OUT; 
