## 12-03-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will parse HLA result 


OUT


die $usage unless @ARGV == 2;

my($f_hla_in,$f_hla_out)=@ARGV;

#open(IN_smut,"<$f_snv_mut");
open(OUT,">$f_hla_out"); 

foreach my $l (`cat $f_hla_in`) 
{ 
 	my $ltr=$l; 
	chomp($ltr); 
	if($ltr=~/Reads/) { next; } 
	else { 
	my @temp=split("\t",$ltr); 
	for(my $i=1;$i<=6;$i++) 
	{ 
	my $hla=$temp[$i]; $hla=~s/\*//g; 
	print OUT "HLA-".$hla,"\t","1000","\n";  }
}
}
close(OUT);
