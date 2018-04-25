## 12-03-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will parse HLA result 


OUT


die $usage unless @ARGV == 2;

my($f_hla_in,$f_hla_out)=@ARGV;

open(IN_hla,"<$f_hla_in");
#open(IN_smut,"<$f_snv_mut");
open(OUT,">$f_hla_out"); 

my $print_out=0;
while(<IN_hla>)
	{
		my $line=$_;
		chomp($line);
		
		if($line=~/Prediction \#1/)
		{
			$print_out=1; 
			next; 
		}
			
		if($print_out==1) 
		{
			if($line=~/\*/ && (! ($line=~/Prediction/))) 
			{ 
				my @temp=split(/\,/,$line);  
				my $hla=$temp[0];
			    $hla=~ s/^\s+|\s+$//g;
				$hla=~ s/\*//g;
				$hla=~ s/P//g;  	
				print OUT "HLA-".$hla,"\t",$temp[3],"\n"; }
			else { $print_out=0; }
		}
		}	

close(IN_hla);
close(OUT);
