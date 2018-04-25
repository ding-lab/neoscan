## 12-03-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will get min result 


OUT


die $usage unless @ARGV == 2;

my($f_sum_in,$f_sum_out)=@ARGV;

open(IN,"<$f_sum_in");
#open(IN_smut,"<$f_snv_mut");
open(OUT,">$f_sum_out"); 
my %min_b=();

while(<IN>)
	{
		my $line=$_;
		chomp($line);
		my @temp=split("\t",$line); 
		my $id=$temp[0]."-".$temp[1]."-".$temp[6]; 
		if(! defined $min_b{$id})
		{
		 $min_b{$id}=$line; 	
		}		
		else 
			{
				my @temp2=split("\t",$min_b{$id});
				#print $temp2[9],"\n";
				if($temp[9]<$temp2[9]) { $min_b{$id}=$line; }
			}
	}	

foreach my $id (sort keys %min_b)
	{
	 print OUT $min_b{$id},"\n"; 

	}
close(IN);
close(OUT);
