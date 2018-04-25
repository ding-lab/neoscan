## 12-09-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will parse mapped read

OUT

die $usage unless @ARGV == 2;

my($f_map_in,$f_fq_out)=@ARGV;
my %num=(); 
my $nmax=0;
open(IN,"<$f_map_in");
open(OUT,">$f_fq_out");
 
while(<IN>)
	{
		my $line=$_; 
		chomp($line); 
		if($line=~/^\@SQ/) { next; }
		my @temp=split("\t",$line); 
		my $cigar=$temp[5]; 
		
		if($cigar=~/(\d+)M/) {  	
		$num{$1}++; }	
	}
close IN;

foreach my $n (sort { $num{$b} <=> $num{$a} } keys %num) 	
	{
		$nmax=$n;
		last;  	
	}

#print $nmax,"\n"; 
my %quality=();
my %seq=(); 

open(IN,"<$f_map_in");

while(<IN>)
	{
		my $line=$_; 
		chomp($line); 
		if($line=~/^\@SQ/) { next; }
        my @temp=split("\t",$line);
        my $cigar=$temp[5];
#		my $id=$temp[0]; 
		my $flag=$temp[1]; 
		my $nm=$temp[12]; 	
#\$r2=\$id; \$r2=~s/\\/2\$//g; \$seq2=\$seq; \$q2=\$q; } if((\$r1 eq \$r2)) { print OUT1 \"\@\",\$r1,\"/1\",\"\\n\"; print OUT1 \$seq1,\"\\n\"; print OUT1 \"+\",\"\\n\"; print OUT1 \$q1,\"\\n\";
        if($cigar=~/^$nmax\M/ && $nm=~/NM:i:0/) 
		{
			my $id=$temp[0];
			if($flag & 0x40) 
			{
				my $rid=$id."#1";  
				$quality{$rid}=$temp[10]; 
				$seq{$rid}=$temp[9];
			}
			if($flag & 0x80)
			{
				my $rid=$id."#2";	
				$quality{$rid}=$temp[10];
            	$seq{$rid}=$temp[9];			
        	}
		}
	}

foreach my $id (sort keys %seq) 
	{
		print OUT "\@",$id,"\n"; 
		print OUT $seq{$id},"\n";
		print OUT "+","\n";
		print OUT $quality{$id},"\n"; 
	}

