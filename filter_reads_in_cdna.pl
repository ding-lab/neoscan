## 12-09-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will parse mapped read

OUT

die $usage unless @ARGV == 3;

my($f_map_in,$f_fq_in,$f_fq_out)=@ARGV;
my %num=(); 
my $nmax=0;
open(IN,"<$f_map_in");
open(IN_fq,"<$f_fq_in"); 
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
my %filter_id=(); 

open(IN,"<$f_map_in");

while(<IN>)
	{
	my $line=$_; 
	chomp($line); 
	if($line=~/^\@SQ/) { next; }
        my @temp=split("\t",$line);
        my $cigar=$temp[5];
	my $flag=$temp[1]; 
	my $nm=$temp[12]; 	
        if($cigar=~/^$nmax\M/ && $nm=~/NM:i:0/) 
		{
		my $id=$temp[0];
            	$filter_id{$id}=1;			
		}
	}

close IN; 
my $write_up=1; 
my $cc=0; 
while(<IN_fq>)
	{
	 my $line=$_; 
	 chomp($line);
	 my $int4=$cc-int($cc/4)*4; 

	 if($line=~/^\@/ && $int4==0)
	 {
	  my $id=$line; 
	  $id=~s/^\@//g;
	  if(defined $filter_id{$id}) { $write_up=0;}
	  else { $write_up=1; print OUT $line,"\n"; }  
	 }
	 else 
	  { 
	   if($write_up==1) { print OUT $line,"\n"; } 
	  }
	  $cc++; 
	}

close IN_fq; 
	
close OUT;	
