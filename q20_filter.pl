#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will extract high quality read q20:

OUT

die $usage unless @ARGV == 4; 

my($f_db,$f_fa,$f_fq,$f_out)=@ARGV;
open(DB,"<$f_db"); 
open(IN1,"<$f_fa"); 
open(IN2,"<$f_fq");
open(OUT,">$f_out");

my %qscore=(); 
my %quality=(); 

while(<DB>)
	{
	    my $line=$_;
        chomp($line);
		my @temp=split("\t",$line); 
		$qscore{$temp[0]}=$temp[1]; 
	}

my $rd_id;
 
while(<IN2>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\@/) { $rd_id=$line; $rd_id=~s/^\@//g; next;  }
        else { $quality{$rd_id}=$line; }
    }
close IN2;


while(<IN1>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else {
		my $print_out=1;  
		my $qstr=$quality{$rd_id}; 
		for(my $i=0;$i<length($qstr);$i++)
		{
			my $q=substr($qstr,$i,1); 
			if(defined $qscore{$q} && $qscore{$q}<20) 
			{
			 $print_out=0; 
			last; 
			}
		}
		if($print_out==1) { print OUT ">$rd_id\n"; print OUT $line,"\n"; }
		}
    }
close IN1;
close OUT;  
