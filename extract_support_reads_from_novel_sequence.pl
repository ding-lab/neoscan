## 12-09-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will extract the supporting reads

OUT

die $usage unless @ARGV == 4;

my($f_trans_fa,$f_novel_fa,$f_pos,$f_out)=@ARGV;

my %novelseq=();
my %pos_first_diff=(); 

open(IN1,"<$f_trans_fa");
open(IN2,"<$f_novel_fa"); 
open(IN3,"<$f_pos"); 
open(OUT,">$f_out");

my $rd_id; 

while(<IN3>)
	{
		my $line=$_; 
        chomp($line); 
		my @temp=split("\t",$line); 
		$pos_first_diff{$temp[0]}=$temp[1]; 	
	}

while(<IN2>)
	{
		my $line=$_; 
		chomp($line); 
		if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
		else { $novelseq{$rd_id}=$line; }
	}

close IN2;

while(<IN1>)
	{
		my $line=$_; 
		chomp($line); 
		if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^\>//g;   }
        else {
				my %support_reads=();
				my %support_pos=(); 
				my $n_count=0; 
				foreach my $id (sort keys %novelseq) 
				{
					my $seq=$novelseq{$id}; 
					my $seq_rc= reverse $seq;
                    $seq_rc=~tr/ATCG/TAGC/;
					if($line=~/$seq/g)
					{
						if(pos($line)-length($seq)<$pos_first_diff{$rd_id} && pos($line)>$pos_first_diff{$rd_id})
                        {
					 	$support_reads{$id}=$seq; 
						$support_pos{$id}=pos($line)-length($seq);  $support_pos{$id}.="-orig";  
						$n_count++; }
					} 			
					else
					{
						if($line=~/$seq_rc/g) 
						{ 
						if(pos($line)-length($seq)<$pos_first_diff{$rd_id} && pos($line)>$pos_first_diff{$rd_id})
						{
						$support_reads{$id}=$seq_rc; $support_pos{$id}=pos($line)-length($seq);  
						$support_pos{$id}.="-revcomp"; $n_count++;
						} 
						}
					}	
				}
				if($n_count>0) 
				{ 
					print OUT $rd_id,"\t",$n_count,"\n"; foreach my $id (sort keys %support_reads) {
					print OUT $id,"\t",$support_reads{$id},"\t",$support_pos{$id},"\n"; }
				}
			}	
	}

close IN1; 
close OUT; 
