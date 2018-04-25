## 12-09-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will extract the supporting reads Step 1:

OUT

die $usage unless @ARGV == 5;

my($f_trans_fa,$f_snv,$f_indel,$f_mapped_fa,$f_out)=@ARGV;

my %mut_seq=();
my %snv_wt_seq=();
my %indel_wt_seq=(); 
my %mapped_seq=(); 
  
open(IN1,"<$f_trans_fa");
open(IN2,"<$f_snv");
if(-f $f_indel)
{
open(IN3,"<$f_indel"); 
}

open(IN4,"<$f_mapped_fa");  

open(OUT,">$f_out");

my $rd_id; 

while(<IN1>)
	{
		my $line=$_; 
		chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else { $mut_seq{$rd_id}=$line;	}
		
	}
close IN1; 


while(<IN2>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else { $snv_wt_seq{$rd_id}=$line; }

    }
close IN2;

if(-f $f_indel) 
{
while(<IN3>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else { $indel_wt_seq{$rd_id}=$line; }

    }
close IN3;
}

my $write_rd=0;
my $cc=0; 
while(<IN4>)
    {
        my $line=$_;
        chomp($line);
		my $int4=$cc-int($cc/4)*4;
        if($line=~/^\@/ && $int4==0) { $write_rd=1; $rd_id=$line; $rd_id=~s/^\@//g; $cc++; next;  }
        if($write_rd==1) { $mapped_seq{$rd_id}=$line; $write_rd=0;  }
	    $cc++;
    }
close IN4;

foreach my $id1 (sort keys %mapped_seq)
	{
		 my $seq=$mapped_seq{$id1};
		 #print $seq,"\n"; 
         my $seq_rc= reverse $seq;
         $seq_rc=~tr/ATCG/TAGC/;
		 my $index_novel=1;
		 #print $id1,"\n";
		  
		 foreach my $id2 (sort keys %snv_wt_seq)
			{
				my $seq2=$snv_wt_seq{$id2}; 
         		if($seq2=~/$seq/ || $seq2=~/$seq_rc/)
           		{
				 $index_novel=0; 
				 last; 
				}
			}
		  if(-f $f_indel) 
			{
		   foreach my $id2 (sort keys %indel_wt_seq)
            {
                my $seq2=$indel_wt_seq{$id2};
                if($seq2=~/$seq/ || $seq2=~/$seq_rc/)
                {
                 $index_novel=0;
                 last;
                }
            }
			}

		if($index_novel==1) 
		{
			$index_novel=0; 
			foreach my $id2 (sort keys %mut_seq)
            {
				#print $seq2,"\n";
                my $seq2=$mut_seq{$id2};
				#print $seq2,"\n"; <STDIN>;
                if($seq2=~/$seq/ || $seq2=~/$seq_rc/)
                {
                 $index_novel=1;
                 last;
                }
            }
		  }
		 # print $id1,"\t",$seq,"\t",$index_novel,"\n";
		  #<STDIN>;
		  if($index_novel==1) { print OUT ">$id1","\n"; print OUT $seq,"\n"; }

	}

close OUT; 
