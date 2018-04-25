
#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will extract the supporting reads Step 1:

OUT

die $usage unless @ARGV == 3;

my($f_ind_wt,$f_ind_mut,$f_out)=@ARGV;

open(IN1,"<$f_ind_wt"); 
open(IN2,"<$f_ind_mut"); 
open(OUT,">$f_out"); 
my %wt_seq=(); 
my $rd_id; 

while(<IN1>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else { $wt_seq{$rd_id}=$line;  }
    }
close IN1; 


while(<IN2>)
	{
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else { 
			my @temp=split("-",$rd_id); my $id=$temp[0];   
			for(my $i=0;$i<length($line);$i++)
			{
				my $r1=substr($line,$i,1); 
				my $r2=substr($wt_seq{$id},$i,1); 
				if($r1 ne $r2) { print OUT $rd_id,"\t",$i,"\n"; last; }
			}
		}
    }

