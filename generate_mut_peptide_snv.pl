#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will generate the peptide database with mutations

perl generate_mut_peptide.pl  f_bed f_snv_wt f_snv_mut f_out

OUT

die $usage unless @ARGV == 4;

my($f_bed,$f_snv_wt,$f_snv_mut,$f_out)=@ARGV;
open(SNVWT,"<$f_snv_wt");
#open(INDELWT,"<$f_indel_wt");
open(BED,"<$f_bed"); 
open(SNVMUT,"<$f_snv_mut");
#open(INDELMUT,"<$f_indel_mut");
open(OUT,">$f_out");

my %wt_seq=();
my %snv_pos=();
my $rd_id; 
my %strandness=(); 

while(<BED>)
{
	my $l=$_; 
	chomp($l); 
	my @temp=split("\t",$l); 
	$strandness{$temp[3]}=$temp[5]; 
}

while(<SNVWT>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else { $wt_seq{$rd_id}=$line;  }
    }
close SNVMT;


while(<SNVMUT>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }
        else {
            my @temp=split("-",$rd_id); my $id=$temp[0];
			print $id,"\n"; 
            for(my $i=0;$i<length($line);$i++)
            {
                my $r1=substr($line,$i,1);
                my $r2=substr($wt_seq{$id},$i,1);
				print $i,"\n";
				print $r1,"\n"; 
				print $r2,"\n";
				<STDIN>;
                if($r1 ne $r2) { $snv_pos{$rd_id}{$i}=1;  }
            }
        }
    }

close SNVMUT; 

open(SNVMUT,"<$f_snv_mut");

while(<SNVMUT>) 
	{
	    my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line; $rd_id=~s/^>//g;  }	
		else { 
		my $pep=$line;
		my @temp=split("chr",$rd_id); 	
		my $n_temp=scalar @temp; 
		my $ind=0; 
		foreach my $p (sort {$a<=>$b} keys %{$snv_pos{$rd_id}}) 	
		{
			#print $p,"\n"; 
			my $start=$p-12; 
			my $end=$p+12; 
			if($start<0) { $start=0; }
			if($end>length($pep)-1) { $end=length($pep)-1; } 	
			my $pepstr=substr($pep,$start,$end-$start); 	
			my $np=$temp[0]; $np=~s/-VAR://g;
			#print $np,"\t",$strandness{$np},"\n"; 
			#<STDIN>;
			if($strandness{$np} eq "-") 
			{
				print OUT ">",$temp[0],"chr",$temp[$n_temp-1-$ind],"\n"; 
			}
			else 
			{
			print OUT ">",$temp[0],"chr",$temp[$ind+1],"\n"; 
			}
			$ind++; 			
			print OUT $pepstr,"\n"; 	
		}												
		}		
	}

close SNVMUT; 
