#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will generate the peptide database with mutations

perl remove_duplicate_mut_peptide_snv.pl f_snv_mut f_out

OUT

die $usage unless @ARGV == 2;

my($f_snv_mut,$f_out)=@ARGV;
open(SNVMUT,"<$f_snv_mut");
open(OUT,">$f_out");

my $rd_id; 
my %peps=();

while(<SNVMUT>)
    {
        my $line=$_;
        chomp($line);
        if($line=~/^\>/) { $rd_id=$line;  }

        else { if(! defined $peps{$line}) 
				{ print OUT $rd_id,"\n";  print OUT $line,"\n"; $peps{$line}=1; }
		  	 }
    }

close SNVMUT;
close OUT; 

