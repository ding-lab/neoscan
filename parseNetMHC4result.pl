## 12-03-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will report neoantigen with strong binding affinity to MHC 


OUT


die $usage unless @ARGV == 4;

my($f_mhc_in,$f_indel_wt,$f_snv_wt,$f_neo)=@ARGV;

open(IN_mhc,"<$f_mhc_in");
open(IN_iwt,"<$f_indel_wt");
#open(IN_imut,"<$f_indel_mut");
open(IN_swt,"<$f_snv_wt");
#open(IN_smut,"<$f_snv_mut");
open(OUT,">$f_neo"); 
my %SB=();
#my %iwt=();
my %wt=();

my $id;

if(-s $f_indel_wt)
{ 
while(<IN_iwt>)
	{
		my $line=$_;
		chomp($line);

		if($line=~/^>/)
		{
			my @temp=split(/\-/,$line);
			$id=$temp[0];
			$id=~s/\>//g; 
		}	
		else 	
		{
			#print "indel","\t",$id,"\n";
			$wt{$id}=$line; 
		}	
	}	
}

while(<IN_swt>)
    {
        my $line=$_;
        chomp($line);
    	#print $line,"\n";
	    if($line=~/^>/)
        {
            #my @temp=split(/\s+/,$line);
            $id=$line;  		
			$id=~s/\>//g; 
        }       
        else
        {
			#<STDIN>;
            $wt{$id}=$line;
        }
    }


while(<IN_mhc>)
    {
        my $line=$_;
        chomp($line);
        my @temp=split /\s+/, $line;
        if($temp[11]=~/^NP/)
        {
            if($temp[13]<=500)
            {

               # print $temp[2],"\t",$temp[10],"\t",$temp[11],"\t",$temp[13],"\n";
                my @temp2=split(/\-/,$temp[11]);
				my $pep=$temp[10];
                my $id=$temp2[0];
				my $wt_type=0;
				## remove if in wild-type
				foreach my $id2 (sort keys %wt) 
				{
				if($wt{$id2}=~/$pep/) { $wt_type=1; }
				}
	#			if($temp2[1]=~/^i/) { if($iwt{$id}=~/$pep/) { $wt_type=1; }}
				if($wt_type==0) 
				{
					$SB{$temp[2]}{$temp[10]}=$temp[13]; 
				}
               #$SB{$id}{$temp[10]}{$temp[2]}=$temp[13];
            }
        }
    }

foreach my $hla (sort keys %SB) 
	{
		foreach my $p (sort keys %{$SB{$hla}}) 
		{
			print OUT $hla, "\t", $p, "\t", $SB{$hla}{$p},"\n"; 
		}
	}	
