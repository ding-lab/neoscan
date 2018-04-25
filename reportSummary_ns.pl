## generate mutant and wild-type protein sequence for snv ##
## 12-03-2016 ##

#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will generate the summary for potentianl neoantigen with strong binding to MHC 


OUT


die $usage unless @ARGV == 6;

my($f_neo,$f_indel_vcf,$f_snv_vcf,$f_indel_fa,$f_snv_fa,$f_sum)=@ARGV;

open(IN_neo,"<$f_neo");
if(-s $f_indel_vcf)
{
open(IN_ivcf,"<$f_indel_vcf");
}

if(-s $f_snv_vcf)
{
open(IN_svcf,"<$f_snv_vcf");
}

if(-s $f_indel_fa)
{
open(IN_ifa,"<$f_indel_fa");
}

if(-s $f_snv_fa)
{ open(IN_sfa,"<$f_snv_fa"); }

open(OUT,">$f_sum");

my %SB=();
my %imut=();
my %smut=();
my %imut_neo=();
my %smut_neo=();

while(<IN_neo>)
	{
	my $line=$_; 
	chomp($line); 
	my @temp=split("\t",$line); 
	$SB{$temp[0]}{$temp[1]}=$temp[2]; 
	}

my $id1;
my $id2; 

if(-s $f_indel_fa)
{
my $np;
my $pos;

while(<IN_ifa>)
	{
		my $l=$_;
		chomp($l);
		if($l=~/^>NP/)
                {
                    my @temp=split("-",$l);
                    $np=$temp[0];
                    $np=~s/^>//g;
                    $np=~s/-indel//g;
                    #my @temp2=split("-",$temp[4]);
                   # $temp[3]+=1;
                    $pos=$temp[2]."-".$temp[3];
                    $pos=~s/^chr//g;
                    $pos=~s/^chr//g;
                }
      else {  $imut{$pos}{$np}=$l;
			  #print $pos,"\t",$np,"\n";
			  #<STDIN>;
               # print "pep indel mut","\t",$np,"\t",$pos,"\n";
                    }
	#	if($line=~/^>/)
	#	{
#
#		    my @temp=split /\s+/, $line;
#			$id1=$temp[4];
#			my @temp2=split("-",$temp[0]); 
#			$id2=$temp2[0];
#			$id2=~s/\>//g; 
#		}	
#		else 	
#		{
#			my @temp3=split("-",$id1); 
#			my $id3=$temp3[0];
#			$id3=~s/\(DEL:chr//g; 	
#			my $id4=$temp3[1]+1;
#			my $id34=$id3."-".$id4;
			#print $id34,"\t",$id2,"\n";
			#<STDIN>;
#			$imut{$id34}{$id2}=$line;
#		}	
	}	
}

if(-s $f_snv_fa)
{
my $np;
my $pos; 

while(<IN_sfa>)
    {
        my $l=$_;
        chomp($l);

                if($l=~/^>NP/)
                {
                    my @temp=split("-",$l);
                    $np=$temp[0];
                    $np=~s/^>//g;
                    $np=~s/-VAR//g;
                   # my @temp2=split("-",$temp[1]);
                   # $temp[2]+=1;
                    $pos=$temp[1]."-".$temp[2];
                    $pos=~s/^chr//g;
                }
                else {
                     $smut{$pos}{$np}=$l;
					 #print $pos,"\t",$np,"\n";
					#<STDIN>;
                  #   print "pep snv mut","\t",$np,"\t",$pos,"\n";
                     }

       # if($line=~/^>/)
       # {

        #    my @temp=split /\s+/, $line;
        #    $id1=$temp[4];
        #    my @temp2=split("-",$temp[0]);
        #    $id2=$temp2[0];
        #    $id2=~s/\>//g;
        #}
        #else
        #{
         #   my @temp3=split("-",$id1);
         #   my $id3=$temp3[0];
         #   $id3=~s/\(VAR:chr//g;
         #   my $id4=$temp3[1]+1;
         #   my $id34=$id3."-".$id4;
         #   $smut{$id34}{$id2}=$line;
        #}
    }

}

foreach my $var (sort keys %imut) 
	{
		foreach my $np (sort keys %{$imut{$var}}) 
		{
			foreach my $hla (sort keys %SB)
			{
				foreach my $pep (sort keys %{$SB{$hla}})
				{
					if($imut{$var}{$np}=~/$pep/) 
					{
					 	$imut_neo{$var}{$np}{$hla}{$pep}=$SB{$hla}{$pep}; 	
					}
				}
			}  
		}
	}


foreach my $var (sort keys %smut)
    {
        foreach my $np (sort keys %{$smut{$var}})
        {
            foreach my $hla (sort keys %SB)
            {
                foreach my $pep (sort keys %{$SB{$hla}})
                {
					#print $hla,"\t",$pep,"\n";
					#print $smut{$var}{$np},"\n";
                    if($smut{$var}{$np}=~/$pep/)
                    {
						#print $var,"\t",$hla,"\t",$pep,"\n";
                    	#print $smut{$var}{$np},"\n";
                        $smut_neo{$var}{$np}{$hla}{$pep}=$SB{$hla}{$pep};
                    }
                }
            }
        }
    }

if(-s $f_indel_vcf)
{
while(<IN_ivcf>)
	{
		my $line=$_; 
		chomp($line); 
		my @temp=split("\t",$line); 
	 	my $id=$temp[0]."-".$temp[1]; 
		if(defined $imut_neo{$id})
		{
			foreach my $np (sort keys %{$imut_neo{$id}})
			{
				foreach my $hla (sort keys %{$imut_neo{$id}{$np}})
				{
					foreach my $pep (sort keys %{$imut_neo{$id}{$np}{$hla}})
					{
					 print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t","INDEL","\t",$np,"\t",$hla,"\t",$pep,"\t",$imut_neo{$id}{$np}{$hla}{$pep},"\n";
					}
				}
			}
		}		
	}

}

if(-s $f_snv_vcf)
{
while(<IN_svcf>)
    {
        my $line=$_;
        chomp($line);
        my @temp=split("\t",$line);
        my $id=$temp[0]."-".$temp[1];
        if(defined $smut_neo{$id})
        {
            foreach my $np (sort keys %{$smut_neo{$id}})
            {
                foreach my $hla (sort keys %{$smut_neo{$id}{$np}})
                {
                    foreach my $pep (sort keys %{$smut_neo{$id}{$np}{$hla}})
                    {
					print OUT $temp[0],"\t",$temp[1],"\t",$temp[2],"\t",$temp[3],"\t",$temp[4],"\t","SNV","\t",$np,"\t",$hla,"\t",$pep,"\t",$smut_neo{$id}{$np}{$hla}{$pep},"\n";
                    }
                }
            }
        }
    }
 }

close OUT;

