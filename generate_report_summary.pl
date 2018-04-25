#!/usr/bin/perl
use strict;

my $usage = "
This script will read corresponding files for each sample and generate a report which 
contains variant, peptides for wild-type and mutants, RNA supporting reads and MHC binding score.  
perl $0 <run folder> <program version>
<run folder> = full path of the folder holding files for this sequence run
";
die $usage unless scalar @ARGV == 2;
my ( $dir, $version ) = @ARGV;

my @temp = split("/", $dir);
my $run_name = pop @temp;
my $outFile = $dir."/Analysis_Report_".$run_name;
my %mut=();
my %pepmut=();
my %pepwt=();
my %suprna=(); 
my %neo=();
my %ref_sup=();

open (OUT, ">$outFile") or die "can not open file $outFile!\n";

my ($wkday,$month,$day,$time,$year) = split(/\s+/, localtime);

#print OUT "Neoantigen-scan ${version}; Processing date: $wkday-$month-$year\n";

#my $c = "**************************************************************************\n";
my $c2 = "#########################################################################\n\n";
#print OUT $c;

#print OUT "Summary:\n\n";
&generate_table($dir);

print OUT "Sample","\t","Refseq_ID","\t","Pos_Mut","\t","Info_Mut","\t","RNA supporting (WT)","\t","RNA supporting (MUT)","\t","MHC binding infor","\t","1stDiff_Mut","\t","Pep_Mut","\t","Pep_WT","\n";
foreach my $s (sort keys %pepmut)
{
	foreach my $np (sort keys %{$pepmut{$s}})
	{
		foreach my $pos (sort keys %{$pepmut{$s}{$np}})
		{
			if(defined $suprna{$s}{$np}{$pos})
			{
				my $p1=$pepmut{$s}{$np}{$pos};
				my $p2=$pepwt{$s}{$np};
				my $diffpos=0; 
				for(my $i=0;$i<length($p1);$i++)
				{
					if(substr($p1,$i,1) ne substr($p2,$i,1))
					{
					$diffpos=$i; last; 
					}
				}

				my $n_ref=0;
			 	if(defined $ref_sup{$s}{$np}{$pos}) { $n_ref=$ref_sup{$s}{$np}{$pos}; }
				else { $n_ref="NA"; }

				if(defined $neo{$s}{$np}{$pos})
				{ 
					print OUT $s,"\t",$np,"\t",$pos,"\t",$mut{$s}{$pos},"\t",$n_ref,"\t",$suprna{$s}{$np}{$pos},"\t",$neo{$s}{$np}{$pos},"\t",$diffpos,"\t",$pepmut{$s}{$np}{$pos},"\t",$pepwt{$s}{$np},"\n";
				}
				else { print OUT $s,"\t",$np,"\t",$pos,"\t",$mut{$s}{$pos},"\t",$n_ref,"\t",$suprna{$s}{$np}{$pos},"\t","NO","\t",$diffpos,"\t",$pepmut{$s}{$np}{$pos},"\t",$pepwt{$s}{$np},"\n"; }
			}
			}
		}
}
##### generate_table####

sub generate_table {
	my ($dir) = @_;
	opendir(DH, $dir) or die "Can not open dir $dir!\n";
	my @files = readdir DH;
	foreach my $name (sort {$a cmp $b} @files) {
		if (!($name =~ /\./) && !($name=~/^Analysis/)) 
		{

		print $name,"\n";
	#	<STDIN>;

		my $full_path = $dir."/".$name;
		my $fIvcf=$full_path."/".$name.".indel.vcf";
		my $fSvcf=$full_path."/".$name.".snp.vcf"; 	
		my $fImut=$full_path."/".$name.".indel.vcf.proteome-indel-mut.fasta";
		my $fIwt=$full_path."/".$name.".indel.vcf.proteome-indel-wt.fasta";
		my $fSmut=$full_path."/".$name.".snp.vcf.proteome-snv-mut.fasta";
		my $fSwt=$full_path."/".$name.".snp.vcf.proteome-snv-wt.fasta";
		my $fsup=$full_path."/".$name.".transcript.novel.support.rna.tsv"; 
		my $fneo=$full_path."/".$name.".neo.summary.min";
		my $fref=$full_path."/".$name.".transcript.fa.supporting.ref";
	
		if(-s $fneo)
		{
			foreach my $l (`cat $fneo`) 
			{
				my $ltr=$l; 
				chomp($ltr); 
				my @temp=split("\t",$ltr); 
				my $pos=$temp[0]."_".$temp[1]; 
				my $np=$temp[6];
				my $pep=$temp[7].";".$temp[8].";".$temp[9];
			    $neo{$name}{$np}{$pos}=$pep;  	
				#print $name,"\t",$np,"\t",$pos,"\n";	
			}	
		}

		if(-s $fref)
        {
            foreach my $l (`cat $fref`)
            {
                my $ltr=$l;
                chomp($ltr);
				if($ltr=~/^>NP/)
				{ 
                	my @temp=split("\t",$ltr);
            		my @temp2=split(/\-/,$temp[0]);
			    	my $pos=$temp2[2];
					my $chr=$temp2[1];
                	my $np=$temp2[0];
					$chr=~s/chr//g; 
					$pos=$chr."_".$pos; 
					$np=~s/\>//g; 
					#print $chr,"\t",$pos,"\t",$np,"\t",$temp[1],"\n";
					#<STDIN>;
                	$ref_sup{$name}{$np}{$pos}=$temp[1]; 
				}
			
            }
        }
	
		if(-s $fIvcf)
		{
			foreach my $l (`cat $fIvcf`)
			{
				chomp($l); 
				my @temp=split("\t",$l);

				if($temp[0]=~/^chr/) { $temp[0]=~s/^chr//g; }
				my $pos=$temp[0]."_".$temp[1]; 
				my $inf=$temp[2]."_".$temp[3]."_".$temp[4];  
				$mut{$name}{$pos}=$inf; 		
				#print $name,"\t",$pos,"\t",$inf,"\n";
			}		
		}
	
		if(-s $fSvcf)
        {
            foreach my $l (`cat $fSvcf`)
            {
                chomp($l);
                my @temp=split("\t",$l);
                my $pos=$temp[0]."_".$temp[1];
                my $inf=$temp[2]."_".$temp[3]."_".$temp[4];
                $mut{$name}{$pos}=$inf;
			  #print $name,"\t",$pos,"\t",$inf,"\n"; 
            }   
        }

		my $np;
		my $pos; 
		### parse the protein fasta file for wild type and mut###

		if(-s $fIwt)
		{
            foreach my $l (`cat $fIwt`)
            {
                chomp($l);
				if($l=~/^>NP/)
				{
				 my @temp=split /\s+/, $l; 
				 #print $l,"\n";
				 #print $temp[0],"\t",$temp[1],"\t",$temp[2],"\n";
				 $np=$temp[0]; 
				 $np=~s/^>//g;
				 $np=~s/-wt$//g;  
				}
				else 
				{ 
				$pepwt{$name}{$np}=$l;
			  	#print "pep indel wt","\t",$name,"\t",$np,"\n";
			 	}
            }
         }
		
		 if(-s $fSwt)
         {
            foreach my $l (`cat $fSwt`)
            {
                chomp($l);
                if($l=~/^>NP/)
                {
                 #my @temp=split /\s+/, $l;   
                 $np=$l; 
                 $np=~s/^>//g;  
                }
                else 
				{ 
					$pepwt{$name}{$np}=$l;
                 	#print "pep snv wt","\t",$name,"\t",$np,"\n";
					#<STDIN>;
				 }
            }
         }

		 if(-s $fSmut)
          {
            foreach my $l (`cat $fSmut`)
            {
                chomp($l);
                if($l=~/^>NP/)
                {
					my @temp=split(":",$l);
					$np=$temp[0];
					$np=~s/^>//g; 
					$np=~s/-VAR//g; 
					my @temp2=split("-",$temp[1]); 
					$temp2[1]+=1;
					$pos=$temp2[0]."_".$temp2[1]; 
					$pos=~s/^chr//g; 
				}
				else {  
					 $pepmut{$name}{$np}{$pos}=$l;
					 #print "pep snv mut","\t",$name,"\t",$np,"\t",$pos,"\n";
			         } 
            }

          }

         #if(-s $fIwt)
         #{
          #  foreach my $l (`cat $fIwt`)
          #  {
           #     chomp($l);
           #     if($l=~/^>NP/)
           #     {
                 #my @temp=split /\s+/, $l;
            #     my $np=$temp[0];
            #     $np=~s/^>//g;
            #     $np=~s/-wt$//g;
            #    }
             #   else { $pepwt{$name}{$np}=$l;
              #    		print "pep indel wt","\t",$name,"\t",$np,"\n";
              #   }
           # }
         	#}
		  if(-s $fImut)
          {
            foreach my $l (`cat $fImut`)
            {
                chomp($l);
                if($l=~/^>NP/)
                {
                    my @temp=split(/\s+/,$l);
                    $np=$temp[0];
					$np=~s/^>//g; 
                    $np=~s/-indel//g;
                    my @temp2=split("-",$temp[4]);
                    $temp2[1]+=1;
                    $pos=$temp2[0]."_".$temp2[1];
                    $pos=~s/^\(INS:chr//g;
                    $pos=~s/^\(DEL:chr//g;
                }
                else {  $pepmut{$name}{$np}{$pos}=$l; 
						 #print "pep indel mut","\t",$name,"\t",$np,"\t",$pos,"\n";
					}
            }
		   
          }

		  if(-s $fsup) 
		  {
			foreach my $l (`cat $fsup`)
            {
			chomp($l);
            if($l=~/^NP\_/)
             {
			    my @temp=split("\t",$l);
			  	my @temp2=split("-",$temp[0]);
				if($temp[0]=~/indel/)
				{
					$np=$temp2[0];
                	$temp2[3]+=1; 
                	$pos=$temp2[2]."_".$temp2[3];
                	$pos=~s/^chr//g;
				}
				else 
				{
					$np=$temp2[0];
					$temp2[2]+=1; 
					$pos=$temp2[1]."_".$temp2[2];
			    	$pos=~s/^chr//g; 
				}
				$suprna{$name}{$np}{$pos}=$temp[1]; 
				#print "supporting","\t",$name,"\t",$np,"\t",$pos,"\n";
			 }
		  }
		 }	
		}
	}
} ## subroutine
