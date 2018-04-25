## 12-09-2016 ##
#!/usr/bin/perl

use strict;

(my $usage = <<OUT) =~ s/\t+//g;
This script will extract reads supporting reference genome.
OUT

die $usage unless @ARGV == 3;


my($f_bam_in,$db_refseq,$f_fa_in)=@ARGV;

my %junc_bed=();
my %mut_pos=();
my $f_out=$f_fa_in.".supporting.ref";
#my $f_out_q20=$f_fa_in.".supporting.ref.q20"; 

#open(IN,"<$f_map_in");
open(DB,"<$db_refseq");
open(FA,"<$f_fa_in"); 

open(OUT,">$f_out");
### refseq is 0-based ###
my $n_exon; 
my $l_exon;
my $p_exon; 

while(<DB>)
	{

		my $line=$_; 
		chomp($line); 

		my @temp=split("\t",$line); 
		my $left=0;
		my $right=0; 
		my $pos=$temp[6];
		my $chr=$temp[0];
		$n_exon=$temp[9];
		$l_exon=$temp[10];
		$p_exon=$temp[11];

		my @temp_l=split(",",$l_exon); 
		my @temp_p=split(",",$p_exon); 
		
		for(my $i=0;$i<$n_exon-1;$i++)
		{
			$right=$pos+$temp_p[$i+1]+1;
			$left=$pos+$temp_p[$i]+$temp_l[$i];
			#if($temp[3] eq "NP_001120800")
			#{ print $chr,"\t",$left,"\t",$right,"\n"; } 
		 	$junc_bed{$temp[3]}{$chr}{$left}{$right}=1; 	
		}			
			
	}
close DB; 
	 
while(<FA>)
	{
		my $line=$_; 
		chomp($line); 
		if($line=~/^\>/) { 
		my @temp=split(/\-/,$line); 
		my $chr;
		my $pos; 
		$temp[0]=~s/^>//g; 
		if($line=~/indel/)
			{
			$chr=$temp[2];
			$pos=$temp[3];
			if(!($pos=~/chr/))
			{
				#print "indel","\t",$temp[0],"\t",$chr,"\t",$pos+1,"\n";
				#<STDIN>;
				$mut_pos{$temp[0]}{$chr}{$pos+1}=1;	
			}
			}
		else 
			{
			$chr=$temp[1];
            $pos=$temp[2];
            if(!($pos=~/chr/))
           	{
				#print "snv","\t",$temp[0],"\t",$chr,"\t",$pos+1,"\n";
				#<STDIN>;
                $mut_pos{$temp[0]}{$chr}{$pos+1}=1;
            }					
			}
			}
	}

close FA;
my %numa; 
my %chra;
my $nmax;
my $chrmax;

foreach my $line (`samtools view $f_bam_in | head -n 100000`)
    {
        chomp($line);
		#print $line,"\n";
        if($line=~/^\@SQ/) { next; }
        my @temp=split("\t",$line);
        my $cigar=$temp[5];
		my $chr=$temp[2];		
        if($cigar=~/(\d+)M/) {
        $numa{$1}++;
		$chra{$chr}++;	
		 }
    
	}
close IN;

foreach my $n (sort { $numa{$b} <=> $numa{$a} } keys %numa) 	
	{
		$nmax=$n;
		last;  	
	}

foreach my $chr (sort { $chra{$b} <=> $chra{$a} } keys %chra)
    {
        $chrmax=$chr;
        last;
    }

#print $nmax,"\t",$chrmax,"\n"; 

foreach my $np (sort keys %mut_pos)
	{
	foreach my $chr (sort keys %{$mut_pos{$np}})
		{
			foreach my $pos (sort keys %{$mut_pos{$np}{$chr}})
			{
				#my $left_pos=$pos-20;
         		#my $right_pos=$pos+20;
				my $chr_=$chr; 	
				if(!($chrmax=~/chr/)) { $chr_=~s/chr//g; }
				my $left_pos=$pos-20;
         		my $right_pos=$pos+20;
         		my $chr_pos=$chr_.":".$left_pos."-".$right_pos;  
				#print $chr_pos,"\t",$pos,"\n";
				#<STDIN>;
				my $f_abs_bam=`readlink -f $f_bam_in`;
				#print $f_abs_bam,"\n";
				chomp($f_abs_bam);
				my $com=`samtools view $f_abs_bam \"$chr_pos\"`;
				my @temp=split("\n",$com);
				my %seq=();
				my %quality=();
				foreach my $t (@temp)
            	{
					#print $t,"\n";
					#<STDIN>;
					my @temp2=split("\t",$t); 
					my $cigar=$temp2[5];	
					my $nm=$temp2[12];
					my $flag=$temp2[1];
					my $leftmostpos=$temp2[3];	 
        			my $mapq=$temp2[4];
					#print $chr,"\t",$pos,"\n";
					#print $t,"\n";
					#print $cigar,"\t",$flag,"\t",$mapq,"\n";
					#<STDIN>;
					if($cigar=~/^$nmax\M/ && $t=~/NM:i:0/ && $mapq>=20) 
					{
						if($pos>=$leftmostpos && $pos<=$leftmostpos+$nmax-1)
						{
						my $id=$temp2[0];
						if(($flag & 0x40) || $id=~/\/1$/) 
						{
						my $rid=$id; 
						$rid=~s/\/1$//g; 
						#my $rid=$id."#1"; 
						$rid=$rid."#1"; 
						$quality{$rid}=$temp2[10]; 
						$seq{$rid}=$temp2[9];
						}
						if($flag & 0x80 || $id=~/\/2$/)
						{
						my $rid=$id;
                        $rid=~s/\/2$//g;
						$rid=$rid."#2";
						$quality{$rid}=$temp2[10];
            			$seq{$rid}=$temp2[9];			
        				}
						}
					}

					elsif($cigar=~/(\d+)M(\d+)N(\d+)M/)
					{
				
				#print $np,"\t",$1,"\t",$2,"\t",$3,"\n"; 
				#<STDIN>;
					my $n1;
					my $n2;
					my $n3;
					$n1=$1; 
					$n2=$2;
					$n3=$3;		
				#print $np,"\t",$n1,"\t",$n2,"\t",$n3,"\n";	
			if($t=~/NM:i:0/ && $mapq>=20 && (defined $junc_bed{$np}{$chr}{$leftmostpos+$n1-1}{$leftmostpos+$n1+$n2}) && $n1+$n3==$nmax && (($pos>=$leftmostpos && $pos<=$leftmostpos+$n1-1) || ($pos>=$leftmostpos+$n1+$n2 && $pos<=$leftmostpos+$n1+$n2+$n3-1)))
						{
                        my $id=$temp2[0];
						#print $id,"\n";
						my $rid;
                        if(($flag & 0x40) || $flag=~/\/1$/)
                        {
                        $rid=$id;
                        $rid=~s/\/1$//g;
                        $rid=$rid."#1";
                        $quality{$rid}=$temp2[10];
                        $seq{$rid}=$temp2[9];
                        }
                        if($flag & 0x80 || $flag=~/\/2$/)
                        {
                         $rid=$id;
						$rid=~s/\/2$//g;
                        $rid=$rid."#2";
                        $quality{$rid}=$temp2[10];
                        $seq{$rid}=$temp2[9];
                        }
						#print $rid,"\t",$quality{$rid},"\t",$seq{$rid},"\n";
						#<STDIN>;
						}	
					}
				} ## foreach t
				my $n_read=keys %seq; 
				print OUT ">",$np,"-",$chr,"-",$pos,"\t",$n_read,"\n";
				foreach my $rd (sort keys %seq)
				{
					print OUT $rd,"\t",$seq{$rd},"\t",$quality{$rd},"\n";
				} 
				#close OUT;
			} ## pos
		} ## chr
	} ## np

close OUT;
#foreach my $id (sort keys %seq) 
#	{
#		print OUT "\@",$id,"\n"; 
#		print OUT $seq{$id},"\n";
#		print OUT "+","\n";
#		print OUT $quality{$id},"\n"; 
#	}

