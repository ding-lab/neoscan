#!/usr/bin/perl

## generate mutant and wild-type protein sequence fir indel ##
## 12-01-2016 ##

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will get the protein sequence for both wide-type and mutants 


OUT

die $usage unless @ARGV == 3;

my $error=0;
#my $filename="";
#my $filename_vcf_germline="";
#my $filename_vcf_somatic="";
#my $dir=".";

my($filename,$dir,$filename_vcf_somatic)=@ARGV;

#if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }
#if ($ARGV[1]=~/\w/) { $dir=$ARGV[1];} else { $dir="."; }
#if ($ARGV[2]=~/\w/) { $filename_vcf_somatic=$ARGV[2];} else { $filename_vcf_somatic=""; }

my @dir_vcf=split(/\//,$filename_vcf_somatic);

#if ($ARGV[3]=~/\w/) { $filename_vcf_somatic=$ARGV[2];} else { $filename_vcf_somatic=""; }

my %mapping = ("TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
				"CTT"=>"L","CTC"=>"L","CTA"=>"L","CTG"=>"L",
				"ATT"=>"I","ATC"=>"I","ATA"=>"I","ATG"=>"M",
				"GTT"=>"V","GTC"=>"V","GTA"=>"V","GTG"=>"V",
				
				"TCT"=>"S","TCC"=>"S","TCA"=>"S","TCG"=>"S",
				"CCT"=>"P","CCC"=>"P","CCA"=>"P","CCG"=>"P",
				"ACT"=>"T","ACC"=>"T","ACA"=>"T","ACG"=>"T",
				"GCT"=>"A","GCC"=>"A","GCA"=>"A","GCG"=>"A",
				
				"TAT"=>"Y","TAC"=>"Y","TAA"=>"*","TAG"=>"*",
				"CAT"=>"H","CAC"=>"H","CAA"=>"Q","CAG"=>"Q",
				"AAT"=>"N","AAC"=>"N","AAA"=>"K","AAG"=>"K",
				"GAT"=>"D","GAC"=>"D","GAA"=>"E","GAG"=>"E",
				
				"TGT"=>"C","TGC"=>"C","TGA"=>"*","TGG"=>"W",
				"CGT"=>"R","CGC"=>"R","CGA"=>"R","CGG"=>"R",
				"AGT"=>"S","AGC"=>"S","AGA"=>"R","AGG"=>"R",
				"GGT"=>"G","GGC"=>"G","GGA"=>"G","GGG"=>"G");
				
if ($error==0)
{

  my $filename_=$filename;

    #$filename_="pr;

    open (OUT,">$filename_vcf_somatic.fa-indel-wt.fasta");
    open (OUT_MOD,">$filename_vcf_somatic.fa-indel-mut.fasta");
    #open (OUT_FS,">$filename-frame-shift.fasta");
    open (LOG,">$filename_vcf_somatic.fa-indel-wt.log");
    open (LOG_MOD,">$filename_vcf_somatic.fa-indel-mut.log");
    open (STAT,">$filename_vcf_somatic.fa-indel.stat");


	my $proteins_count=0;
	my $proteins_modified=0;
	my $count_stop_removed=0;
	my $count_stop_introduced=0;
	my $count_stop_removed_somatic=0;
	my $count_stop_introduced_somatic=0;
	my $protein_modifications_count=0;
	my @protein_modifications_distr=();
	my $protein_modifications_count_somatic=0;
	my @protein_modifications_distr_somatic=();
	my $count_variant_in_exon=0;
	my $count_variant_in_exon_somatic=0;
	my $count_variant_in_exon_old_error=0;
	my $count_variant_in_exon_old_error_somatic=0;
	my $count_variant_in_exon_nonsyn=0;
	my $count_variant_in_exon_nonsyn_somatic=0;
	my $line="";
	my %chr=();
	my %bed=();
	my %seq=();

	if (open (IN,"$filename"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
			{
				my $chr=$1;
				my $name=$4;
				$bed{$name}=$line;
				$chr{$chr}.="#$name#";
			} else { print LOG qq!Error parsing: $line!; }
		}
		close(IN);
	}
	my %descriptions=();

	#if (open (IN,"$filename_-descriptions.txt"))
	#{
	#	while ($line=<IN>)
	#	{
	#		chomp($line);
	#		if ($line=~/^([^\t]+)\t([^\t]+)/)
	#		{
	#			$descriptions{$1}=$2;
	#		}
	#	}
	#	close(IN);
	#}

	my %vcf_old=();
	my %vcf_new=();
	my %vcf_type=();
	my %vcf_anno=();
	my %vcf_gene=(); 
	my $germline_vcf=0;
	my $germline_only=0;
	my $somatic_only=0;
	my $both_vcf=0;
	my $both_vcf_differ=0;
	#print "open $filename_vcf_somatic\n"; <STDIN>;

    if(!-e $filename_vcf_somatic)
    {
        print STAT "no vcf file\n";
        last;
    }

    open (IN,"<$filename_vcf_somatic");

	while ($line=<IN>)
		{
			chomp($line);
			print $line,"\n";
              my @temp=split("\t",$line);
                my $chr="chr".$temp[0];
                my $pos=$temp[1];
                my $old=$temp[2];
                my $new=$temp[3];
                my $gene=$temp[4];
                my $anno=$temp[5];
                my $type=$temp[6];		
				$pos--;
				$new=~s/\,.*$//; #####
				print "pos=$pos","\t","old=$old", "\t", "new=$new","\t","anno=$anno","\t","type=$type","\n"; 
				#<STDIN>;
				#if ($vcf_old{"$chr#$pos"}=~/\w/)
				#{
				#	if ($vcf_new{"$chr#$pos"}!~/^$new$/)
				#	{
						#print LOG_MOD qq!$chr $pos: germline:$vcf_old{"$chr#$pos"}->$vcf_new{"$chr#$pos"} somatic:$old->$new\n!;
						$vcf_old{"$chr#$pos"}=$old;
						$vcf_new{"$chr#$pos"}=$new;
						if($type eq "SOMATIC") 
						{ $vcf_type{"$chr#$pos"}="S";}
						else { $vcf_type{"$chr#$pos"}="G"; }

						$vcf_anno{"$chr#$pos"}=$anno;
						#print OUT_ONLY qq!$line\n!;
						$somatic_only++;
						$vcf_gene{"$chr#$pos"}=$gene;
		}
		close(IN);
		#close(OUT_ONLY);
	
	foreach my $chr (sort keys %chr)
	{
		#$chr="chr9"; 
		print qq!$chr\n!;
        my $chr2=$chr;
        $chr2=~s/chr//g;
        my $f_in=$dir."/Homo_sapiens.GRCh38.dna.chromosome.".$chr2.".fa";
        if (open (IN,"$f_in"))
#		if (open (IN,"$dir/$chr.fa"))
		{
			print qq!opened $chr\n!;
			my $sequence="";
			$line=<IN>;
			chomp($line);
			my $chr_=$chr;
			$chr_=~s/^chr//i;
			if ($line=~/^>$chr\s*/ or $line=~/^>$chr_\s*/)
			{
				while ($line=<IN>)
				{
					chomp($line);
					if ($line=~/^>/)
					{
						print LOG qq!Error: > not expected: $line\n!;
					}
					else
					{
						$line=~s/\s+//g;
						if ($line!~/^[atcgATCGnN]+$/)
						{
							print LOG qq!Error: unexpected character: $line\n!;
						}
						else
						{
							$sequence .= "\U$line";
						}
					}
				}
				my $temp=$chr{$chr};
				my %seq_chr_pos; 

				while ($temp=~s/^#([^#]+)#//)
				{
					my $name=$1;
					my $modified=0;
					#print LOG qq!\n$name: $bed{$name}\n!;
					#print qq!\n$name: $bed{$name}\n!;
					if ($bed{$name}=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)$/)
					{
						my $start=$2;
						my $end=$3;
						my $strand=$6;
						my $num=$10;
						my $segment_lengths="$11,";
						my $segment_starts="$12,";
						my $segment_lengths_extended="$11,";
						my $segment_starts_extended="$12,";

						#if ($strand=~/\+/) { if ($segment_lengths_extended=~s/([0-9]+)\,$//) { my $temp=$1; $temp+=100; $segment_lengths_extended.="$temp,"; } }
						#else 
						#{
						#	if ($segment_lengths_extended=~s/^([0-9]+)\,//) { my $temp=$1; $temp+=100; $segment_lengths_extended="$temp,$segment_lengths_extended,"; }
						#	if ($segment_starts_extended=~s/^([0-9]+)\,//) { my $temp=$1; $temp-=100; $segment_starts_extended="$temp,$segment_starts_extended,"; }
						#}

						#print LOG qq!$segment_lengths $segment_lengths_extended   $segment_starts $segment_starts_extended\n!;

						my $seq_original="";
						my $segment_starts_=$segment_starts;
						my $segment_lengths_=$segment_lengths;
						my $accu_c=0; 
						while ($segment_starts_=~s/^([0-9\-]+)\,//)
						{
							my $segment_start=$1;
							if ($segment_lengths_=~s/^([0-9]+)\,//)
							{
								my $segment_length=$1;
									
								my $seq_=substr $sequence,$start+$segment_start,$segment_length;
								for(my $pos=1;$pos<=$segment_length;$pos++) { $seq_chr_pos{$accu_c}=$start+$segment_start+$pos-1; $accu_c++;}
								$seq_original.=$seq_; 
								#print LOG qq!$start+$segment_start,$segment_length: $seq_\n!;
							} else { print LOG qq!Error parsing $bed{$name}\n!; }
						}
						
						my $seq=$seq_original;
						my $seq_original_rc= reverse $seq_original;
                                                $seq_original_rc=~tr/ATCG/TAGC/;

						my $description_="";
						my %variants=();
						my %variants_=();
						my %variants2=();
						my %variants2_=();
						my $variant_loc=""; 
						my $seqment_lengths_sum=0;
						my $seqment_count=0;
						my $num_indel=0; 
						my $segment_start_first=0;
						#my $find_indel; 
						$segment_start_first=0;

						#$segment_starts_=$segment_starts_extended;
						#$segment_lengths_=$segment_lengths_extended;
						#make modification for the reference genome based on the indel and snps

						for(my $i=length($seq_original)-1;$i>=0;$i--) 
							{
							  my $var_pos=$seq_chr_pos{$i};							  
							  #print $chr,"\t",$var_pos,"\n"; 
							  if(defined $vcf_old{"$chr#$var_pos"}) { 
								#print "finding indel\n"; 
								print $name,"\n";		
								print $var_pos,"\n";
								print $vcf_old{"$chr#$var_pos"},"\n"; 
								print $vcf_new{"$chr#$var_pos"},"\n";
								#<STDIN>; 
								#print substr($seq_original,$i-1,3),"->",$vcf_old{"$chr#$var_pos"},"->",$vcf_new{"$chr#$var_pos"},"\n";
								#<STDIN>;
								my $leno; 
								my $lenn;
								if($vcf_old{"$chr#$var_pos"} eq "-") { $leno=0; }
								else { $leno=length($vcf_old{"$chr#$var_pos"}); }	
								if($vcf_new{"$chr#$var_pos"} eq "-") { $lenn=0; }
                                                                else { $lenn=length($vcf_new{"$chr#$var_pos"}); } 
								$num_indel++; 	
								my $inframe;
								if(abs($leno-$lenn)%3==0) { $inframe=1; }
								else { $inframe=0; }												
								if($strand=~/\-/) { 

										    if($vcf_new{"$chr#$var_pos"} eq "-") 
												{ 
												my $rev_pos=length($seq_original)-$i-$leno; 
												my $left_3=($rev_pos+1)%3;	
												my $int_3s=int (($rev_pos+1)/3)+1;
												my $int_3e=int (($rev_pos+$leno+1)/3)+1; 
												if($inframe==1) { $description_.=$chr."-".$var_pos;} 
												else { $description_.=$chr."-".$var_pos; }
												}

										     else {
											        my $rev_pos=length($seq_original)-$i-$leno;
                                                                                                my $left_3=($rev_pos+1)%3;
                                                                                                my $int_3s=int (($rev_pos+1)/3)+1;
                                                                                                my $int_3e=int (($rev_pos+$lenn+1)/3)+1;
                                                                                                if($inframe==1) { $description_.=$chr."-".$var_pos; }
                                                                                                else { $description_.=$chr."-".$var_pos; }
											  }

											}
										
								else {
									if ($vcf_new{"$chr#$var_pos"} eq "-")
                                                                                                {
                                                                                                my $left_3=($i+1)%3;
												#print "i=$i\n"; 
												#print substr($seq_original,$i,40),"\n"; 
                                                                                                my $int_3s=int (($i+1)/3)+1;
												#print "int_3s=$int_3s\n";
												#<STDIN>; 
                                                                                                my $int_3e=int (($i+$leno+1)/3)+1;
                                                                                                if($inframe==1) { $description_.=$chr."-".$var_pos;}
                                                                                                else { $description_.=$chr."-".$var_pos; }
                                                                                                }
                                                                         else {
                                                                                                my $left_3=($i+1)%3;
                                                                                                my $int_3s=int (($i+1)/3)+1;
                                                                                                my $int_3e=int (($i+$lenn+1)/3)+1;
                                                                                                if($inframe==1) { $description_.=$chr."-".$var_pos;}
                                                                                                else { $description_.=$chr."-".$var_pos; }
                                                                                          }
									}

								#print $description_,"\n"; 
								#if($leno != $lenn) { $num_fs++; }		

								my $seql=substr($seq,0,$i); 
								print $seql,"\n";
								my $seqr=substr($seq,$i+$leno,length($seq)-length($seql)-$leno);

								if($lenn>=1) 
								{
								 $seq=$seql.$vcf_new{"$chr#$var_pos"}.$seqr; 
								} 
								else { $seq=$seql.$seqr; }
								}  	
							}

						my $name_=$name; $name_=~s/\-[^\-]+$//;
						#print "orignal sequence\n"; 
						#print $seq_original,"\n"; 
						#print "modified sequence\n"; 
						#print $seq,"\n";  
						#my $gene=$name; $gene=~s/^[^\-]+\-//;

						my $protein="";
						my $protein_original="";
					
						my $lenb_added=length($seq_original);


                        if($strand=~/\-/) {
                        my $seq_added=substr($sequence,$start-10000,10000);
                        $seq_original=$seq_added.$seq_original;
                        }

                        else {
                        my $seq_added=substr($sequence,$end,10000);
                        $seq_original=$seq_original.$seq_added;
                        }


                        if($strand=~/\-/) {
                        my $seq_added=substr($sequence,$start-10000,10000);
                        $seq=$seq_added.$seq;
                        }
                        else {
                        my $seq_added=substr($sequence,$end,10000);
                        $seq=$seq.$seq_added;
                        }
	
						if ($strand=~/\-/)
						{
							my $seq_ = reverse $seq;
							$seq=$seq_;
							$seq=~tr/ATCG/TAGC/;
							$seq_ = reverse $seq_original;
							$seq_original=$seq_;
							$seq_original=~tr/ATCG/TAGC/;
						}
					
						my $stop_found=0;

                        for(my $n=0;$n<length($seq_original) and $stop_found==0;$n=$n+3)
                        {
                            my $triplet = substr($seq_original, $n, 3);
                            if (length($triplet)==3)
                            {
                                if ($mapping{$triplet}=~/\*/)
                                {
                                    $stop_found=1;
                                    if($n+3>$lenb_added+1)
                                    {
                                        print LOG qq!Error: Stop codon not found within the bed file:$name \n$protein_original\n!;
                                        $protein_original="";
                                    }
                                }
                                if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
                                $protein_original.=$triplet;
                            }
                        }
	
						#for(my $n=0;$n<length($seq_original);$n=$n+3)
						#{
						#	my $triplet = substr($seq_original, $n, 3);
						#	if (length($triplet)==3)
						#	{
						#		if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
						#		$protein_original.=$triplet; 
						#	}
						#}

						#$modified=0;
						#my $modified_somatic=0;
						$stop_found=0;
						my $triplet_count=0;
						#my $frame_shift=0;
						#if ($description_=~/\w/) { $description_.="; "; }
						for(my $n=0;$n<length($seq) and $stop_found==0;$n=$n+3)
						{
							my $n_=$n+2;
							#my $triplet_count_=$triplet_count+1;
							my $triplet = substr($seq, $n, 3);
							if (length($triplet)==3)
							{

								if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
								#print LOG_MOD qq!$name $triplet_type $n-$n_: $triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet}\n!;
								if ($mapping{$triplet}=~/\*/)
								{
									$stop_found=1;
								}
									$protein.=$triplet;
							}

							$triplet_count++;
						}
						#$proteins_count++;
						#$protein_modifications_count+=$modified;
						#$protein_modifications_distr[$modified]++;
						#$protein_modifications_count_somatic+=$modified_somatic;
						#$protein_modifications_distr_somatic[$modified_somatic]++;
						#$protein_original=~s/\*$//;
						#print "protein original\n"; 
						#print $protein_original,"\n"; 
						#<STDIN>;

						if ($protein_original=~/^([^\*]+)\*.*$/)
							{
							 print LOG qq!Error: No Stop codon found in middle of sequence:$name \n$protein_original\n!;
							 $protein_original="";
							}

						#if (length($protein_original)>6 and $protein_original!~/^\*/)
						#	{
						#	print OUT qq!>$name (MAP:$chr:$start$strand $segment_lengths $segment_starts)\n$protein_original\n!;
						#	}
						
						 #$protein=~s/\*$//;
						 #$protein=~s/^([^\*]+)\*.*$//;

					         if(($protein ne $protein_original) && ($num_indel==1))
                                                	{

							  #$proteins_modified++;
							 # print "original protein\n"; 
							 # print $protein_original,"\n"; 
							 # print "modified protein\n"; 
							 # print $protein,"\n"; 

							  if (length($protein)>6 && $protein!~/^\*/ && length($protein_original)>6)
							   {
								print OUT qq!>$name-wt\n$protein_original\n!;
								#$description_.=")";
								#print LOG_MOD qq!>$name-indel $descriptions{$name} (MAP:$chr:$start$strand $segment_lengths $segment_starts) $description_\n$protein\n!;
								print OUT_MOD qq!>$name-indel-$description_\n$protein\n!;
							   }

						         }

					} else { print LOG qq!Error parsing $bed{$name}\n!; }
				}
			} else { print LOG qq!Error in name $chr: $line\n!; }
			
			close(IN);
		}
	}
	print STAT qq!
		done; 
		!;
	print STAT qq!\nnumber of modifications\tnumber of proteins\tnumber of proteins (somatic variants)\n!;
	#for (my $k=0;$k<=100;$k++) { print STAT qq!$k\t$protein_modifications_distr[$k]\t$protein_modifications_distr_somatic[$k]\n!; }
	close(OUT);
	close(OUT_MOD);
	#close(OUT_FS);
	close(LOG);
	close(LOG_MOD);
	close(STAT);
}
