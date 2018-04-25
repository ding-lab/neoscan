#
use strict;

#!/usr/bin/perl

## generate mutant and wild-type protein sequence for snv ##
## 12-01-2016 ##

use strict;

(my $usage = <<OUT) =~ s/\t+//g;

This script will get the protein sequence for both wide-type and mutants 


OUT

my $error=0;

die $usage unless @ARGV == 3;


my($filename,$dir,$filename_vcf_somatic)=@ARGV;

#if ($ARGV[0]=~/\w/) { $filename=$ARGV[0];} else { $error=1; }
#if ($ARGV[1]=~/\w/) { $dir=$ARGV[1];} else { $dir="."; }
#if ($ARGV[2]=~/\w/) { $filename_vcf_germline=$ARGV[2];} else { $filename_vcf_germline=""; }
#if ($ARGV[3]=~/\w/) { $filename_vcf_somatic=$ARGV[3];} else { $filename_vcf_somatic=""; }


my %mapping = (	"TTT"=>"F","TTC"=>"F","TTA"=>"L","TTG"=>"L",
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
	$filename_=~s/\.bed$//;
	open (OUT,">$filename_vcf_somatic.proteome-snv-wt.fasta");
    open (OUT_MOD,">$filename_vcf_somatic.proteome-snv-mut.fasta");
    open (OUT_FS,">$filename_vcf_somatic.proteome-snv-fs.fasta");
    open (LOG,">$filename_vcf_somatic.proteome-snv-wt.log");
    open (LOG_MOD,">$filename_vcf_somatic.proteome-snv-mut.log");
    open (STAT,">$filename_vcf_somatic.proteome-snv-mod.stat");
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

	if (open (IN,"$filename_-descriptions.txt"))
	{
		while ($line=<IN>)
		{
			chomp($line);
			if ($line=~/^([^\t]+)\t([^\t]+)/)
			{
				$descriptions{$1}=$2;
			}
		}
		close(IN);
	}

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
#	if (open (IN,"$filename_vcf_germline"))
#	{
#		while ($line=<IN>)
#		{
#			chomp($line);
#			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
#			{
#				my $chr="chr$1";
#				my $pos=$2;
#				my $id=$3;
#				my $old=$4;
#				my $new=$5;
#				my $qaul=$6;
#				$pos--;
#				$new=~s/\,.*$//; #####
#				$vcf_old{"$chr#$pos"}=$old;
#				$vcf_new{"$chr#$pos"}=$new;
#				$vcf_type{"$chr#$pos"}="G";
#				$germline_vcf++;
				#print qq!$chr#$pos#$new\n!;
#			} else { print LOG_MOD qq!Error parsing: $line!; }
#		}
#		close(IN);
#	}

#	if (open (IN,"$filename_vcf_somatic"))
#	{
		#open (OUT_ONLY,">$filename_vcf_somatic-somatic_only.vcf");
#		while ($line=<IN>)
#		{
#			chomp($line);
#			if ($line=~/^([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)\t([^\t]*)/)
#			{
#				my $chr="chr$1";
#				my $pos=$2;
#				my $id=$3;
#				my $old=$4;
#				my $new=$5;
#				my $qaul=$6;
#				$pos--;
#				$new=~s/\,.*$//; #####
#				if ($vcf_old{"$chr#$pos"}=~/\w/)
#				{
#					if ($vcf_new{"$chr#$pos"}!~/^$new$/)
#					{
						#print LOG_MOD qq!$chr $pos: germline:$vcf_old{"$chr#$pos"}->$vcf_new{"$chr#$pos"} somatic:$old->$new\n!;
#						$vcf_old{"$chr#$pos"}=$old;
#						$vcf_new{"$chr#$pos"}=$new;
#						$vcf_type{"$chr#$pos"}="S";
						#print OUT_ONLY qq!$line\n!;
#						$somatic_only++;
#						$both_vcf_differ++;
#					}
#					$both_vcf++;
#				}
#				else
#				{
#					$vcf_old{"$chr#$pos"}=$old;
#					$vcf_new{"$chr#$pos"}=$new;
#					$vcf_type{"$chr#$pos"}="S";
					#print OUT_ONLY qq!$line\n!;
#					$somatic_only++;
#				}
				#print qq!$chr#$pos#$new\n!;
#			} else { print LOG_MOD qq!Error parsing: $line!; }
#		}
#		close(IN);
		#close(OUT_ONLY);
#	}
#	$germline_only=$germline_vcf-$both_vcf_differ;
#	print qq!Somatic: $somatic_only\nGermline: $germline_only\nIn both VCF: $both_vcf (differ: $both_vcf_differ)\n!;

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
                #   if ($vcf_new{"$chr#$pos"}!~/^$new$/)
                #   {
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
	
	foreach my $chr (sort keys %chr)
	{
		print qq!$chr\n!;
		if (open (IN,"$dir/$chr.fa"))
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
				while ($temp=~s/^#([^#]+)#//)
				{
					my $name=$1;
					my $modified=0;
					#print LOG qq!\n$name: $bed{$name}\n!;
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
						if ($strand=~/\+/) { if ($segment_lengths_extended=~s/([0-9]+)\,$//) { my $temp=$1; $temp+=100; $segment_lengths_extended.="$temp,"; } }
						else 
						{
							if ($segment_lengths_extended=~s/^([0-9]+)\,//) { my $temp=$1; $temp+=100; $segment_lengths_extended="$temp,$segment_lengths_extended,"; }
							if ($segment_starts_extended=~s/^([0-9]+)\,//) { my $temp=$1; $temp-=100; $segment_starts_extended="$temp,$segment_starts_extended,"; }
						}
					#	print LOG qq!$segment_lengths $segment_lengths_extended   $segment_starts $segment_starts_extended\n!;
						my $seq_original="";
						my $segment_starts_=$segment_starts;
						my $segment_lengths_=$segment_lengths;
						while ($segment_starts_=~s/^([0-9\-]+)\,//)
						{
							my $segment_start=$1;
							if ($segment_lengths_=~s/^([0-9]+)\,//)
							{
								my $segment_length=$1;
								my $seq_=substr $sequence,$start+$segment_start,$segment_length;
								$seq_original.=$seq_; 
							#	print LOG qq!$start+$segment_start,$segment_length: $seq_\n!;
							} else { print LOG qq!Error parsing $bed{$name}\n!; }
						}
						my $seq="";
						my $description_="(VAR:";
						my %variants=();
						my %variants_=();
						my %variants2=();
						my %variants2_=();
						my $variant_loc=""; 
						my $seqment_lengths_sum=0;
						my $seqment_count=0;
						my $segment_start_first=0;
						$segment_start_first=0;
						#$segment_starts_=$segment_starts_extended;
						#$segment_lengths_=$segment_lengths_extended;
						$segment_starts_=$segment_starts;
						$segment_lengths_=$segment_lengths;
						if ($segment_starts_=~/^([0-9\-]+)\,/) { $segment_start_first=$1; }
						while ($segment_starts_=~s/^([0-9\-]+)\,//)
						{
							my $segment_start=$1;
							if ($segment_lengths_=~s/^([0-9]+)\,//)
							{
								my $segment_length=$1;
								$seqment_lengths_sum+=$segment_length;
								#print qq!$name $seqment_count. $seqment_lengths_sum\n!;
								my $seq_=substr $sequence,$start+$segment_start,$segment_length;
								my $seq__="";
								my $changed=0;
								for(my $i=0,my $j=$start+$segment_start;$i<$segment_length;$i++,$j++)
								{
									my $temp="";
									if ($vcf_new{"$chr#$j"}=~/\w/ and $vcf_new{"$chr#$j"}!~/^$vcf_old{"$chr#$j"}$/ and length($vcf_new{"$chr#$j"})==length($vcf_old{"$chr#$j"}))
									{
										my $n=substr $seq_,$i,length($vcf_old{"$chr#$j"});
										print qq!$chr#$j: $n $vcf_new{"$chr#$j"}\n!;
										if ($n!~/^$vcf_old{"$chr#$j"}$/) 
										{ 
											$temp=qq! (Warning: $n\!\=$vcf_old{"$chr#$j"}) - Ignored!; 
											my $i_=$i+$segment_start_first+length($seq);
						#					print LOG_MOD qq!Variant $chr $j: $n,$vcf_old{"$chr#$j"}$temp->$vcf_new{"$chr#$j"} $i_\n!;
											$count_variant_in_exon_old_error++;
											if ($vcf_type{"$chr#$j"}=~/S/) 
											{
												$count_variant_in_exon_old_error_somatic++;
											}
										}
										else
										{
											my $i_=$i+$segment_start_first+length($seq);
											if ($i_>=0)
											{
												$modified++;
												$seq__.=$vcf_new{"$chr#$j"};
												$changed=length($vcf_old{"$chr#$j"});
					#							print LOG_MOD qq!Variant $chr $j: $vcf_type{"$chr#$j"} $n,$vcf_old{"$chr#$j"}->$vcf_new{"$chr#$j"} $i_\n!;
												#$description_.=qq!$chr-$j($i_)-$vcf_type{"$chr#$j"}:$vcf_old{"$chr#$j"}->$vcf_new{"$chr#$j"}, !;
												$count_variant_in_exon++;
												if ($vcf_type{"$chr#$j"}=~/S/) 
												{ 
													$count_variant_in_exon_somatic++; 
												}
												$variants{$i_}=qq!$vcf_old{"$chr#$j"}#$vcf_new{"$chr#$j"}#$vcf_type{"$chr#$j"}!;
												$variants2{$i_}=qq!$chr-$j!; 
				#								print LOG_MOD qq!-$i_. $variants{$i_} $variants2{$i_}\n!;
											}
											else
											{
												print LOG_MOD qq!Ignored: Variant $chr $j: $vcf_type{"$chr#$j"} $n,$vcf_old{"$chr#$j"}->$vcf_new{"$chr#$j"} $i_\n!;
											}
										}
									}
									if ($changed<=0) { $seq__.=substr $seq_,$i,1; $changed=0; } else { $changed--; }
								}
			#					print LOG_MOD qq!$start+$segment_start,$segment_length: $seq__\n!;
								$seq.=$seq__; 
								$seqment_count++;
							} else { print LOG_MOD qq!Error parsing $bed{$name}\n!; }
						}
						my $name_=$name; $name_=~s/\-[^\-]+$//;
						#my $gene=$name; $gene=~s/^[^\-]+\-//;
						my $protein="";
						my $protein_original="";
						if ($strand=~/\-/)
						{
							my $seq_ = reverse $seq;
							$seq=$seq_;
							$seq=~tr/ATCG/TAGC/;
							$seq_ = reverse $seq_original;
							$seq_original=$seq_;
							$seq_original=~tr/ATCG/TAGC/;
							foreach my $key (keys %variants) 
							{ 
								my $key_=length($seq)+$segment_start_first-$key-1; 
								$variants_{$key_}=$variants{$key};
								$variants2_{$key_}=$variants2{$key};
								#print qq!$key->$key_\n!;
							}
						}
						else
						{
							foreach my $key (keys %variants)
							{
								$variants_{$key}=$variants{$key};
								$variants2_{$key}=$variants2{$key};
							}
						}
						for(my $n=0;$n<length($seq_original);$n=$n+3)
						{
							my $triplet = substr($seq_original, $n, 3);
							if (length($triplet)==3)
							{
								if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
								$protein_original.=$mapping{$triplet}; 
							}
						}
						$modified=0;
						my $modified_somatic=0;
						my $stop_found=0;
						my $triplet_count=0;
						my $frame_shift=0;
						#if ($description_=~/\w/) { $description_.="; "; }
						for(my $n=0;$n<length($seq) and $stop_found==0;$n=$n+3)
						{
							my $n_=$n+2;
							my $triplet_count_=$triplet_count+1;
							my $triplet = substr($seq, $n, 3);
							if (length($triplet)==3)
							{
								if ($mapping{$triplet}!~/[\w\*]/) { $mapping{$triplet}="X"; }
								my $triplet_old=$triplet;
								my $triplet_type="";
								my $frame_shift_this=0;
								for (my $i=$n, my $j=0; $i<=$n_; $i++,$j++) 
								{
									#print LOG_MOD qq!$i. $variants_{$i}\n!;
									if ($variants_{$i}=~/^([^#]+)#([^#]+)#([^#]+)$/)
									{
										my $old=$1;
										my $new=$2;
										$frame_shift += abs(length($old)-length($new)) % 3;
										$frame_shift_this = abs(length($old)-length($new)) % 3;
										$triplet_type=$3;
										if ($strand=~/\-/) { $old=~tr/ATCG/TAGC/; } 
										substr $triplet_old,$j,1,$old;
										$variant_loc=$variants2_{$i};
										#print LOG_MOD qq!$i. $triplet $triplet_old $old $variant_loc $variants2_{$i}\n!;
									}
								}

								if ($mapping{$triplet_old}!~/[\w\*]/) { $mapping{$triplet_old}="X"; }
								#print LOG_MOD qq!$name $triplet_type $n-$n_: $triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet}\n!;

								if ($mapping{$triplet}=~/\*/)
								{
									$stop_found=1;
								}

								if ($mapping{$triplet}=~/\*/ and $mapping{$triplet_old}!~/\*/ and $frame_shift==0)
								{
									$protein.=$mapping{$triplet_old}; 
		#							print LOG_MOD      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet}\n!;
									#$description_.=qq!$variant_loc-$triplet_type:$mapping{$triplet_old}$triplet_count_$mapping{$triplet},!;
									$modified++;
									$count_variant_in_exon_nonsyn++;
									$count_stop_introduced++;
									if ($triplet_type=~/S/) 
									{
										$modified_somatic++;
										$count_variant_in_exon_nonsyn_somatic++;
										$count_stop_introduced_somatic++;
									}
								}
								else
								{
									if ($mapping{$triplet}!~/\*/ and $mapping{$triplet_old}=~/\*/ and $frame_shift==0)
									{
									
										$protein.=$mapping{$triplet};	
											 
	#									print LOG_MOD      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet}\n!;
										#$description_.=qq!$variant_loc-$triplet_type:$mapping{$triplet_old}$triplet_count_$mapping{$triplet},!;
										$modified++;
										$count_variant_in_exon_nonsyn++;
										$count_stop_removed++;
										if ($triplet_type=~/S/) 
										{
											$modified_somatic++;
											$count_variant_in_exon_nonsyn_somatic++;
											$count_stop_removed_somatic++;
										}
									}
									else
									{
										if ($mapping{$triplet}!~/^$mapping{$triplet_old}$/ and $frame_shift==0)
										{
											$protein.=$mapping{$triplet}; 
#											print LOG_MOD      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet}\n!;
											#$description_.=qq!$variant_loc-$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet},!;
											$description_.=qq!$variant_loc-$triplet_type:$mapping{$triplet_old}$triplet_count_$mapping{$triplet},!;
											$modified++;
											$count_variant_in_exon_nonsyn++;
											if ($triplet_type=~/S/) 
											{
												$modified_somatic++;
												$count_variant_in_exon_nonsyn_somatic++;
											}
										}
										else
										{
											if ($frame_shift_this!=0)
											{
												$protein.=$mapping{$triplet_old};
												$frame_shift_this=0;
#												print LOG_MOD      qq!$triplet_type $n-$n_:$triplet_old->$triplet - $triplet_count_:$mapping{$triplet_old}->$mapping{$triplet}\n!;
												#$description_.=qq!$variant_loc-$triplet_type:$triplet_count_-frame-shift,!;
												$modified++;
												$count_variant_in_exon_nonsyn++;
												if ($triplet_type=~/S/) 
												{
													$modified_somatic++;
													$count_variant_in_exon_nonsyn_somatic++;
												}
											}
											else
											{
												$protein.=$mapping{$triplet_old};
											}
										}
									}
								}
							}
							$triplet_count++;
						}
						$proteins_count++;
						$protein_modifications_count+=$modified;
						$protein_modifications_distr[$modified]++;
						$protein_modifications_count_somatic+=$modified_somatic;
						$protein_modifications_distr_somatic[$modified_somatic]++;
						$protein_original=~s/\*$//;
						if ($protein_original=~/^([^\*]+)\*.*$/)
						{
							print LOG qq!Error: Stop codon found in middle of sequence:$name \n$protein_original\n!;
							$protein_original="";
						}
					#	if (length($protein_original)>6 and $protein_original!~/^\*/)
					#	{
					#		print OUT qq!>$name $descriptions{$name} (MAP:$chr:$start$strand $segment_lengths $segment_starts)\n$protein_original\n!;
					#	}
						if ($modified!=0)
						{
							$protein=~s/\*$//;
							$protein=~s/^([^\*]+)\*.*$//;
							$proteins_modified++;
							if (length($protein)>6 and $protein!~/^\*/)
							{
								$description_.=")";
#								print LOG_MOD qq!>$name-variant $descriptions{$name} (MAP:$chr:$start$strand $segment_lengths $segment_starts) $description_\n$protein\n!;
								print OUT_MOD qq!>$name-variant $descriptions{$name} (MAP:$chr:$start$strand $segment_lengths $segment_starts) $description_\n$protein\n!;
							}
						if (length($protein_original)>6 and $protein_original!~/^\*/)
                        {
                            print OUT qq!>$name $descriptions{$name} (MAP:$chr:$start$strand $segment_lengths $segment_starts)\n$protein_original\n!;
                        }
						}
					} else { print LOG qq!Error parsing $bed{$name}\n!; }
				}
			} else { print LOG qq!Error in name $chr: $line\n!; }
			
			close(IN);
		}
	}
	print STAT qq!
		$proteins_count proteins
		$proteins_modified proteins modified
		$count_stop_removed stop codons removed
		$count_stop_introduced stop codons introduced
		Somatic: $count_stop_removed_somatic stop codons removed
		Somatic: $count_stop_introduced_somatic stop codons introduced
		$protein_modifications_count total modifications
		$protein_modifications_count_somatic somatic modifications
		Somatic: $somatic_only
		Germline: $germline_only
		In both VCF: $both_vcf (differ: $both_vcf_differ)
		$count_variant_in_exon variants in exons
		$count_variant_in_exon_old_error variants in exons where old does not match genome
		$count_variant_in_exon_nonsyn non-synonymous variants
		Somatic: $count_variant_in_exon_somatic variants in exons
		Somatic: $count_variant_in_exon_old_error_somatic variants in exons where old does not match genome
		Somatic: $count_variant_in_exon_nonsyn_somatic non-synonymous variants
	!;
	print STAT qq!\nnumber of modifications\tnumber of proteins\tnumber of proteins (somatic variants)\n!;
	#for (my $k=0;$k<=100;$k++) { print STAT qq!$k\t$protein_modifications_distr[$k]\t$protein_modifications_distr_somatic[$k]\n!; }
	close(OUT);
	close(OUT_MOD);
	close(OUT_FS);
	close(LOG);
	close(LOG_MOD);
	close(STAT);
}
