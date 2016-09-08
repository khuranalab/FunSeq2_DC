package Funseq_SNV;

use 5.000;
use strict;
our $VERSION = '0.1';
use List::Util qw[min max];

=head1 NAME

=head1 SYNOPSIS

   use Funseq_SNV;
   Funseq_SNV::function

=head1 METHOD

	format_check()
		Input file format checking - BED or VCF format.
		
	new()
		Object constructor.
		
	snv_filter()
		Filter variants based on minor allele frequency threshold. This function also filter out structural variants.
		
	gerp_score()
		Obtain Gerp score for each variant. 
		
	conserved()
		Locate in conserved regions or not.
		
	hot_region()
		Locate in HOT regions (Kevin Yip, Genome Biology, 2012). 
		
	motif_break()
		Check whether variants break bound motifs. 
		
	motif_gain()
		Check whether variants create novel motifs.
		
	sensitive()
		Locate in sensitive/ultra-sensitive regions. 
		
	gene_link()
		Link non-coding mutations with genes (promoters & regulatory elements).
		
	coding()
		Coding pipe - with VAT (variant annotation tool) pipe
		
	intergrate()
		Interage coding & noncoding outputs
	
	recur()
		Recurrent analysis (recurrent in the same position, same gene or same regulatory elements)
	
	fit()
		transformation of continous score

	weight()
		Assign weight for each feature	
=head1 AUTHOR
	Yao Fu; Gerstein Lab; Yale University

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2013 by  Gerstein Lab

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself, either Perl version 5.10.1 or,
at your option, any later version of Perl 5 you may have available.

=cut

## Input format checking 
sub format_check{

	my ($input,$format,$err) = @_;
	my $err_bed = <<EOF;
Error: incorrect bed format; please prepare your file as : chromosome	start	end	reference_allele	alternative_allele
* Note: the 6th column is saved for sample names if having multiple samples.

e.g.,
chr1  213941196  213941197	G	T	sample1	  ...
chr1  213942363  213942364	A	C	sample1	  ...
chr1  213943530  213943531	T	A	sample1	  ...
EOF

	my $err_vcf = <<EOF;
Error: incorrect vcf format; please prepare your file as :

e.g.,
#CHROM POS    ID        REF  ALT     QUAL FILTER INFO ...
2      4370   rs6057    G    A       29   .      ...   ...
EOF

	my $status = 0;
	
	if ($format =~ /bed/i){
		open(IN, $input);
		while(<IN>){
			if (!/\S+\s+\d+\s+\d+\s+\S+\s+\S+/){
				print $_;
				$status = 1;
			}
		}
		if ($status == 1){
			print $err_bed;
			print $err $err_bed;
		}
		return $status;
	}
	
	if ($format =~ /vcf/i){
		open(IN, $input);
		while(<IN>){
			if (!/\S+\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+/ && !/^#/){
				print $_;
				$status=1;
			}
		}
		if ($status == 1){
			print $err_vcf;
			print $err $err_vcf;
		}
		return $status;
	}
}

## Constructor 
sub new {
    my $class=shift;
    my $self ={};
    bless($self,$class);
    return $self;
}
    
## Filter variants against polymorphisms data with minor allele frequency threshold (Note: structural variants are filtered out too)...
sub snv_filter{
# Arguments taken : input file ,input format ,1000 genome file, minor allele frequency input ,output
	my $self = shift;
	my @ARGV = @_;
	my $input = $ARGV[0];   	   #Input File
	my $input_format = $ARGV[1];   #Input File Format
	my $polymorphism = $ARGV[2];   #Polymorphisms data
	my $maf_cut_off  = $ARGV[3];   #Minor allele frequency cut-off
	my $output = $ARGV[4];         #output
	my $out_indel = $ARGV[5];	   #indel output
	my $sv_length_cut = $ARGV[6]; 	   #indel length cut-offs	
	

	my $input_line; 
	my @snp_tabix;
	my $snp;
	my $id;
	my %saw;
	
	open(OUT,">$output")||die;
	
	my $tmp_file = $output.".tmp";
	`sed 's/^chr//gI' $input | sed -e '/^#/! s/^/chr/' > $tmp_file`;
	
	if ($maf_cut_off != 1){

		##
		# use tabix to get subset of variants from input that overlap with 1KG.genome file
		# minor allel frequency is between [1,0], thus, mat_cut_off=0 is removing all 1KG variants, mat_cut_off=1 is keeping all 1KG variants     
		##
		@snp_tabix = split /\n/, `tabix $polymorphism -B $tmp_file| awk '{FS="\t";OFS="\t"} \$4 >= $maf_cut_off'`;
		foreach $snp (@snp_tabix){
			$id = join("\t",(split /\t/,$snp)[0,1]);
			$saw{$id} =1;
		}
	}
	
	if ($input_format =~ /vcf/i){
		open(IN,$tmp_file);
		while(<IN>){
		   	if (!/^#/){
		   		chomp $_;
		   		my @tmp = split /\t/,$_;
		   		# only SNVs
		   		if (length($tmp[3])==1 && length($tmp[4])==1 && $tmp[3] =~ /[ATCG]/i && $tmp[4] =~ /[ATCG]/i){
		   			my $id = join("\t",$tmp[0],$tmp[1]-1);
		   			if(defined $saw{$id}){
		   			}else{
		   				$self->{VCF} ->{$id} = join("\t",@tmp[0..6]);
		   				$self->{DES} ->{$id} = join("\t",$tmp[0],$tmp[1]-1,@tmp[1,3,4]);
		   				print OUT $self -> {DES} -> {$id},"\n";
		   			}
		   		}
		   	}
		}
		close IN;
	}elsif($input_format =~ /bed/i){
		open(INDEL,">$out_indel")||die;
		open(IN,$tmp_file);
		while(<IN>){
			chomp $_;
		   	my @tmp = split /\t/,$_;
		   	# only SNVs
		   	if (length($tmp[3])==1 && length($tmp[4])==1 && $tmp[3] =~ /[ATCG]/i && $tmp[4] =~ /[ATCG]/i){
		   		my $id = join("\t",@tmp[0,1]);
		   		if(defined $saw{$id}){
		   		}else{
		   			$self->{DES} -> {$id} = join("\t",@tmp[0..4]);
		   			print OUT $self -> {DES} -> {$id},"\n";
		   		}
		   	}else{
				if ($sv_length_cut eq 'inf'){
					 print INDEL join("\t",@tmp[0..4]),"\n";
				}elsif(max(length($tmp[3]),length($tmp[4])) <= $sv_length_cut){
		   			print INDEL join("\t",@tmp[0..4]),"\n";
		   		}
			}
		}
		close IN;
		close INDEL;
	}
	close OUT;
	unlink($tmp_file);
}

## Gerp score of mutations (from UCSC)
sub gerp_score{
	my $self = shift;
	my ($gerp_file,$infile) = @_;
	my $id;
	
	foreach $id (sort keys %{$self->{DES}}){
		$self->{GERP}->{$id} = ".";
	}
	if (-s $gerp_file && -f $gerp_file){		
		my $out = `awk '{FS="\t";OFS="\t"}{print \$1,\$2,\$2+1,\$1":"\$2}' $infile |sort -k 1,1 -k 2,2n |uniq| bigWigAverageOverBed $gerp_file stdin stdout`;
		foreach my $line (split /\n/,$out){
			my @tmp = split /:|\s+/,$line;
			$id = join("\t",@tmp[0,1]);
			if ($tmp[4] != 0){
				$self -> {GERP} -> {$id} = $tmp[4];
			}
		}	
	}
}

## conserved regions
sub conserved {
	my $self = shift;
	my ($input,$file) = @_;
	open(IN, "intersectBed -u -a $input -b $file |");
	while(<IN>){
		my $id =join("\t",(split /\s+/,$_)[0,1]);
		$self -> {CONS} ->{$id} =1;
	}
	close IN;
}

## Hot regions (Kevin Yip et al., Genome Biology, 2012)
sub hot_region{
	my $self = shift;
	my ($input,$file) = @_;
	open(IN, "intersectBed -a $input -b $file -wo | ");
	while(<IN>){
		my @tmp = split /\s+/,$_;
		my $id = join("\t",@tmp[0,1]);
		$self -> {HOT} -> {$id} -> {$tmp[8]} =1;
	}
	close IN;
} 

## Motif Analysis 
# 1. Motif Breaking (bound motifs under peak regions)
sub motif_break{
# Arguments taken : input_file, ancestral_file, encode_motif_data, somatic/germline_mode, motif_pfm
	my $self = shift;
	local @ARGV=@_;
	unless (scalar @ARGV==5){
	print <<EOF;
		motif_break: Missing argument...
		out_annotated_nc, informat, ancestral_path, bound_motif, genome_mode, motif_pfm, out_motif_break
EOF
	exit;
	}

	$| = 1;
	
	my $input_file = $ARGV[0];
	my $ancestral_file = $ARGV[1];
	my $bound_motif = $ARGV[2];
	my $mode = $ARGV[3];
	my $pfm = $ARGV[4];

	my $input_file_new;
	my $id;
	my %ancestral;
	my $chr;
	my $pos;
	my $line;
	my @info;
	my $temp;
	my $prev_name;
	my %motif;
	my $A = 1;
	my $C = 2;
	my $G = 3;
	my $T = 4;
	my $ref;
	my $alt;
	my $AA;
	my $der_al;
	my $motif_len;
	my $factor;
	my $pos_in_motif;
	my $motif_name;
	my $derived_allele_freq;
	my $AA_freq;
	my $interval;
	my $diff_freq;
	
# Read in PFM file 
	open MOTIF, $pfm or die;
	while(<MOTIF>){
		chomp $_;
		if(/^>/){
			$prev_name = (split/>|\s+/,$_)[1];
		}else{
			@info = split/\s+/,$_;
			if(not exists $motif{$prev_name}){
				$motif{$prev_name}->[0] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};
			}else{
				$temp = $motif{$prev_name};
				$motif{$prev_name}->[scalar(@$temp)] = {(A=>$info[$A], T=>$info[$T], C=>$info[$C], G=>$info[$G])};			
			}
		}
	}
	close MOTIF;


	$input_file_new = `intersectBed -a $bound_motif -b $input_file -wo | sort -k 1,1 -k 2,2n |uniq | awk '{OFS="\t"}{print \$8,\$9,\$10,\$11,\$12,\$1,\$2,\$3,\$4,\$5,\$6,\$7}'`;
	
	if ($mode == 2){
		my @line = split /\n+/, `printf "$input_file_new" | awk '{FS="\t";OFS="\t"}{gsub("chr","",\$1); print \$1,\$2,\$3}' | fastaFromBed -fi $ancestral_file -bed stdin -fo stdout`;
		foreach my $line (@line){
			if ($line =~ />/){
				($chr, $pos) = (split />|:|-/,$line)[1..2];
				$id = join("","chr",$chr,"\t",$pos);
			}else{
				chomp $line;
				$ancestral{$id} =$line;
			}
		}
	}
	
	foreach my $input_line (split /\n+/,$input_file_new){
		chomp $input_line;
		@info = split /\t+/,$input_line;
		$ref = uc($info[3]); $alt = uc($info[4]);
		$id = join("\t",$info[0],$info[1]);
			
		if($mode == 1){
			$AA = $ref;
			$der_al = $alt;
		}else{
			$AA = uc($ancestral{$id});			
			if($AA eq $ref){
				$der_al = $alt;
			}elsif($AA eq $alt){
				$der_al = $ref;
			}else{
				next;
			}
		}
		
		$motif_len = $info[7]-$info[6];
		$factor = $info[11];
		if($info[10] eq "+"){
			$pos_in_motif = ($info[1]-$info[6]+1);
		}elsif($info[10] eq "-"){
			$pos_in_motif = $info[7]-$info[1];
			$AA =~ tr/ACGTacgt/TGCAtgca/;
			$der_al =~ tr/ACGTacgt/TGCAtgca/;
		}
		
		$motif_name = (split/_\d+mer/,$info[8])[0];
		if(not exists $motif{$motif_name}){
			#print "Motif_Name_Not_Found\t$line\n";
			next;
		}
		
		my $ref = $motif{$motif_name};
		if($motif_len != scalar(@$ref)){
			#print "Motif_Len_Not_Matched\t$line\n";
			next;
		}
		
		$derived_allele_freq = $motif{$motif_name}->[$pos_in_motif-1]->{$der_al};
		$AA_freq = $motif{$motif_name}->[$pos_in_motif-1]->{$AA};
		$interval = join("",$info[5],":",$info[6],"-",$info[7]);
		$self -> {ANNO} ->{$id} -> {"TFM($factor|$info[8]|$interval)"}=1;
		if ($derived_allele_freq < $AA_freq){
			$diff_freq = $AA_freq-$derived_allele_freq;
			if (defined $self -> {MOTIFBR}->{$id}){
				$self -> {MOTIFBR} -> {$id} = join(",",$self->{MOTIFBR}->{$id}, "$factor#$info[8]#$info[6]#$info[7]#$info[10]#$pos_in_motif#$derived_allele_freq#$AA_freq");
				my $max_v = max($diff_freq,$self->{BR_PROB} -> {$id});
				$self->{BR_PROB} -> {$id} = $max_v;
			}else{
				$self->{MOTIFBR}->{$id} = "$factor#$info[8]#$info[6]#$info[7]#$info[10]#$pos_in_motif#$derived_allele_freq#$AA_freq"; 
				$self->{BR_PROB} -> {$id} = $diff_freq;
			}
		}		
	}
}

# 2. Motif gain 
sub motif_gain{
	my $self = shift;
	my ($pfm_file,$reference_file,$score_file,$p_cut,$out_tmp) = @_; 
	
	my $id;
	my $line;
	my $snp_file = "";
	my @des;
	my @alt;
	my @ref;
	my @id;
	my @start;
	
	my @info;
	my %ref = ();
	my %motif = ();
	my $prev_name;
	my $A = 1;
	my $C = 2;
	my $G = 3;
	my $T = 4;
	my $temp;
	my $extract_seq;
	my @extract_seq;
	my $chr;
	my $start; 
	my $end;
	my $alt;
	my $ref;
	my $alt_pwm_score;
	my $ref_pwm_score;
	my $motif_length;
	my $pwm;
	my $i;
	my $j;
	my $alt_p;
	my $ref_p;
	my $strand;
	
	$|++;
	
		
	# Read in PFM and transform it to PWM using 20 sequences. 

	open MOTIF, "$pfm_file" or die;
	while($line = <MOTIF>){
		chomp($line);
		if($line =~ /^>/){
			$prev_name = (split/>|\s+/,$line)[1];
		}else{
			@info = split/\s+/,$line;
			if(not exists $motif{$prev_name}){
				$motif{$prev_name}->[0] = {(A=>log(($info[$A]*20+0.25)/21)-log(0.25), T=>log(($info[$T]*20+0.25)/21)-log(0.25), C=>log(($info[$C]*20+0.25)/21)-log(0.25), G=>log(($info[$G]*20+0.25)/21)-log(0.25))};
			}else{
				$temp = $motif{$prev_name};
				$motif{$prev_name}->[scalar(@$temp)] = {(A=>log(($info[$A]*20+0.25)/21)-log(0.25), T=>log(($info[$T]*20+0.25)/21)-log(0.25), C=>log(($info[$C]*20+0.25)/21)-log(0.25), G=>log(($info[$G]*20+0.25)/21)-log(0.25))};			
			}
		}
	}
	close MOTIF;

	# Read in Score file corresponding.
	
	my %score_lower;
	my %score;
	my %score_upper;
	
	open SCORE,"$score_file";
	while(<SCORE>){
		my ($prev_name,$cut_off,$p) = (split /\s+/,$_)[0..2];
		if ($p < $p_cut){
			$score_upper{$prev_name} = $cut_off;	
		}elsif ($p == $p_cut){
			$score{$prev_name} = $cut_off;
		}elsif($p > $p_cut){
			$score_lower{$prev_name} = $cut_off;
		}
	}
	close SCORE;
	

	
	# retrieve + & - 29bp around the SNP & gain of motif calculation; 
	
	foreach $id(sort keys %{$self->{NGENE}}){
		my $switch = 0;
		foreach my $gene(keys %{$self -> {NGENE}->{$id}}){
			if (defined $self ->{NGENE} ->{$id} ->{$gene} -> {"Promoter"} || defined $self ->{NGENE}->{$id}->{$gene}->{"Distal"}){
				$switch = 1;
			}
		}
		if ($switch == 1){
			@des = split /\s+/,$self->{DES}->{$id};
			push @id,$id;
			$chr = $des[0];
			$chr =~ s/chr//g;
			$start = $des[1];
			$snp_file .= join("",$chr,"\t",$start-29,"\t",$start+30,"\n");
			push @ref,uc($des[3]);
			push @alt,uc($des[4]);
			push @start,$start;
		}
	}

	open(O,">$out_tmp")||die;
	print O $snp_file;
	close O;
	
	if (-s $out_tmp){
	
		@des = split /\n/, `fastaFromBed -fi $reference_file -bed $out_tmp -fo stdout`;
	
		if (scalar @des >0){
			for $i (0 .. (scalar @des/2)-1){
				$ref = $ref[$i]; $alt = $alt[$i];
				$id = $id[$i]; $start = $start[$i];
				$extract_seq=uc($des[$i*2+1]);
				@extract_seq = split //,$extract_seq;
				$extract_seq[29] = $alt;	
			
			#positive strand
				&seq_scan(\@extract_seq,"+");
		
			#negative strand
				$extract_seq =~ tr/ATCGatcg/TAGCTAGC/;
				@extract_seq = reverse(split //,$extract_seq);
				$ref =~ tr/ATCGatcg/TAGCTAGC/; 
				$alt =~ tr/ATCGatcg/TAGCTAGC/;
				$extract_seq[29] = $alt;
				&seq_scan(\@extract_seq,"-");
			}
		}
	}
	
	
# sub-routine for the sequence scanning ... 
	
	sub seq_scan{
		my $seq;
		($seq, $strand) = @_;
		my @seq;

		
		@seq = @{$seq};
		
		foreach $prev_name(sort keys %motif){
			#print $prev_name,"\n";

			$motif_length = scalar @{$motif{$prev_name}};
	
			unless (-d "pwm"){
				mkdir "pwm";
			}
			
			if (-e "pwm/$prev_name"){
			}else{
				open(TMP,">pwm/$prev_name")||die;
				$pwm = "";
				for $i(0 .. $motif_length-2){
					$pwm .= join("",$motif{$prev_name}->[$i]->{"A"},"\t");
				}
				$pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"A"},"\n");
				for $i(0 .. $motif_length-2){
					$pwm .= join("",$motif{$prev_name}->[$i]->{"C"},"\t");
				}
				$pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"C"},"\n");
				for $i(0 .. $motif_length-2){
					$pwm .= join("",$motif{$prev_name}->[$i]->{"G"},"\t");
				}
				$pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"G"},"\n");
				for $i(0 .. $motif_length-2){
					$pwm .= join("",$motif{$prev_name}->[$i]->{"T"},"\t");
				}
				$pwm .= join("",$motif{$prev_name}->[$motif_length-1]->{"T"},"\n");	
				print TMP $pwm;
				close TMP;
			}			
					
			# sequence scanning ....
			
			for ($i=30 - $motif_length; $i < 30; $i ++){
				$alt_pwm_score = 0;
				$ref_pwm_score = 0;
				
				# calculate reference & alternative PWM score;
				for ($j = 0; $j < $motif_length; $j++){
					$alt_pwm_score += $motif{$prev_name}->[$j]->{$seq[$i+$j]};
				}
				$ref_pwm_score = $alt_pwm_score - $motif{$prev_name}->[29-$i]->{$alt} + $motif{$prev_name}->[29-$i]->{$ref};
				
				if ($alt_pwm_score > $ref_pwm_score){
					if (defined $score{$prev_name}){
						if ($alt_pwm_score > $score{$prev_name} && $ref_pwm_score < $score{$prev_name}){
							&out_print($prev_name,$i,$motif_length,$alt_pwm_score,$ref_pwm_score,$strand,$id,$start);
						}else{
							next;
						}
					}elsif ((defined $score_upper{$prev_name} && $ref_pwm_score >= $score_upper{$prev_name}) || (defined $score_lower{$prev_name} && $alt_pwm_score <= $score_lower{$prev_name})){
						next;
					}elsif(defined $score_upper{$prev_name} && defined $score_lower{$prev_name} && $alt_pwm_score >= $score_upper{$prev_name} && $ref_pwm_score <= $score_lower{$prev_name}){
						&out_print($prev_name,$i,$motif_length,$alt_pwm_score,$ref_pwm_score,$strand,$id,$start);	
					}else{
						if (defined $score_upper{$prev_name} && $alt_pwm_score > $score_upper{$prev_name}){
								$ref_p = (split /\s+/,`TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $ref_pwm_score -m pwm/$prev_name -w`)[2];
								if ($ref_p > $p_cut){
									&out_print($prev_name,$i,$motif_length,$alt_pwm_score,$ref_pwm_score,$strand,$id,$start);	
								}
						}else{	
								$alt_p = (split /\s+/,`TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $alt_pwm_score -m pwm/$prev_name -w`)[2];
								$ref_p = (split /\s+/,`TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $ref_pwm_score -m pwm/$prev_name -w`)[2];
								if ($alt_p < $p_cut && $ref_p > $p_cut){
									&out_print($prev_name,$i,$motif_length,$alt_pwm_score,$ref_pwm_score,$strand,$id,$start);	
								}
						}
					}	
				}
			}					
		}		
	}

#end of routine
	
	sub out_print{
		my ($name,$i,$motif_length,$alt_pwm_score,$ref_pwm_score,$strand,$id,$start) = @_;
		my $diff_score = $alt_pwm_score - $ref_pwm_score;
		if ($strand eq "+"){
			my $tmp = join("","$name#",$start-29+$i,"#",$start+$i+$motif_length-29, "#$strand#",30-$i,"#",sprintf('%.3f', $alt_pwm_score),"#",sprintf('%.3f', $ref_pwm_score));
			if (defined $self->{MOTIFG}->{$id}){
				$self->{MOTIFG} ->{$id} = join(",",$self->{MOTIFG}->{$id},$tmp);
				my $max_v = max($diff_score,$self -> {G_PROB} -> {$id});
				$self -> {G_PROB} -> {$id} = $max_v;
			}else{
				$self -> {MOTIFG}->{$id} = $tmp;
				$self -> {G_PROB} -> {$id} = $diff_score;
			}
		}else{
			my $tmp = join("", "$name#",$start-$motif_length-$i+30,"#",$start-$i+30, "#$strand#",30-$i,"#",sprintf("%.3f", $alt_pwm_score),"#",sprintf("%.3f", $ref_pwm_score));
			if (defined $self ->{MOTIFG}->{$id}){
				$self->{MOTIFG} ->{$id} = join(",",$self->{MOTIFG}->{$id},$tmp);
				my $max_v = max($diff_score,$self -> {G_PROB} -> {$id});
				$self -> {G_PROB} -> {$id} = $max_v;
			}else{
				$self -> {MOTIFG}->{$id} = $tmp;
				$self -> {G_PROB} -> {$id} = $diff_score;
			}
		}
	}
	
}

## Sensitive
sub sensitive{
	my $self =shift;
	my ($input,$sensitive) = @_;
	my $id;
	
	open(IN, "intersectBed -u -a  $input -b $sensitive |");
	while(<IN>){
		$id = join("\t",(split /\t+/,$_)[0,1]);
		$self->{SEN}->{$id} =1;			
	}
	close IN;
	
	open(IN, "grep 'Ultra' $sensitive| intersectBed -u -a  $input -b stdin |");
	while(<IN>){
		$id = join("\t",(split /\t+/,$_)[0,1]);
		$self -> {USEN}->{$id}=1;
	}
	close IN;
}

## ENCODE Annotation
sub read_encode{
	my $self = shift;
	my ($out_nc,$encode_annotation)= @_;
	open(IN,"intersectBed -a $out_nc -b $encode_annotation -wo | ");
	while(<IN>){
		my @tmp = split /\t+/,$_;
		my $id = join("\t",@tmp[0,1]);
		my $interval = join("",$tmp[5],":",$tmp[6],"-",$tmp[7]);			
		my $tag = "";
		my @anno = split /\./,$tmp[8];
		if ($anno[-1] eq 'u'){
			$tag = $anno[0].'('.join('.',@anno[1..$#anno-1]).')';
		}else{
			$tag = $anno[0].'('.join('.',@anno[1..$#anno])."|$interval)";
		}
		$self->{ANNO}->{$id} -> {$tag}=1;
	}
	close IN;
}

## User Annotation
sub user_annotation{
	my $self = shift;
	my ($snv_input,$anno_dir)= @_;
	
	opendir (DIR, $anno_dir);
	if( grep ! /^\.\.?$/, readdir DIR){
		open(A,"ls $anno_dir/* |");
		while(<A>){
			chomp $_;
			my $file = $_;
			my $cate = uc((split /\./, ((split /\//,$file)[-1]))[0]);
			open(IN,"intersectBed -a $snv_input -b $file -wo | ");
			while(<IN>){
				my @tmp = split /\t+/,$_;
				my $id = join("\t",@tmp[0,1]);
				my $interval = join("",$tmp[5],":",$tmp[6],"-",$tmp[7]);			
				my $tag;
				if(scalar @tmp > 9){			
					$tag = $cate.'('.$tmp[8]."|$interval)";
				}else{
					$tag = $cate."($interval)";
				}
				$self->{USER} -> {$id} -> {$tag}=1;
			}
			close IN;
		}
		close A;
	}
	closedir(DIR);
}

## Promoter,Intron,UTR,Enhancer  (network centrality ... )
sub gene_link{
	my $self = shift;
	my ($input,$promoter,$distal,$intron,$utr,$network) = @_;
	my %network;
	my %net_degree = ();
	
	my @pro =split /\n/, `intersectBed -a $input -b $promoter -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @dis = split /\n/, `intersectBed -a $input -b $distal -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @intron = split /\n/, `intersectBed -a $input -b $intron -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @utr = split /\n/, `intersectBed -a $input -b $utr -wo | sort -k 1,1 -k 2,2n |uniq`;
	
	open(A,"ls $network/* |");
	while(<A>){
		chomp $_;
		my $file = $_;
		my $net = (split /\./, ((split /\//,$file)[-1]))[0];
		open(NET,$file)||die;
		while(<NET>){
			if (!/GENE_NAME/){
				my ($gene,$degree) = (split /\s+/,$_)[0,1];
				$network{$gene}{$net} =$degree ;
				if(not exists $net_degree{$net}){
					$net_degree{$net}[0] = $degree;
				}else{
					push @{$net_degree{$net}},$degree; 
				}
			}
		}
		close NET;
	}
	close A;
	
	&func(\@intron,"Intron");
	&func(\@utr,"UTR");
	&func(\@pro,"Promoter");
	&func(\@dis,"Distal");

	
	sub func{
		my ($cate,$tag) = @_;
		foreach my $line(@$cate){
			my @tmp = split /\s+/,$line;
			my $id = join("\t",@tmp[0,1]);
			my $gene = $tmp[8];
			
			$self->{NGENE}->{$id}->{$gene}->{$tag}=1;
			
			if ($tag eq "Promoter" && scalar @tmp > 11){
				$self -> {LINK} -> {$id} -> {$gene} -> {$tag} -> {$tmp[10]} =1;
			}
			if ($tag eq "Distal" && scalar @tmp > 10){
				$self -> {LINK} -> {$id} -> {$gene} -> {$tag} -> {$tmp[9]} =1;
			}
		 
			if (defined $network{$gene}){
				my @prob_gene = ();
				foreach my $net (sort keys %{$network{$gene}}){
					my $degree = $network{$gene}{$net};
					my $greater = scalar grep { $_ > $degree } @{$net_degree{$net}};					
					my $prob = sprintf("%.3f", 1- $greater/scalar @{$net_degree{$net}});	
					push @prob_gene, "$net($prob)";
					if (defined $self -> {NET_PROB} -> {$id}){
						$self -> {NET_PROB} -> {$id} = max($prob,$self -> {NET_PROB}->{$id}); 
					}else{
						$self -> {NET_PROB} -> {$id} = $prob;
					}
				}
				my $tmp_hub = $gene.":".join('',@prob_gene);
				$self->{NHUB} -> {$id} -> {$tmp_hub}=1;
			}
		}
	}
}

## Coding : VAT (variant annotation tool) analysis
sub coding{
	my ($snp_input,$file_interval,$file_fasta,$network,$file_selection,$out,$out_path,$nc_mode,$cds)=@_;
	my @tmp;
	my %hub;
	my %selection;
	my $gene;
	my @vat_out;
	my $vat_out;
	my @nvat_out;
	my $nvat_out;
	
# run VAT
	
	if($nc_mode==0){
		
		@vat_out = split /\n/, `awk '{FS="\t";OFS="\t"}BEGIN{print "##fileformat=VCFv4.0\\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}{print \$1,\$3,".",\$4,\$5,".","PASS","."}' $snp_input | snpMapper $file_interval $file_fasta`;
		
		# Read in networks
		my %network;
		my %net_degree = ();
		open(A,"ls $network/* |");
		while(<A>){
			chomp $_;
			my $file = $_;
			my $net = (split /\./, ((split /\//,$file)[-1]))[0];
			open(NET,$file)||die;
			while(<NET>){
				if (!/GENE_NAME/){
					my ($gene,$degree) = (split /\s+/,$_)[0,1];
					$network{$net}{$gene} =$degree ;
					if(not exists $net_degree{$net}){
						$net_degree{$net}[0] = $degree;
					}else{
						push @{$net_degree{$net}},$degree; 
					}
				}
			}
			close NET;
		}
		close A;
	
		foreach my $net (keys %net_degree){
			my @tmp_degree = sort {$b<=>$a} @{$net_degree{$net}};
			my $cut = $tmp_degree[int(scalar @tmp_degree*0.25)-1];
			foreach my $gene(keys %{$network{$net}}){
				if ($network{$net}{$gene} >=$cut){
					$hub{$gene}{$net} = 1;
				}
			}
		}



		open(IN,$file_selection)||die;
		while(<IN>){
			if (!/GENE_NAME/){
				$gene = (split /\s+/,$_)[0];
				$selection{$gene}=1;
			}
		}
		close IN;


		if(scalar @vat_out ==0){
			@nvat_out = split /\n/, `intersectBed -a $snp_input -b $cds -wo | cut -f 1,3,9 | sort |uniq`;
			open(O,">$out")||die;
			foreach $nvat_out (@nvat_out){
				@tmp = split /\s+/,$nvat_out;
				print O $tmp[0],"\t",$tmp[1],"\t",".\t",$tmp[2],"\t";
				$gene = $tmp[2];
				if (defined $hub{$gene}){
					print O join('&',sort keys %{$hub{$gene}}),"\t";
				}else{
					print O ".\t";
				}
				if (defined $selection{$gene}){
					print O "NegativeSelection\n";
				}else{
					print O "\n";
				}
			}
			close O;
			
			return "1";
		}else{
			open(O,">$out")||die;
			foreach $vat_out(@vat_out){
				if ($vat_out !~ /^#/){
					@tmp =split /\s+/,$vat_out;
					print O $tmp[0],"\t",$tmp[1],"\t",$tmp[7],"\t";
					my $info = (split /,/,$tmp[7])[0];
					if ($info =~ /VA=\d+:(.*?):/){
						$gene =$1;
						print O $gene,"\t";
						if (defined $hub{$gene}){
							print O join('&',sort keys %{$hub{$gene}}),"\t";
						}else{
							print O ".\t";
						}
						if (defined $selection{$gene}){
							print O "NegativeSelection\n";
						}else{
							print O "\n";
						}
					}
				}
			}
			close O;
			
			if ((split /\s+/, `awk 'END{print FNR}' $snp_input`)[0] > scalar @vat_out){
				@nvat_out = split /\n/, `awk '{OFS="\t"}{print \$1,\$2-1,\$2}' $out | intersectBed -a $snp_input -b stdin -v | intersectBed -a stdin -b $cds -wo | cut -f 1,3,9 | sort | uniq`;
				open(O,">>$out")||die;
				foreach $nvat_out (@nvat_out){
					@tmp = split /\s+/,$nvat_out;
					print O $tmp[0],"\t",$tmp[1],"\t",".\t",$tmp[2],"\t";
					$gene = $tmp[2];
					if (defined $hub{$gene}){
						print O join('&',sort keys %{$hub{$gene}}),"\t";
					}else{
						print O ".\t";
					}
					if (defined $selection{$gene}){
						print O "NegativeSelection\n";
					}else{
						print O "\n";
					}
				}
			}
			close O;
			return "0";

		}
	}else{
		return "0";
	}
}

## Integration of all outputs
sub intergrate{
	my $self=shift;
	local @ARGV= @_;
	unless (scalar @ARGV == 10){
print <<EOF;
	intergrate.pl : Missing arguments ...
EOF
exit;		
}

	$| = 1; #Auto Flush

	my $output_format = $ARGV[0]; # Output Format
	my $sample = $ARGV[1]; 		  # sample name 
	my $nc_snp = $ARGV[2];        # noncoding_snp 
	my $coding_info = $ARGV[3];   # coding output
	my $gene_dir = $ARGV[4];      # Prior knowledge of genes
	my $reg_net = $ARGV[5];       # regulatory network
	my $out = $ARGV[6];           # output
	my $de_data = $ARGV[7];       # deferentially expression output 
	my $weight_mode = $ARGV[8]; 
	my $weight_file = $ARGV[9];
	  
	open(OUT,">$out")||die;
	my %nc_score;
	my %cds_score;
	my %gene_info;
	
	
	## Read in gene info 
	open(A,"ls $gene_dir/* |");
	while(<A>){
		chomp $_;
		my $file = $_;
		my $gene_cate = join("",'[',(split /\./, ((split /\//,$file)[-1]))[0],']');
		open(GENE,$file)||die;
		while(<GENE>){
			my $gene = (split /\s+/,$_)[0];
			$gene_info{$gene}{$gene_cate} = 1;
		}
		close GENE;
	}
	close A;
 	
 	if (-e $de_data && -s $de_data ){
 		open(DE,$de_data);
 		while(<DE>){
 			chomp $_;
 			my ($gene_li,$cls) = (split /\t+/,$_)[0,1];
 			foreach my $gene(split /\|/,$gene_li){
 				if ($cls =~ /benign/i){
 					$gene_info{$gene}{"[down_regulated]"}=1;
 				}else{
 					$gene_info{$gene}{"[up_regulated]"}=1;
 				}
 			}
 		}
 	}
 	close DE;
 	
 	
 	my %cancer_reg;
 	open(REG,$reg_net);
 	while(<REG>){
 		my ($tf,$gene) = (split /\s+/,$_)[0,1];
 		if (defined $gene_info{$gene} && defined $gene_info{$gene}{"[cancer]"}){
 			$cancer_reg{$tf}{$gene} =1;
 		}
 	}
 	close REG;
 	
 	foreach my $tf(keys %cancer_reg){
 		my $info = join("",'[TF_regulating_known_cancer_gene:',join(",",sort keys %{$cancer_reg{$tf}}),']');
 		$gene_info{$tf}{$info} = 1;
 	}
 	
 	undef %cancer_reg;
 	# done 
 	
 	
 	&read_nc($nc_snp);
	&read_cds($coding_info);

	if ($weight_mode ==1){
		&nc_score_scheme();
	}else{
		&nc_score_scheme_uw();
	}

	&rank_output($output_format);
 
 
 	sub read_nc{
		my ($file) = @_;
		open(IN,$file);
		while(<IN>){	
			my $id = join("\t",(split /\s+/,$_)[0,1]);
			$self ->{NC} ->{$id} =1;
		}
		close IN;
	}
   

	sub read_cds{
		my ($file) = @_;
		open(IN,$file);
		while(<IN>){
			chomp $_;
			my @tmp = split /\t+/,$_;
			my $id = join("\t",$tmp[0],$tmp[1]-1);
			my $vat = (split /;/,$tmp[2])[-1];
			$vat =~ s/nonsynonymous/missense_variant/g;
			$vat =~ s/synonymous/synonymous_variant/g;
			$vat =~ s/prematureStop/stop_gained/g;
			$vat =~ s/removedStop/stop_lost/g;
			$vat =~ s/spliceOverlap/splice_variant/g;
			
			my $score = 0;
			$self -> {VAT} ->{$id}= $vat;
			$self ->{GENE} -> {$id} = $tmp[3];

			# add spliceOverlap into CDS scoring scheme
			if (/nonsynonymous/ || /prematureStop/ || /spliceOverlap/){
			#if (/nonsynonymous/ || /prematureStop/){
				$score ++;
				
				if (/prematureStop/ || /spliceOverlap/){
				#if (/prematureStop/){
					$score ++;
				}

				if ($tmp[4] ne "."){
					$self -> {HUB} ->{$id} = $tmp[4];
					$score ++;
				}

				if (/NegativeSelection/){
					$self -> {SELECTION} ->{$id} = 1;
					$score ++;
				}

				if (defined $self ->{CONS} ->{$id}){
    				$score++;
    			}

    			if ($self ->{GERP} -> {$id} ne "." && $self ->{GERP} ->{$id} > 2){
    				$score++;    # gerp greater than 2
    			}
			}
			$cds_score{$score}{$id} = 1;
		}
		close IN;
	}


	
	
	sub nc_score_scheme{
		my $score;
		my $id;
		my $return_v = &weight($weight_file);
		my %weight = %{$return_v};
		
		foreach $id (keys %{$self -> {NC}}){
			$score = 0;
			if (defined $self->{ANNO}->{$id} && (defined $self -> {SEN}->{$id})!=1 && (defined $self -> {MOTIFBR}->{$id})!=1 && (defined $self -> {HOT}->{$id})!=1){
				$score = $weight{ANNO};
			}

			if (defined $self -> {USEN} ->{$id}){
				$score = $weight{USEN};
			}elsif (defined $self -> {SEN} -> {$id}){
				$score = $weight{SEN};
			}

			if (defined $self -> {MOTIFBR} -> {$id}){
				my $tmp_score = $weight{MOTIFBR};
				$tmp_score =~ s/value/$self->{BR_PROB}->{$id}/;	
				$score +=eval($tmp_score);	
			}

			if (defined $self -> {HOT} -> {$id}){
    			$score += $weight{HOT};
			}
			
			if (defined $self ->{CONS} ->{$id}){
    			$score += $weight{CONS};
    		}elsif ($self ->{GERP} -> {$id} ne "."){
				my $tmp_score = $weight{GERP};
				$tmp_score =~ s/value/$self->{GERP}->{$id}/;	
				$score +=eval($tmp_score);
    		}

    		if (defined $self ->{NGENE} -> {$id} && (defined $self ->{NHUB} -> {$id}) != 1 && (defined $self ->{MOTIFG}->{$id}) !=1){
    			$score += $weight{NGENE};
    		}

			if (defined $self ->{NHUB} ->{$id}){
				my $tmp_score = $weight{NHUB};
				$tmp_score =~ s/value/$self->{NET_PROB}->{$id}/;	
				$score +=eval($tmp_score);
			}

			if (defined $self -> {MOTIFG}->{$id}){
				if ($self->{G_PROB}->{$id} < (log(20.25)-log(0.25))){
					my $tmp_score = $weight{MOTIFG};
					$tmp_score =~ s/value/$self->{G_PROB}->{$id}/;
					$score +=eval($tmp_score);
				}else{
					$score += 0.99629183; 
				}
			}
			
			$nc_score{$score}{$id} = 1;
		}
	}
	sub nc_score_scheme_uw{
		my $score;
		my $id;
		
		foreach $id (keys %{$self -> {NC}}){
			$score = 0;
			if (defined $self->{ANNO}->{$id}){
				$score ++ ;
			}
			if (defined $self -> {SEN} -> {$id}){
				$score ++;
			}
			if (defined $self -> {USEN} ->{$id}){
				$score ++;
			}
			if (defined $self ->{CONS} ->{$id}){
    			$score ++;
    		}
    		if ($self ->{GERP} -> {$id} ne "." && $self -> {GERP} -> {$id} > 2){
    			$score ++;
    		}
    		if (defined $self -> {HOT} -> {$id}){
    			$score ++;
			}
			if (defined $self ->{NHUB} ->{$id}){
				if ($self ->{NET_PROB} -> {$id} > 0.75){
					$score ++;
				}
			}
			if (defined $self -> {MOTIFBR} -> {$id}){
					$score ++;
			}
			if (defined $self -> {MOTIFG}->{$id}){
				$score ++;
			}
			if (defined $self -> {NGENE} -> {$id}){
				$score ++;
			}
			$nc_score{$score}{$id} = 1;
		}
	}
	
	sub rank_output{
		# rank coding infront of noncoding 
		my ($format) = @_;
		if ($format =~ /bed/i){
    		&print_bed_format();
    	}elsif($format =~ /vcf/i){
			&print_vcf_format();
   		}else{
    		return "Please specify proper output format";	
    	}
	}
	
	
	sub print_bed_format{
		# * Coding bed format ***
		
		my $score;
		my $id;
    	foreach $score (sort {$b<=>$a} keys %cds_score){
    		foreach $id (sort keys %{$cds_score{$score}}){
    			print OUT $self -> {DES} -> {$id},"\t",$sample,"\t";
    			print OUT $self -> {GERP} -> {$id},";";
    			print OUT "Yes;";
    			print OUT $self -> {VAT} ->{$id},";";
    			
    			if (defined $self ->{HUB} ->{$id}){
    				#print OUT "HUB=",join(",", sort keys %{$self -> {HUB} -> {$id}}),";";
                                print OUT $self ->{HUB} ->{$id},";";
    			}else{
    				print OUT ".;";
    			}

    			if (defined $self ->{SELECTION} ->{$id}){
    				print OUT "Yes;";
    			}else{
    				print OUT ".;";
    			}
    			
    			
    			if (defined $self -> {ANNO} -> {$id}){
    				print OUT join(",",sort keys %{$self -> {ANNO} -> {$id}}),";";
    			}else{
    				print OUT ".;";
    			}

    			if (defined $self -> {HOT} ->{$id}){
    				print OUT join(",",sort keys %{$self ->{HOT} -> {$id}}),";";
    			}else{
    				print OUT ".;";
    			}
    			
    			if (defined $self -> {MOTIFBR} -> {$id} || defined $self -> {MOTIFG} -> {$id}){
    				if (defined $self -> {MOTIFBR} -> {$id}){
    					print OUT "MOTIFBR=",$self -> {MOTIFBR} -> {$id};
    					if (defined $self -> {MOTIFG}-> {$id}){
    						print OUT ",MOTIFG=",$self -> {MOTIFG} -> {$id}, ";";
    					}else{
    						print OUT ";";
    					}
    				}else{
    					print OUT "MOTIFG=",$self -> {MOTIFG}->{$id},";";
    				}
    			}else{
    				print OUT ".;";
    			}

    			print OUT ".;"x2;
    			
    			if (defined $self ->{CONS} ->{$id}){
    				print OUT "Yes;";
    			}else{
    				print OUT ".;";
    			}

    			my $gene = $self -> {GENE} -> {$id};
    			
    			if (defined $gene_info{$gene}){
    				print OUT $gene,join("",sort keys %{$gene_info{$gene}}),";";
    			}else{
    				print OUT $gene,";";
    			}

    			if (defined $self -> {USER}){
    				if (defined $self ->{USER} ->{$id}){
    					print OUT join(",",sort keys %{$self ->{USER}->{$id}}),";";
    				}else{
    					print OUT ".;";
    				}
    			}

    			print OUT $score,";";
    			print OUT ".\n";
    		}
    	}
   
    	# * Non-coding bed format **** 
    
    	foreach $score (sort {$b<=>$a} keys %nc_score){
    		foreach $id (sort keys %{$nc_score{$score}}){
    			print OUT $self -> {DES} ->{$id},"\t",$sample,"\t";
    			print OUT $self -> {GERP} -> {$id},";";
    			print OUT "No;.;";

    			if (defined $self -> {VAT} -> {$id}){
    				print OUT $self -> {VAT} -> {$id},";";
				}

    			if (defined $self -> {NHUB} -> {$id}){
    				print OUT join(",",sort keys %{$self -> {NHUB} -> {$id}}),";";
    			}else{
    				print OUT ".;";
    			}

    			print OUT ".;";   # negative selection
    		
    			if (defined $self -> {ANNO} -> {$id}){
    				print OUT join(",",sort keys %{$self -> {ANNO} -> {$id}}),";";
    			}else{
    				print OUT ".;";
    			}

    			if (defined $self -> {HOT} ->{$id}){
    				print OUT join(",",sort keys %{$self ->{HOT} -> {$id}}),";";
    			}else{
    				print OUT ".;";
    			}
    			
    			if (defined $self -> {MOTIFBR} -> {$id} || defined $self -> {MOTIFG} -> {$id}){
    				if (defined $self -> {MOTIFBR} -> {$id}){
    					print OUT "MOTIFBR=",$self -> {MOTIFBR} -> {$id};
    					if (defined $self -> {MOTIFG}-> {$id}){
    						print OUT ",MOTIFG=",$self -> {MOTIFG} -> {$id}, ";";
    					}else{
    						print OUT ";";
    					}
    				}else{
    					print OUT "MOTIFG=",$self -> {MOTIFG}->{$id},";";
    				}
    			}else{
    				print OUT ".;";
    			}
    		
    			if (defined $self -> {SEN} -> {$id}){
    				print OUT "Yes;";
    			}else{
    				print OUT ".;";
    			}

    			if (defined $self -> {USEN} -> {$id}){
    				print OUT "Yes;";
    			}else{
    				print OUT ".;";
    			}

    			if (defined $self ->{CONS} ->{$id}){
    				print OUT "Yes;";
    			}else{
    				print OUT ".;";
    			}
    			
    			if (defined $self -> {NGENE} -> {$id}){
    				my $gene_info = "";
    				foreach my $gene (sort keys %{$self -> {NGENE} -> {$id}}){
    					
    					my $info = "$gene(";
						foreach my $tag (sort keys %{$self -> {NGENE} -> {$id}->{$gene}}){
							if (defined $self -> {LINK} ->{$id} && defined $self -> {LINK} ->{$id} -> {$gene} && defined $self -> {LINK} ->{$id} ->{$gene} -> {$tag}){
								$info = $info.$tag.'['.join(',',sort keys %{$self -> {LINK} ->{$id} ->{$gene} -> {$tag}}).']&';
							}else{
								$info = $info.$tag.'&';
							}
						}
    					$info =~ s/&$/)/;
    					
    					if (defined $gene_info{$gene}){
    						$info = $info.join("",sort keys %{$gene_info{$gene}});
    					}
    					$gene_info = join("",$gene_info,$info,",");
    				}
    				
    				$gene_info =~ s/,+$//;
    				print OUT $gene_info,";";
    			}else{
    				print OUT ".;";
    			}

    			if (defined $self -> {USER}){
    				if (defined $self ->{USER} ->{$id}){
    					print OUT join(",",sort keys %{$self ->{USER}->{$id}}),";";
    				}else{
    					print OUT ".;";
    				}
    			}

    			print OUT ".;";
    			print OUT $score,"\n";
    		}
    	}
	}

	sub print_vcf_format{
	# * Coding vcf format **** 
		my $score;
		my $id;

    	foreach $score (sort {$b<=>$a} keys %cds_score){
    		foreach $id (sort keys %{$cds_score{$score}}){
    			if (defined $self -> {VCF} ->{$id}){
    				print OUT $self -> {VCF} ->{$id},"\t";
    			}else{
    				my @tmp = split /\t+/,$self->{DES}->{$id};
    				print OUT $tmp[0],"\t",$tmp[2],"\t.\t",$tmp[3],"\t",$tmp[4],"\t",".\t.\t";
    			}
    			print OUT "SAMPLE=$sample;";
    			print OUT "GERP=",$self -> {GERP}->{$id},";";
    			
    			print OUT "CDS=Yes",";";
    			print OUT $self -> {VAT} -> {$id},";";
    			if (defined $self -> {HUB}->{$id}){
                                #print OUT "HUB=",join(",", sort keys %{$self -> {HUB} -> {$id}}),";";
    				print OUT "HUB=",$self -> {HUB}->{$id},";";
    			}
    			if (defined $self -> {SELECTION} -> {$id}){
    				print OUT "GNEG=Yes;";
    			}
    			
    			if (defined $self -> {ANNO} -> {$id}){
    				print OUT "NCENC=",join(",", sort keys %{$self -> {ANNO}->{$id}}),";";
    			}
    			if (defined $self -> {HOT} ->{$id}){
    				print OUT "HOT=",join(",",sort keys %{$self -> {HOT} ->{$id}}),";";
    			}
    			if (defined $self ->{MOTIFBR} ->{$id}){
    				print OUT "MOTIFBR=",$self ->{MOTIFBR}->{$id},";";
    			}
    			if (defined $self ->{MOTIFG} ->{$id}){
    				print OUT "MOTIFG=",$self ->{MOTIFG}->{$id},";";
    			}
    			
    			if (defined $self ->{CONS} ->{$id}){
    				print OUT "UCONS=Yes;";
    			}
    			print OUT "GENE=",$self -> {GENE} -> {$id},";";
    			
    			if (defined $gene_info{$self->{GENE}->{$id}}){
    				print OUT "CANG=",join("",sort keys %{$gene_info{$self->{GENE}->{$id}}}),";";
    			}
    			if (defined $self -> {USER} && defined $self ->{USER} ->{$id}){
    					print OUT "USER_ANNO=",join(",",sort keys %{$self ->{USER}->{$id}}),";";
    			}
    			print OUT "CDSS=",$score,"\n";
    		}
    	}
    
   
    	# *  Non-coding vcf format **** 
    
    	foreach $score (sort {$b<=>$a} keys %nc_score){
    		foreach $id (sort keys %{$nc_score{$score}}){
    			if (defined $self -> {VCF} ->{$id}){
    				print OUT $self -> {VCF} ->{$id},"\t";
    			}else{
    				my @tmp = split /\t+/,$self -> {DES} ->{$id};
    				print OUT $tmp[0],"\t",$tmp[2],"\t.\t",$tmp[3],"\t",$tmp[4],"\t",".\t.\t";
    			}

    			print OUT "SAMPLE=$sample;";
    			print OUT "GERP=",$self -> {GERP}->{$id},";";
    			
    			print OUT "CDS=No;";
    			
    			if (defined $self -> {VAT} -> {$id}){
    				print OUT $self -> {VAT} -> {$id},";";
				}

    			if (defined $self -> {NHUB} -> {$id}){
    				print OUT "HUB=",join(",", sort keys %{$self -> {NHUB} -> {$id}}),";";
    			}

    			if (defined $self -> {ANNO} -> {$id}){
    				print OUT "NCENC=",join(",", sort keys %{$self -> {ANNO}->{$id}}),";";
    			}

    			if (defined $self -> {HOT} ->{$id}){
    				print OUT "HOT=",join(",",sort keys %{$self -> {HOT} ->{$id}}),";";
    			}

    			if (defined $self ->{MOTIFBR} ->{$id}){
    				print OUT "MOTIFBR=",$self ->{MOTIFBR}->{$id},";";
    			}

    			if (defined $self ->{MOTIFG} ->{$id}){
    				print OUT "MOTIFG=",$self ->{MOTIFG}->{$id},";";
    			}

    			if (defined $self ->{SEN}->{$id}){
    				print OUT "SEN=Yes;";
    			}

    			if (defined $self ->{USEN}->{$id}){
    				print OUT "USEN=Yes;";
    			}

    			if (defined $self ->{CONS} ->{$id}){
    				print OUT "UCONS=Yes;";
    			}

    			if (defined $self ->{NGENE}->{$id}){    				
    				my $gene_info = "";
    				my $info = "";
    				foreach my $gene (sort keys %{$self -> {NGENE} -> {$id}}){
    				
    				    $info = $info."$gene(";
						foreach my $tag (sort keys %{$self -> {NGENE} -> {$id}->{$gene}}){
							if (defined $self -> {LINK} ->{$id} && defined $self -> {LINK} ->{$id} -> {$gene} && defined $self -> {LINK} ->{$id} ->{$gene} -> {$tag}){
								$info = $info.$tag.'['.join(',',sort keys %{$self -> {LINK} ->{$id} ->{$gene} -> {$tag}}).']&';
							}else{
								$info = $info.$tag.'&';
							}
						}
    					$info =~ s/&$/),/;
    					    				     
    				    if (defined $gene_info{$gene}){
    						$gene_info = join("",$gene_info,$gene,join("",sort keys %{$gene_info{$gene}}),",");
    					}
    				}
    				
    				$info =~ s/,+$//;
    				print OUT "GENE=$info;";
    				$gene_info =~ s/,+$//;

    				if ($gene_info ne ""){
    					print OUT "CANG=$gene_info;";
    				}
    			}

    			if (defined $self -> {USER} && defined $self ->{USER} ->{$id}){
    					print OUT "USER_ANNO=",join(",",sort keys %{$self ->{USER}->{$id}}),";";
    			}

    			print OUT "NCDS=",$score,"\n";
    		}
    	}
	}

}

## Recurrence Analysis
sub recur{
	my ($out_detail,$out_recur,$out_driver,$outformat,$cancer_dir,$cancer_type,$score_cut,$anno_dir,$weight_file,$recur_db_score) = @_;
	my $cancer_recur;
	my $return_v = &weight($weight_file);
	my %weight = %{$return_v};


	# read in recurrence data 
	my $type_check = 0;
	sub read_cancer{
		open(A,"ls $cancer_dir/* |");
		while(<A>){
			chomp $_;
			my $file = $_;
			my $type = (split /\./, ((split /\//,$file)[-1]))[0];
			if ($type eq $cancer_type || $cancer_type eq "all"){
				$type_check =1;
			}
			open(CAN,"$file")||die;
			while(<CAN>){
				chomp $_;
				my($elm,$info) = (split /\t/,$_)[0,1];
				$cancer_recur ->{$elm} -> {$type} = $info;
			}
			close CAN;
		}
		close A;
	}
	
	if (-d $cancer_dir){
		&read_cancer();
	}else{
		print "WARNING: File directory $cancer_dir not exists ...\n";
	}
	
	if ($type_check == 0){
		print "WARNING: $cancer_type not found ...\n";
	}
	# done 
	
	my $vcf_header = <<EOF;
##fileformat=VCFv4.0
##INFO=<ID=OTHER,Number=.,Type=String, Description = "Other Information From Original File">
##INFO=<ID=SAMPLE,Number=.,Type=String,Description="Sample id">
##INFO=<ID=CDS,Number=.,Type=String,Description="Coding Variants or not">
##INFO=<ID=VA,Number=.,Type=String,Description="Coding Variant Annotation">
##INFO=<ID=HUB,Number=.,Type=String,Description="Network Hubs, PPI (protein protein interaction network), REG (regulatory network), PHOS (phosphorylation network)...">
##INFO=<ID=GNEG,Number=.,Type=String,Description="Gene Under Negative Selection">
##INFO=<ID=GERP,Number=.,Type=String,Description="Gerp Score">
##INFO=<ID=NCENC,Number=.,Type=String,Description="NonCoding ENCODE Annotation">
##INFO=<ID=HOT,Number=.,Type=String,Description="Highly Occupied Target Region">
##INFO=<ID=MOTIFBR,Number=.,Type=String,Description="Motif Breaking">
##INFO=<ID=MOTIFG,Number=.,Type=String,Description="Motif Gain">
##INFO=<ID=SEN,Number=.,Type=String,Description="In Sensitive Region">
##INFO=<ID=USEN,Number=.,Type=String,Description="In Ultra-Sensitive Region">
##INFO=<ID=UCONS,Number=.,Type=String,Description="In Ultra-Conserved Region">
##INFO=<ID=GENE,Number=.,Type=String,Description="Target Gene (For coding - directly affected genes ; For non-coding - promoter, enhancer regulatory module)">
##INFO=<ID=CANG,Number=.,Type=String,Description="Prior Gene Information, e.g.[cancer][TF_regulating_known_cancer_gene][up_regulated][actionable]...";
##INFO=<ID=CDSS,Number=.,Type=String,Description="Coding Score">
##INFO=<ID=NCDS,Number=.,Type=String,Description="NonCoding Score">
##INFO=<ID=RECUR,Number=.,Type=String,Description="Recurrent elements / variants">
##INFO=<ID=DBRECUR,Number=.,Type=String,Description="Recurrence database">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
EOF
	
	my $bed_header;
	opendir (DIR, $anno_dir);
	if( grep ! /^\.\.?$/, readdir DIR){
		$bed_header=join("\t","#chr","start","end","ref","alt","sample","gerp;cds;variant.annotation.cds;network.hub;gene.under.negative.selection;ENCODE.annotated;hot.region;motif.analysis;sensitive;ultra.sensitive;ultra.conserved;target.gene[known_cancer_gene/TF_regulating_known_cancer_gene,differential_expressed_in_cancer,actionable_gene];user.annotations;coding.score;noncoding.score;recurrence;recurrence.cancer.dataset");
	}else{
		$bed_header=join("\t","#chr","start","end","ref","alt","sample","gerp;cds;variant.annotation.cds;network.hub;gene.under.negative.selection;ENCODE.annotated;hot.region;motif.analysis;sensitive;ultra.sensitive;ultra.conserved;target.gene[known_cancer_gene/TF_regulating_known_cancer_gene,differential_expressed_in_cancer,actionable_gene];coding.score;noncoding.score;recurrence;recurrence.cancer.dataset");
	}
	closedir(DIR);
	
	
	my $id;
	my @temp;
	my $tag;
	my $id_info=();
	my $anno;
	my $elm_info;
	my %saw_tag;
	
	open(IN,$out_detail);
	if ($outformat =~ /bed/i){
		while(<IN>){
			@temp = split /\s+/,$_;
			$id = join(":",@temp[0,2]);
			$tag = $temp[5];
			$saw_tag{$tag}=1;
			$id_info->{$id}->{$tag} =1;
			my @feature = split /;/,$temp[6];
			if ($feature[1] eq "Yes" && (/missense_variant/ || /stop_gained/)){
				$anno = $feature[11];
				$elm_info -> {$anno}->{$tag}->{$id}=1;
			}elsif($feature[1] eq "No"){
				foreach $anno((split /,/,$feature[5])){
					if ($anno ne "."){
						$elm_info -> {$anno} -> {$tag} -> {$id}=1;
					}
				}
			}
		}
	}else{
		while(<IN>){
			@temp = split /\s+/,$_;
			$id = join(":",@temp[0,1]);
			
			if (/SAMPLE=(.*?);/){
				$tag = $1;
			}
			$saw_tag{$tag}=1;
			$id_info->{$id}->{$tag} =1;
			if (/CDS=Yes/ && (/missense_variant/ || /stop_gained/)){
				if (/GENE=(.*?);/){
					$anno = $1;
					$elm_info -> {$anno} -> {$tag} -> {$id}= 1;
				}
			}elsif(/CDS=No/ && /NCENC=/){
				if (/NCENC=(.*?);/){
					foreach $anno((split /,/,$1)){
						if ($anno ne "."){
							$elm_info -> {$anno} -> {$tag} -> {$id}=1;
						}
					}
				}
			}
		}
	}

	
	open(JOIN,">$out_recur")||die;
	my $recur_score;
	my %recur_elm;
	my %recur_id;
	my $score;
	
	foreach $anno (sort keys %{$elm_info}){
		$score = scalar keys %{$elm_info -> {$anno}};
		if ($score > 1){
			$recur_score -> {$score} -> {$anno} =1;
		}
	}
	foreach $id (keys %{$id_info}){
		$score = scalar keys %{$id_info -> {$id}};
		if ($score >1){
			$recur_score -> {$score} -> {$id} =1;
		}
	}
	
	foreach $score (sort {$b <=> $a} keys %{$recur_score}){
		foreach $anno (sort keys %{$recur_score -> {$score}}){

			if(defined $elm_info -> {$anno}){
				print JOIN $anno,"\t";
				print JOIN "Altered in ",scalar keys %{$elm_info -> {$anno}},'/',scalar keys %saw_tag, '(', sprintf("%.2f", ((scalar keys %{$elm_info -> {$anno}})/scalar keys %saw_tag)*100),"%) samples.\t"; 
				my $line = "";
				foreach $tag(sort keys %{$elm_info -> {$anno}}){
					$line = $line."$tag(";
					foreach $id (sort keys %{$elm_info->{$anno}->{$tag}}){
						if (scalar keys %{$id_info->{$id}} > 1){
							$line = $line."$id*,";
						}else{
							$line = $line."$id,";
						}
					}
					$line =~ s/,$/\),/;
				}
				$line =~ s/,$//;
				print JOIN $line,"\n";
				$recur_elm{$anno} = join("",$anno,":",join('&',sort keys %{$elm_info->{$anno}}));
			}else{
				$id = $anno;
				print JOIN "$id\t";
				print JOIN "Altered in ",scalar keys %{$id_info -> {$id}},'/',scalar keys %saw_tag, '(',sprintf("%.2f", ((scalar keys %{$id_info->{$id}})/scalar keys %saw_tag)*100),"%) samples.\t"; 
				print JOIN "$id*(",join(",",keys %{$id_info->{$id}}),")\n";
				$recur_id{$id} = join ("", "$id:",join('&',sort keys %{$id_info->{$id}}));
			}		
		}
	}
	close JOIN;
	
	#Done 

	my @db_recur;
	my %recur_cancer;
	sub db_recur{
		my ($x) = @_;
		if (defined $cancer_recur -> {$x}){
			if ($cancer_type eq "all"){
				my $tmp_line = "$x:";
				foreach my $type(sort keys %{$cancer_recur -> {$x}}){
					$tmp_line = join("",$tmp_line,$type,"(",$cancer_recur->{$x}->{$type},")|");
					$recur_cancer{$type} = 1;
				}
				$tmp_line=~ s/\|$//;
				push @db_recur,$tmp_line;
			}else{
				if (defined $cancer_recur->{$x}->{$cancer_type}){
					push @db_recur, join("","$x:",$cancer_type,"(",$cancer_recur->{$x}->{$cancer_type},")");
					$recur_cancer{$cancer_type} = 1;
				}
			}
		}
	}
	
	
	open(DRIVER, ">$out_driver");	
	open(OUT,">$out_detail.process");
	
	open(IN,$out_detail);
	if ($outformat =~ /bed/i){
		print OUT $bed_header,"\n";
		while(<IN>){
			chomp $_;
			@temp = split /\s+/,$_;
			$id = join(":",@temp[0,2]);
			$tag = $temp[5];
			
			###
			@db_recur = ();
			%recur_cancer = ();
			my $db_recur = ".";
			my $recur_cancer = ".";
			###
			
			my @feature = split (/;/,$temp[6]);
			if ($feature[1] eq "Yes"){
				$anno = $feature[11];
				
				###
				&db_recur($anno);
				&db_recur($id);
				if(scalar @db_recur >0){
					$db_recur=join(",",@db_recur);
				}
				if(scalar keys %recur_cancer >0){
					$recur_cancer=join(",",sort keys %recur_cancer);
				}
				###
				
				if (defined $recur_elm{$anno}){
					print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-2]),";",$feature[$#feature-1]+1,";",$feature[$#feature],";",$recur_elm{$anno},";",$db_recur,"\n";
					if (/missense_variant/ || /stop_gained/){
						print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-2]),";",$feature[$#feature-1]+1,";",$feature[$#feature],";Yes;$recur_cancer\n";
					}
				}elsif ($db_recur ne "." && $recur_db_score ==1){
					print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-2]),";",$feature[$#feature-1]+1,";",$feature[$#feature],";.;",$db_recur,"\n";
					if (/missense_variant/ || /stop_gained/){
						print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-2]),";",$feature[$#feature-1]+1,";",$feature[$#feature],";.;$recur_cancer\n";
					}
				}else{
					print OUT $_,";.;$db_recur\n";
					if (/missense_variant/ || /stop_gained/){
						print DRIVER $_,";.;$recur_cancer\n";
					}
				}
			}else{
				my $hot_simp;
				if ($feature[6] eq "."){
					$hot_simp=".";
				}else{
					$hot_simp="Yes";
				}
				
				my $motif_simp;
				if ($feature[7] eq "."){
					$motif_simp = ".";
				}else{
					if ($feature[7] =~ /MOTIFBR=(.*?),MOTIFG=(.*?)$/){
						my @motif = split /,/,$1;
						my %motif_simp = ();
						foreach my $motif(@motif){
							$motif_simp{(split /#/,$motif)[0]}=1;
						}
						my @motif = split /,/,$2;
						$motif_simp = join("","MOTIFBR:",join(",",sort keys %motif_simp),",MOTIFG:");
						my %motif_simp = ();
						foreach my $motif(@motif){
							$motif_simp{(split /#/,$motif)[0]}=1;
						}
						$motif_simp = $motif_simp.join(",",sort keys %motif_simp);
						
					}elsif ($feature[7]=~ /MOTIFBR=(.*?)$/){
						my @motif = split /,/,$1;
						my %motif_simp = ();
						foreach my $motif(@motif){
							$motif_simp{(split /#/,$motif)[0]}=1;
						}
						$motif_simp = "MOTIFBR:".join(",",sort keys %motif_simp);
						
					}elsif ($feature[7]=~ /MOTIFG=(.*?)$/){
						my @motif = split /,/,$1;
						my %motif_simp = ();
						foreach my $motif(@motif){
							$motif_simp{(split /#/,$motif)[0]}=1;
						}
						$motif_simp = "MOTIFG:".join(",",sort keys %motif_simp);
					}
				}
				
				if ($feature[5] eq "."){
				
				###
					&db_recur($id);
					if(scalar @db_recur >0){
						$db_recur=$db_recur[0];
					}
					if(scalar keys %recur_cancer >0){
						$recur_cancer=join(",",sort keys %recur_cancer);
					}
				###
				
					if(defined $recur_id{$id}){	
						print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-1]),";", $feature[$#feature],':',$feature[$#feature]+$weight{RECUR},";",$recur_id{$id},";$db_recur\n";
						if (($feature[13]+$weight{RECUR})>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..5]),";$hot_simp;$motif_simp;",join(";",@feature[8 ..$#feature-1]),";",$feature[$#feature]+$weight{RECUR},";Yes;$recur_cancer\n";
						}
						
					}elsif($db_recur ne "." && $recur_db_score ==1){
						print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-1]),";", $feature[$#feature],':',$feature[$#feature]+$weight{RECUR},";.;$db_recur\n";
						if (($feature[13]+$weight{RECUR})>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..5]),";$hot_simp;$motif_simp;",join(";",@feature[8 ..$#feature-1]),";",$feature[$#feature]+$weight{RECUR},";.;$recur_cancer\n";
						}
					}else{
						if (($feature[13])>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..5]),";$hot_simp;$motif_simp;",join(";",@feature[8 ..$#feature]),";.;$recur_cancer\n";
						}
						print OUT $_,";.;$db_recur\n";
					}
					
				}else{
					my @encode_simp = ();
					
					if ($feature[5] =~ /DHS\(/){
						push @encode_simp,"DHS";
					}
					if ($feature[5] =~ /TFP\(/){
						push @encode_simp,"TFP";
					}
					if ($feature[5] =~ /TFM\(/){
						push @encode_simp, "TFM";
					}
					if ($feature[5] =~ /Pseudogene\(/){
						push @encode_simp,"Pseudogene";
					}
					if ($feature[5] =~ /Enhancer\(/){
						push @encode_simp,"Enhancer";
					}
					if ($feature[5] =~ /RNA\(/){
						push @encode_simp,"ncRNA";
					}
					
					my $line = "";
					foreach $anno((split /,/,$feature[5])){
						###
						&db_recur($anno);
						###
						if (defined $recur_elm{$anno}){
							$line = join("",$line,$recur_elm{$anno},",");
						}
					}
					if (defined $recur_id{$id}){
                                                $line = join("",$line,$recur_id{$id},",");
                                        }       
	
					###
					&db_recur($id);
					if(scalar @db_recur >0){
						$db_recur=join(",",@db_recur);
					}
					if(scalar keys %recur_cancer >0){
						$recur_cancer=join(",",sort keys %recur_cancer);
					}
					###
					
					if (length($line )>0){
						$line =~ s/,$//;
						print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-1]),";",$feature[$#feature],':',$feature[$#feature]+$weight{RECUR},";",$line,";$db_recur\n";
						if (($feature[13]+$weight{RECUR})>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..4]),";",join(",",@encode_simp),";$hot_simp;$motif_simp;",join(";", @feature[8..$#feature-1]),";",$feature[$#feature]+$weight{RECUR},";Yes;$recur_cancer\n";
						}
					}elsif(defined $recur_id{$id}){	
						if (($feature[13]+$weight{RECUR})>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..4]),";",join(",",@encode_simp),";$hot_simp;$motif_simp;",join(";", @feature[8..$#feature-1]),";",$feature[$#feature]+$weight{RECUR},";Yes;$recur_cancer\n";
						}
						print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-1]),";",$feature[$#feature],':',$feature[$#feature]+$weight{RECUR},";",$recur_id{$id},";$db_recur\n";
					}elsif($db_recur ne "." && $recur_db_score ==1){
						if (($feature[13]+$weight{RECUR})>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..4]),";",join(",",@encode_simp),";$hot_simp;$motif_simp;",join(";", @feature[8..$#feature-1]),";",$feature[$#feature]+$weight{RECUR},";.;$recur_cancer\n";
						}
						print OUT join("\t",@temp[0..5]),"\t",join(";",@feature[0..$#feature-1]),";",$feature[$#feature],':',$feature[$#feature]+$weight{RECUR},";.;$db_recur\n";
					}else{
						if (($feature[13])>=$score_cut || /\[known_cancer_gene\]/){
							print DRIVER join("\t",@temp[0..5]),"\t",join(";",@feature[0..4]),";",join(",",@encode_simp),";$hot_simp;$motif_simp;",join(";", @feature[8..$#feature]),";.;$recur_cancer\n";
						}									
						print OUT $_,";.;$db_recur\n";
					}
				}
			}
		}
	}else{

		## re-processing VCF output with cancer gene and db recurrence annotation

		print OUT $vcf_header;

		while(<IN>){
			@temp = split /\s+/,$_;
			$id = join(":",$temp[0],$temp[1]);
			if (/SAMPLE=(.*?);/){
				$tag = $1;
			}
			chomp $_;
			###
			@db_recur = ();
			my $db_recur = "";
			%recur_cancer = ();
			my $recur_cancer = "";
			###
			
			if (/;CDSS=(\d+)$/){
				$score = $1;
				if (/GENE=(.*?);/){
					$anno = $1;
					
					###
					&db_recur($anno);
					&db_recur($id);
					if(scalar @db_recur >0){
						$db_recur=join(",",@db_recur);
						$db_recur = ";DBRECUR=".$db_recur;
					}
					if(scalar keys %recur_cancer >0){
						$recur_cancer=join(",",sort keys %recur_cancer);
						$recur_cancer = ";DBRECUR=".$recur_cancer;
					}
					###
					
					if (defined $recur_elm{$anno}){							
						s/;CDSS=\d+$//g;
						#print OUT $_,";CDSS=",$score+1,";RECUR=",$recur_elm{$anno},"$db_recur\n";
						print OUT $_,";CDSS=",$score,':',$score+1,";RECUR=",$recur_elm{$anno},"$db_recur\n";
						
						if (/missense_variant/ || /stop_gained/){
							print DRIVER $_,";CDSS=",$score+1,";RECUR=Yes$recur_cancer\n";
						}
						
					}elsif(scalar @db_recur >0 && $recur_db_score ==1){
						s/;CDSS=\d+$//g;
						#print OUT $_,";CDSS=",$score+1,"$db_recur\n";
						print OUT $_,";CDSS=",$score,':',$score+1,"$db_recur\n";
						if (/missense_variant/ || /stop_gained/){
							print DRIVER $_,";CDSS=",$score+1,";$recur_cancer\n";
						}
					}else{
						print OUT $_,"$db_recur\n";
						
						if (/missense_variant/ || /stop_gained/){
							print DRIVER $_,"$recur_cancer\n";
						}
					}
				}
			}elsif(/NCDS=(\S+)$/){
				$score = $1;
				my $motif_simp = "";
				my $motifbr_simp="";
				my $motifg_simp="";
				my $hot_simp="";
			
				if (/HOT=/){
					$hot_simp = "HOT=Yes;";
				}
				
				if(/MOTIFBR=(.*?);/){
					my @motif = split /,/,$1;
					my %motif_simp = ();
					foreach my $motif(@motif){
						 $motif_simp{(split /#/,$motif)[0]}=1;
					}
					$motifbr_simp = join("","MOTIFBR=",join(",",sort keys %motif_simp),";");
				}
				if(/MOTIFG=(.*?);/){
					my @motif = split /,/,$1;
					my %motif_simp = ();
					foreach my $motif(@motif){
						$motif_simp{(split /#/,$motif)[0]}=1;
					}
					$motifg_simp = join("","MOTIFG=",join(",",sort keys %motif_simp),";");
				}
				$motif_simp = join("",$motifbr_simp,$motifg_simp);
				
				if (/NCENC=(.*?);/){
					my $line = ""; 
					my @encode_simp = ();
					my $encode_info = $1;
					if ($encode_info =~ /DHS\(/){
						push @encode_simp,"DHS";
					}
					if ($encode_info =~ /TFP\(/){
						push @encode_simp,"TFP";
					}
					if ($encode_info =~ /TFM\(/){
						push @encode_simp, "TFM";
					}
					if ($encode_info =~ /Pseudogene\(/){
						push @encode_simp,"Pseudogene";
					}
					if ($encode_info =~ /Enhancer\(/){
						push @encode_simp,"Enhancer";
					}
					if ($encode_info =~ /RNA\(/){
						push @encode_simp,"ncRNA";
					}
					
					
					foreach $anno((split /,/,$encode_info)){
					###
					&db_recur($anno);
					###
						if (defined $recur_elm{$anno}){
							$line = join("", $line,$recur_elm{$anno},",");
						}
					}
				   	
					if (defined $recur_id{$id}){
                                                $line = join("",$line,$recur_id{$id},",");
                                        }
	
					###
					&db_recur($id);
					if(scalar @db_recur >0){
						$db_recur=join(",",@db_recur);
						$db_recur = ";DBRECUR=".$db_recur;
					}
					if(scalar keys %recur_cancer >0){
						$recur_cancer=join(",",sort keys %recur_cancer);
						$recur_cancer = ";DBRECUR=".$recur_cancer;
					}
					###
					
					if (length($line) >0 ){
						s/;NCDS=\S+$//g;	
						$line =~ s/,$//;
						print OUT $_,";NCDS=",$score,':',$score+$weight{RECUR},";RECUR=",$line,"$db_recur\n";
						if($score+$weight{RECUR} >=$score_cut || /\[known_cancer_gene\]/){
							s/MOTIFG=.*?;|MOTIFBR=.*?;|NCENC=.*?;|HOT=.*?;|MOTIFG=.*?$|MOTIFBR=.*?$|NCENC=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";NCENC=",join(",",@encode_simp),";",$hot_simp,$motif_simp,"NCDS=",$score+$weight{RECUR},";RECUR=Yes$recur_cancer\n";
						}
					}elsif(defined $recur_id{$id}){
						s/;NCDS=\S+$//g;	
						print OUT $_,";NCDS=",$score,':',$score+$weight{RECUR},";RECUR=",$recur_id{$id},"$db_recur\n";
						if($score+$weight{RECUR} >=$score_cut || /\[known_cancer_gene\]/){
							s/MOTIFG=.*?;|MOTIFBR=.*?;|NCENC=.*?;|HOT=.*?;|MOTIFG=.*?$|MOTIFBR=.*?$|NCENC=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";NCENC=",join(",",@encode_simp),";",$hot_simp,$motif_simp,"NCDS=",$score+$weight{RECUR},";RECUR=Yes$recur_cancer\n";
						}
					}elsif(scalar @db_recur >0 && $recur_db_score ==1){
						s/;NCDS=\S+$//g;	
						print OUT $_,";NCDS=",$score,':',$score+$weight{RECUR},"$db_recur\n";
						if($score+$weight{RECUR} >=$score_cut || /\[known_cancer_gene\]/){
							s/MOTIFG=.*?;|MOTIFBR=.*?;|NCENC=.*?;|HOT=.*?;|MOTIFG=.*?$|MOTIFBR=.*?$|NCENC=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";NCENC=",join(",",@encode_simp),";",$hot_simp,$motif_simp,"NCDS=",$score+$weight{RECUR},";$recur_cancer\n";
						}
					}else{							
						print OUT $_,"$db_recur\n";
						if($score >=$score_cut || /\[known_cancer_gene\]/){
							s/;NCDS=\S+$//g;	
							s/MOTIFG=.*?;|MOTIFBR=.*?;|NCENC=.*?;|HOT=.*?;|MOTIFG=.*?$|MOTIFBR=.*?$|NCENC=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";NCENC=",join(",",@encode_simp),";",$hot_simp,$motif_simp,"NCDS=",$score,"$recur_cancer\n";
						}
					}
				}else{
					###
					&db_recur($id);
					if(scalar @db_recur >0){
						$db_recur=";DBRECUR=".$db_recur[0];
					}
					if(scalar keys %recur_cancer >0){
						$recur_cancer=join(",",sort keys %recur_cancer);
						$recur_cancer = ";DBRECUR=".$recur_cancer;
					}
					###
					
					if(defined $recur_id{$id}){
						s/;NCDS=\S+$//g;	
						print OUT $_,";NCDS=",$score,':',$score+$weight{RECUR},";RECUR=",$recur_id{$id},"$db_recur\n";
						if($score+$weight{RECUR} >=$score_cut || /\[known_cancer_gene\]/){
							s/MOTIFG=.*?;|HOT=.*?;|MOTIFG=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";",$hot_simp,$motif_simp,"NCDS=",$score+$weight{RECUR},";RECUR=Yes$recur_cancer\n";
						}
					}elsif(scalar @db_recur >0 && $recur_db_score ==1){	
						s/;NCDS=\S+$//g;	
						print OUT $_,";NCDS=",$score,':',$score+$weight{RECUR},"$db_recur\n";
						if($score+$weight{RECUR} >=$score_cut || /\[known_cancer_gene\]/){
							s/MOTIFG=.*?;|HOT=.*?;|MOTIFG=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";",$hot_simp,$motif_simp,"NCDS=",$score+$weight{RECUR},";$recur_cancer\n";
						}
					}else{					
						print OUT $_,"$db_recur\n";
						if($score >=$score_cut || /\[known_cancer_gene\]/){
							s/;NCDS=\S+$//g;
							s/MOTIFG=.*?;|HOT=.*?;|MOTIFG=.*?$|HOT=.*?$//g;
							s/;+$//;
							print DRIVER $_,";",$hot_simp,$motif_simp,"NCDS=",$score,"$recur_cancer\n";
						}
					}
				}
			}
		}				
	}
	close DRIVER;
	close OUT;
	close IN;
	`mv $out_detail.process $out_detail`;	
}


## Weight Assign
sub weight{
	my %weight;
	my ($file) = (@_);
	open(A,$file)||die;
	while(<A>){
		chomp $_;
		my ($cate,$func) = (split /\t/,$_)[0,1];
		$weight{$cate} = $func;
	} 
	return \%weight;
}



1;

