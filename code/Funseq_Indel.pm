package Funseq_Indel;

use 5.000;
use strict;
our $VERSION = '0.2';
use List::Util qw[min max];


## Constructor 
sub new {
    my $class=shift;
    my $self ={};
    bless($self,$class);
    return $self;
}

## Gerp score of mutations (from UCSC)
sub gerp_score{
	my $self = shift;
	my ($gerp_file,$infile) = @_;
	my $id;

	open(IN,$infile)||die;
	while(<IN>){
		$id = join("\t",(split /\s+/,$_)[0..4]);
		$self -> {GERP} -> {$id}  = ".";
	}
	close IN;
	
	if (-s $gerp_file && -f $gerp_file){		
		my $out = `awk '{FS="\t";OFS="\t"}{print \$1,\$2,\$3,\$1":"\$2":"\$3":"\$4":"\$5}' $infile |sort -k 1,1 -k 2,2n |uniq| bigWigAverageOverBed $gerp_file stdin stdout`;
		foreach my $line (split /\n/,$out){
			my @tmp = split /:|\s+/,$line;
			$id = join("\t",@tmp[0..4]);
			if ($tmp[-1] != 0){
				$self -> {GERP} -> {$id} = $tmp[-1];
			}
		}	
	}
}

## annotations
sub annotations {
	my $self = shift;
	my ($input,$conserved,$hot,$sensitive,$encode,$anno_dir,$bound_motif) = @_;
	my $id;
	my @tmp;
	
	open(IN, "intersectBed -u -a $input -b $conserved |");
	while(<IN>){
		$id =join("\t",(split /\s+/,$_)[0..4]);
		$self -> {CONS} ->{$id} =1;
	}
	close IN;

	open(IN, "intersectBed -a $input -b $hot -wo | ");
	while(<IN>){
		@tmp = split /\s+/,$_;
		$id = join("\t",@tmp[0..4]);
		$self -> {HOT} -> {$id} -> {$tmp[8]} =1;
	}
	close IN;
	
	open(IN, "intersectBed -u -a  $input -b $sensitive |");
	while(<IN>){
		$id = join("\t",(split /\s+/,$_)[0..4]);
		$self->{SEN}->{$id} =1;			
	}
	close IN;
	
	open(IN, "grep 'Ultra' $sensitive| intersectBed -u -a $input -b stdin |");
	while(<IN>){
		$id = join("\t",(split /\s+/,$_)[0..4]);
		$self -> {USEN}->{$id}=1;
	}
	close IN;

	open(IN,"intersectBed -a $input -b $encode -wo | ");
	while(<IN>){
		my @tmp = split /\s+/,$_;
		my $id = join("\t",@tmp[0..4]);
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

	opendir (DIR, $anno_dir);
	if( grep ! /^\.\.?$/, readdir DIR){
		open(A,"ls $anno_dir/* |");
		while(<A>){
			chomp $_;
			my $file = $_;
			my $cate = uc((split /\./, ((split /\//,$file)[-1]))[0]);
			open(IN,"intersectBed -a $input -b $file -wo | ");
			while(<IN>){
				my @tmp = split /\t+/,$_;
				my $id = join("\t",@tmp[0..4]);
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
	
	open(IN, "intersectBed -a $bound_motif -b $input -wo | sort -k 1,1 -k 2,2n |uniq | awk '{OFS=\"\t\"}{print \$8,\$9,\$10,\$11,\$12,\$1,\$2,\$3,\$4,\$5,\$6,\$7}' | ");
	while(<IN>){
		chomp $_;
		my @tmp = split /\t+/,$_;
		my $id = join("\t",@tmp[0..4]);
		my $interval = join("",$tmp[5],":",$tmp[6],"-",$tmp[7]);
		$self -> {ANNO} ->{$id} -> {"TFM($tmp[11]|$tmp[8]|$interval)"}=1;
		
		if (defined $self -> {MOTIFBR}->{$id}){
			$self -> {MOTIFBR} -> {$id} = join(",",$self->{MOTIFBR}->{$id}, "$tmp[11]#$tmp[8]#$tmp[6]#$tmp[7]#$tmp[10]");
		}else{
			$self->{MOTIFBR}->{$id} = "$tmp[11]#$tmp[8]#$tmp[6]#$tmp[7]#$tmp[10]"; 
		}
	}
}

## Promoter,Intron,UTR,Enhancer  (network centrality ... )
sub gene_link{
	my $self = shift;
	my ($input,$promoter,$distal,$intron,$utr,$network,$cds,$selection) = @_;
	my %network;
	my %net_degree = ();
	my %selection;
	
	my @pro =split /\n/, `intersectBed -a $input -b $promoter -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @dis = split /\n/, `intersectBed -a $input -b $distal -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @intron = split /\n/, `intersectBed -a $input -b $intron -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @utr = split /\n/, `intersectBed -a $input -b $utr -wo | sort -k 1,1 -k 2,2n |uniq`;
	my @cds = split /\n/, `intersectBed -a $input -b $cds -wo | sort -k 1,1 -k 2,2n |uniq`;
	
	
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
	
	open(IN,$selection)||die;
	while(<IN>){
		if (!/GENE_NAME/){
			my $gene = (split /\s+/,$_)[0];
			$selection{$gene}=1;
		}
	}
	close IN;



	&func(\@intron,"Intron");
	&func(\@utr,"UTR");
	&func(\@pro,"Promoter");
	&func(\@dis,"Distal");
	&func(\@cds,"Coding");
	
	sub func{
		my ($cate,$tag) = @_;
		foreach my $line(@$cate){
			my @tmp = split /\s+/,$line;
			my $id = join("\t",@tmp[0..4]);
			my $gene = $tmp[8];
			
			$self->{GENE}->{$id}->{$gene}->{$tag}=1;
			
			if ($tag eq "Promoter" && scalar @tmp > 11){
				$self -> {LINK} -> {$id} -> {$gene} -> {$tag} -> {$tmp[10]} =1;
			}
			if ($tag eq "Distal" && scalar @tmp > 10){
				$self -> {LINK} -> {$id} -> {$gene} -> {$tag} -> {$tmp[9]} =1;
			}
		 
			if (defined $selection{$gene}){
				$self ->{SELECTION} -> {$id} ->{$gene} =1;
			}
			
			if (defined $network{$gene}){
				my @prob_gene = ();
				foreach my $net (sort keys %{$network{$gene}}){
					my $degree = $network{$gene}{$net};
					my $greater = scalar grep { $_ > $degree } @{$net_degree{$net}};					
					my $prob = sprintf("%.3f", 1- $greater/scalar @{$net_degree{$net}});	
					push @prob_gene, "$net($prob)";
				}
				my $tmp_hub = $gene.":".join('',@prob_gene);
				$self->{HUB} -> {$id} -> {$tmp_hub}=1;
			}
		}
	}
}

## Coding : VAT (variant annotation tool) analysis
sub coding{
	my $self = shift;
	my ($snp_input,$file_interval,$file_fasta,$nc_mode)=@_; 
	my @tmp;
	my %hub;
	my %selection;
	my $gene;
	my @vat_out;
	my $vat_out;
	my $id;
	
# run VAT
	
	if($nc_mode==0){
		
		@vat_out = split /\n/, `awk '{FS="\t";OFS="\t"}BEGIN{print "##fileformat=VCFv4.0\\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"}{print \$1,\$3,\$2,\$4,\$5,".","PASS","."}' $snp_input | indelMapper $file_interval $file_fasta`;
		
		foreach $vat_out(@vat_out){
			if ($vat_out !~ /^#/){
				@tmp =split /\s+/,$vat_out;
				$id = join("\t",@tmp[0,2,1,3,4]);
				$tmp[7] =~ s/insertionFS/frameshift_variant/g;
				$tmp[7] =~ s/insertionNFS/inframe_insertion/g;
				$tmp[7] =~ s/deletionFS/frameshift_variant/g;
				$tmp[7] =~ s/deletionNFS/inframe_deletion/g;
				$tmp[7] =~ s/startOverlap/start_lost/g;
				$tmp[7] =~ s/endOverlap/stop_lost/g;
				$tmp[7] =~ s/spliceOverlap/splice_variant/g;
				
				$self ->{VAT} ->{$id} = $tmp[7];
			}
		}
	}
}

## 2. Motif gain
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
	my @end;
	
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
	
	my $ref_seq;
	my @ref_seq;
	
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
	
	foreach $id(sort keys %{$self->{GENE}}){
		my $switch = 0;
		foreach my $gene(keys %{$self -> {GENE}->{$id}}){
			if (defined $self ->{GENE} ->{$id} ->{$gene} -> {"Promoter"} || defined $self ->{GENE}->{$id}->{$gene}->{"Distal"}){
				$switch = 1;
			}
		}
		if ($switch == 1){
			@des = split /\t+/,$id;
			push @id,$id;
			$chr = $des[0];
			$chr =~ s/chr//g;
			$start = $des[1];
			$end = $des[2];
			$snp_file .= join("",$chr,"\t",$start-29,"\t",$end+29,"\n");
			push @ref,uc($des[3]);
			push @alt,uc($des[4]);
			push @start,$start;
			push @end,$end;
		}
	}
	
	open(O,">$out_tmp")||die;
	print O $snp_file;
	close O;

	if (-s $out_tmp){

	@des = split /\n/, `fastaFromBed -fi $reference_file -bed $out_tmp -fo stdout`;	
	
	if (scalar @des >0){
		for $i (0 .. (scalar @des/2)-1){
			$ref = $ref[$i]; 
			$alt = $alt[$i];
			$alt =~ s/-//g;
			$ref =~ s/-//g;

			$id = $id[$i]; 
			$start = $start[$i]; 
			$end = $end[$i];
			$extract_seq=uc($des[$i*2+1]);
			$ref_seq = $extract_seq;
			
			substr($extract_seq,29,length($ref)) = $alt;
			
			#positive strand
			@extract_seq = split //,$extract_seq;
			@ref_seq = split //,$ref_seq;
			
			&seq_scan(\@extract_seq,\@ref_seq, "+");
		
		#negative strand
			$extract_seq =~ tr/ATCGatcg/TAGCTAGC/;
			@extract_seq = reverse(split //,$extract_seq);
			
			$ref_seq =~ tr/ATCGatcg/TAGCTAGC/;
			@ref_seq = reverse(split //,$ref_seq);
			
			&seq_scan(\@extract_seq,\@ref_seq,"-");
		}
	}
}

# sub-routine for the sequence scanning ... 
	
	sub seq_scan{
		my ($seq,$refseq, $strand) = @_;
		
		my @seq = @{$seq};
		my @refseq = @{$refseq};
		undef my %refpwm;
		
		
		foreach $prev_name(sort keys %motif){
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
			
					
			# sequence scanning .... reference 
			for ($i=30 - $motif_length; $i < 29 +length($ref); $i ++){
				$ref_pwm_score = 0;
				
				# calculate reference & alternative PWM score;
				for ($j = 0; $j < $motif_length; $j++){
					$ref_pwm_score += $motif{$prev_name}->[$j]->{$refseq[$i+$j]};
				}
				
				if (defined $score_upper{$prev_name} && $ref_pwm_score >= $score_upper{$prev_name}){
						$refpwm{$prev_name}=1
				}elsif(defined $score_lower{$prev_name} && $ref_pwm_score < $score_lower{$prev_name}){	
				}else{	
					$ref_p = (split /\s+/,`TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $ref_pwm_score -m pwm/$prev_name -w`)[2];
					if ($ref_p < $p_cut){
						$refpwm{$prev_name}=1;
					}
				}
			}			
		}
		#print keys %refpwm,"\n";
		
		foreach $prev_name(sort keys %motif){
			if (defined $refpwm{$prev_name}){
			
			}else{
			$motif_length = scalar @{$motif{$prev_name}};
			
	
			# sequence scanning....  alternative 		
			for ($i=30 - $motif_length; $i < 29 + length($alt); $i ++){
				$alt_pwm_score = 0;
				
				# calculate reference & alternative PWM score;
				for ($j = 0; $j < $motif_length; $j++){
					$alt_pwm_score += $motif{$prev_name}->[$j]->{$seq[$i+$j]};
				}
				
				if (defined $score_upper{$prev_name} && $alt_pwm_score >= $score_upper{$prev_name}){
						&out_print($prev_name,join('',@seq[$i..$i+$motif_length-1]), $motif_length,$alt_pwm_score,$strand,$id);	
				}elsif(defined $score_lower{$prev_name} && $alt_pwm_score < $score_lower{$prev_name}){	
				}else{	
					$alt_p = (split /\s+/,`TFMpvalue-sc2pv -a 0.25 -c 0.25 -g 0.25 -t 0.25 -s $alt_pwm_score -m pwm/$prev_name -w`)[2];
					if ($alt_p < $p_cut){
						&out_print($prev_name,join('',@seq[$i..$i+$motif_length-1]),$motif_length,$alt_pwm_score,$strand,$id);	
					}
				}
			}					
		}
		}		
	}

#end of routine
	
	sub out_print{
		my ($name,$motif_seq,$motif_length,$alt_pwm_score,$strand,$id) = @_;
		my $tmp = join("", "$name#",$motif_seq,'#',$motif_length,"#$strand#",sprintf("%.3f", $alt_pwm_score));
		if (defined $self ->{MOTIFG}->{$id}){
			$self->{MOTIFG} ->{$id} = join(",",$self->{MOTIFG}->{$id},$tmp);
		}else{
			$self -> {MOTIFG}->{$id} = $tmp;
		}
	}
	
}

## Integration of all outputs
sub intergrate{
	my $self=shift;
	
	$| = 1; #Auto Flush

	my ($output_format,$sample,$gene_dir,$reg_net,$out,$de_data) = @_;
	
	open(OUT,">$out")||die;
	
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
 	

	&print_output($output_format);

	sub print_output{
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
		
		my $id;
    	foreach $id (sort keys %{$self->{GERP}}){
    		print OUT $id,"\t",$sample,"\t";
    		print OUT $self -> {GERP} -> {$id},";";
    		
    		if (defined $self -> {VAT} -> {$id}){
    			print OUT $self -> {VAT} -> {$id},";";
    		}else{
    			print OUT ".;";
    		}

    		if (defined $self -> {HUB} -> {$id}){
    			print OUT join(",",sort keys %{$self -> {HUB} -> {$id}}),";";
    		}else{
    			print OUT ".;";
    		}
    		if (defined $self ->{SELECTION} ->{$id}){
    			print OUT join(",",sort keys %{$self -> {SELECTION} -> {$id}}),";";
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
    		
    		if (defined $self -> {GENE} -> {$id}){
    			my $gene_info = "";
    			foreach my $gene (sort keys %{$self -> {GENE} -> {$id}}){
    			    my $info = "$gene(";
						foreach my $tag (sort keys %{$self -> {GENE} -> {$id}->{$gene}}){
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
    		
    		if (defined $self ->{USER} ->{$id}){
    			print OUT join(",",sort keys %{$self ->{USER}->{$id}}),"\n";
    		}else{
    			print OUT ".\n";
    		}
    	
    	}
	}

	
	sub print_vcf_format{

		my $id;
    	foreach $id (sort keys %{$self->{GERP}}){
    		my @tmp = split /\t+/,$id;
    		print OUT $tmp[0],"\t",$tmp[1]+1,"\t.\t",$tmp[3],"\t",$tmp[4],"\t",".\t.\t";
    		print OUT "SAMPLE=$sample;";
    		print OUT "GERP=",$self -> {GERP}->{$id},";";
    		if (defined $self ->{VAT}->{$id}){
    			print OUT $self -> {VAT} -> {$id},";";
    		}
    		if (defined $self -> {HUB} -> {$id}){
    			print OUT "HUB=",join(",", sort keys %{$self -> {HUB} -> {$id}}),";";
    		}
    		if (defined $self -> {SELECTION} -> {$id}){
    			print OUT "GNEG=",join(",", sort keys %{$self -> {SELECTION} -> {$id}}),";";
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
    		if (defined $self ->{GENE}->{$id}){    				
    			my $gene_info = "";
    			my $info = "";
    			foreach my $gene (sort keys %{$self -> {GENE} -> {$id}}){
    			        $info = $info."$gene(";
						foreach my $tag (sort keys %{$self -> {GENE} -> {$id}->{$gene}}){
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
    		if (defined $self ->{USER} ->{$id}){
    				print OUT "USER_ANNO=",join(",",sort keys %{$self ->{USER}->{$id}}),"\n";
    		}else{
    			print OUT "\n";
    		}
    	}
    
   }

}

## Recurrence Analysis
sub combine{
	my ($file,$format) = @_;
	my $vcf_header = join("\\n",'##fileformat=VCFv4.0','##INFO=<ID=SAMPLE,Number=.,Type=String,Description="Sample id">','##INFO=<ID=VA,Number=.,Type=String,Description="Coding Variant Annotation">','##INFO=<ID=HUB,Number=.,Type=String,Description="Network Hubs, PPI (protein protein interaction network), REG (regulatory network), PHOS (phosphorylation network)...">','##INFO=<ID=GNEG,Number=.,Type=String,Description="Gene Under Negative Selection">','##INFO=<ID=GERP,Number=.,Type=String,Description="Gerp Score">','##INFO=<ID=NCENC,Number=.,Type=String,Description="NonCoding ENCODE Annotation">','##INFO=<ID=HOT,Number=.,Type=String,Description="Highly Occupied Target Region">','##INFO=<ID=MOTIFBR,Number=.,Type=String,Description="Motif Breaking">','##INFO=<ID=MOTIFG,Number=.,Type=String,Description="Motif Gain">','##INFO=<ID=SEN,Number=.,Type=String,Description="In Sensitive Region">','##INFO=<ID=USEN,Number=.,Type=String,Description="In Ultra-Sensitive Region">','##INFO=<ID=UCONS,Number=.,Type=String,Description="In Ultra-Conserved Region">','##INFO=<ID=GENE,Number=.,Type=String,Description="Target Gene (For coding - directly affected genes ; For non-coding - promoter, distal (>10kb from TSS) or Medial (within 10kb) regulatory module)">','##INFO=<ID=CANG,Number=.,Type=String,Description="Prior Gene Information, e.g.[cancer][TF_regulating_known_cancer_gene][up_regulated][actionable]...";','##INFO=<ID=USER_ANNO,Number=.,Type=String,Description="Annotations from user-input";','#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO');

	my $bed_header=join("\t","#chr","start","end","ref","alt","sample","gerp;variant.annotation.cds;network.hub;gene.under.negative.selection;ENCODE.annotated;hot.region;motif.analysis;sensitive;ultra.sensitive;ultra.conserved;target.gene[known_cancer_gene/TF_regulating_known_cancer_gene,differential_expressed_in_cancer,actionable_gene];user.annotations");

	if ($format =~ /bed/i){
		`sed -i '1i $bed_header' $file`;
	}else{
		`sed -i '1i $vcf_header' $file`;
	}
}


1;

