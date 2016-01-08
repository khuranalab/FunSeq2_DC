#!/usr/bin/perl

use strict;
use IO::Handle;
use Parallel::ForkManager;
use code::Funseq_SNV;
use code::Funseq_Indel;

##   Obtain basic input parameters ... 
$| = 1; 
my $infile = $ARGV[0];        # input variants file 
my $maf = $ARGV[1];           # minor allele frequency threshold
my $genome_mode = $ARGV[2];   # 1: somatic ; 2: germline | personal genome
my $informat = $ARGV[3];      # input format [bed | vcf]
my $outformat = $ARGV[4];     # output format [bed| vcf] 
my $nc_mode=$ARGV[5];         # only do non-coding if nc_mode=1
my $output_path = $ARGV[6];   # output path
my $num_per_run = $ARGV[7];   # number of genome per run.
my $cancer_type = $ARGV[8];   # cancer type retrieved.
my $score_cut = $ARGV[9];	  # Non-coding candidate score cut 
my $weight_mode = $ARGV[10];  # weighted / unweighted scoring scheme. 0: unweighted; 1: weighted. 
my $user_anno_dir = $ARGV[11];  #directory for user-specific annotations
my $recur_db_use = $ARGV[12];  
my $sv_length_cut = $ARGV[13];
my $gene_list = "";	          # gene list file
my $expression = "";	      # expression file
my $class = "";		          # sample class file
my $exp_format = "";		  # expression format : rpkm / raw read count
my $motif_p_value_cut = 4e-8;


##   Get additional parameters ...
if (scalar @ARGV == 15){
	$gene_list = $ARGV[14];
}
if (scalar @ARGV==17){
	$expression = $ARGV[14];
	$class = $ARGV[15];
	$exp_format = $ARGV[16];
}
if (scalar @ARGV==18){
	$gene_list = $ARGV[14];
	$expression = $ARGV[15];
	$class = $ARGV[16];
	$exp_format = $ARGV[17];
}


##   Check existence of Ouput directory 
unless (-d $output_path){
	mkdir $output_path;
}


##   Input parameter checking 
die "Error. Please enter proper MAF...\n" unless ($maf =~ /\d+/ && $maf >=0 && $maf <= 1);
die "Error. Please specify correct Genome Mode...\n" unless ($genome_mode == 1 || $genome_mode ==2); 
die "Error. Input format should be bed or vcf...\n" unless ($informat =~ /bed|vcf/i);
die "Error. Output format should be bed or vcf...\n" unless ($outformat =~ /bed|vcf/i);
die "Error. Indel length cutoff should be integer or 'inf'...\n" unless($sv_length_cut =~ /^\d+$/ || $sv_length_cut eq 'inf');

if ($exp_format ne ""){
	die "Error. Expression format should be rpkm or raw...\n" unless ($exp_format =~ /rpkm|raw/i);
}


#################### read parameters ###########
##   Required files
my %variable;
my $file_path;

open(PA,"config.txt")||die;
while(<PA>){
	if (/^file_path=(.+)$/){
		$file_path=$1;
	}elsif (/^(\w+) *= *(.+)$/) {
		$variable{$1}=$file_path.'/'.$2;
  	}
}


my $tgp_snp=$variable{'tgp_snp'} if defined $variable{'tgp_snp'};					# 1000 genomes SNP file
my $gerp_file=$variable{'gerp_file'} if defined $variable{'gerp_file'};										# Gerp score file 
my $sensitive=$variable{'sensitive'} if defined $variable{'sensitive'};										# Sensitive Region
my $conserved=$variable{'conserved'} if defined $variable{'conserved'};                           			# Ultra conserved region

# Annotations (e.g.ENCODE )
my $encode_annotation=$variable{'encode_annotation'} if defined $variable{'encode_annotation'}; 								# ENCODE  Annotation (except Motif info)
my $bound_motif=$variable{'bound_motif'} if defined $variable{'bound_motif'};  									# ENCODE TF bound Motif 
my $hot_file=$variable{'hot_file'} if defined $variable{'hot_file'};										# Highly occupied regions
my $enhancer=$variable{'enhancer'} if defined $variable{'enhancer'};		 								# Enhancer-gene pairs; 
my $motif_pfm=$variable{'motif_pfm'} if defined $variable{'motif_pfm'};  									# Motif PFM files 
my $score_file=$variable{'score_file'} if defined $variable{'score_file'};  									# Rough motif score (corresponding to ~4e-8) produced by TFM-Pvalue 

# Networks
my $network_dir=$variable{'network_dir'} if defined $variable{'network_dir'};									# Network hubs
my $reg_net=$variable{'reg_net'} if defined $variable{'reg_net'};										# Regulatory Network
my $gencode=$variable{'gencode'} if defined $variable{'gencode'};

# Genome sequences
my $reference_file=$variable{'reference_file'} if defined $variable{'reference_file'}; 								# Reference genome
my $ancestral_file=$variable{'ancestral_file'} if defined $variable{'ancestral_file'};									# hg19 ancestral allele

# Gene info & recurrence
my $gene_info_dir=$variable{'gene_info_dir'} if defined $variable{'gene_info_dir'};                                  # Prior knowledge of genes, such as cancer / actionable genes.
my $cancer_dir=$variable{'cancer_dir'} if defined $variable{'cancer_dir'};										# Public Recurrent data
my $selection=$variable{'selection'} if defined $variable{'selection'};  								    # Gene under negative selection 

# weighted scoring scheme
my $weight_file=$variable{'weight_file'} if defined $variable{'weight_file'};                                    # weighted scoring scheme

# GENCODE
my $gencode_v = (split /\./,(split /\//,`ls $gencode/*.promoter.bed`)[-1])[1]; 
my $cds = "$file_path/gencode/gencode.$gencode_v.cds.bed";  							# GENCODE CDS
my $promoter = "$file_path/gencode/gencode.$gencode_v.promoter.bed"; 					# Promoter Region; upstream -2.5kb
my $intron = "$file_path/gencode/gencode.$gencode_v.intron.bed";                        # Intronic Region
my $utr = "$file_path/gencode/gencode.$gencode_v.utr.bed";                              # UTR Region
my $coding_interval = "$file_path/gencode/gencode.$gencode_v.cds.interval";				# VAT
my $coding_fasta = "$file_path/gencode/gencode.$gencode_v.cds.fa"; 	 					# VAT


########   Check required files
die "Error: $tgp_snp not found or empty...\n" unless (-f $tgp_snp && -s $tgp_snp);
die "Error: $encode_annotation not found or empty...\n" unless (-f $encode_annotation && -s $encode_annotation);
die "Error: $cds not found or empty...\n" unless (-f $cds && -s $cds);
die "Error: $bound_motif not found or empty...\n" unless (-f $bound_motif && -s $bound_motif);
die "Error: $promoter not found or empty...\n" unless (-f $promoter && -s $promoter);
die "Error: $enhancer not found or empty...\n" unless (-f $enhancer && -s $enhancer);
if ($genome_mode==2){
	die "Error: $ancestral_file not found or empty...\n" unless (-f $ancestral_file && -s $ancestral_file);
}
if ($nc_mode==0){
	die "Error: $coding_interval not found or empty...\n" unless (-f $coding_interval && -s $coding_interval);
	die "Error: $coding_fasta not found or empty...\n" unless (-f $coding_fasta && -s $coding_fasta);
}
die "Error: $motif_pfm not found or empty...\n" unless (-f $motif_pfm && -s $motif_pfm);
die "Error:	$reference_file not found or empty...\n" unless (-f $reference_file && -s $reference_file);
die "Error: $score_file not found or empty...\n" unless (-f $score_file && -s $score_file);
die "Error: $weight_file not found or empty...\n" unless (-f $weight_file && -s $weight_file);

##   Main output files
my $file_detail = "$output_path/Output.$outformat";
my $file_recur = "$output_path/Recur.Summary";
my $file_driver = "$output_path/Candidates.Summary";
my $file_err = "$output_path/Error.log";
my $file_indel = "$output_path/Output.indel.$outformat";

unlink("$file_detail");
unlink("$file_recur");
unlink("$file_driver");
unlink("$file_err");
unlink("$file_indel");

## 	 Expression analysis
if (-f "$output_path/DE.gene.txt"){
	unlink "$output_path/DE.gene.txt";
}
if (-f "$output_path/DE.pdf"){
	unlink "$output_path/DE.pdf";
}
if ($expression ne "" && $class ne ""){
	print "Differential Gene Expression Analysis ...\n";
	`Rscript scripts/differential_gene_expression.r $expression $class $exp_format $output_path $file_err`;
}


##   If number of input files larger than 1, then do parallel analysis. 
my @input = split /,/,$infile;
my @samples;

if($informat =~ /bed/i && scalar @input ==1){
	@samples = split /\n/, `cut -f 6 $infile | sort|uniq`;
	if (scalar @samples <=1){
		&main($input[0]);
		my $tag = (split /\./,((split /\//,$input[0])[-1]))[0];
		my $tmp_err = "$output_path/$tag.err";
		if (-f $tmp_err){
			`mv $tmp_err $file_err`; 
		}
		my $tmp_in = "$output_path/$tag.result.$outformat";
		my $tmp_indel = "$output_path/$tag.indel.$outformat";

		if (-f $tmp_in){
			`mv $tmp_in $file_detail`;
		}		
		
		if (-f $tmp_indel){
			`mv $tmp_indel $file_indel`;
		}
	}else{
    	my $pm = new Parallel::ForkManager($num_per_run);
    	my $data;
    	my $i;
    	my $j;
    	my $num = int(scalar @samples / $num_per_run);
		
		for $j (1 .. $num){
			my @sample_tag = ();
			for $i (($j-1)*$num_per_run  .. $j*$num_per_run - 1){
       			sleep(3);     			# to avoid memory overuse 
       			my $sample_id = $samples[$i];
       			my $tag = (split /\./,((split /\//,$sample_id)[-1]))[0];
				push @sample_tag,$tag;
				
       			my $pid = $pm -> start and next;
       			my $sample_file = "$output_path/$sample_id";
				`awk '{FS="\t";OFS="\t"} \$6 == "$sample_id"' $infile > $sample_file`;
				&main($sample_file);
				unlink($sample_file);
        		$SIG{INT} = sub {kill 9, $pid;};
				$pm -> finish;
			}
			$pm -> wait_all_children;
			
			foreach my $tag(@sample_tag){
				my $tmp_in = "$output_path/$tag.result.$outformat";
				my $tmp_err = "$output_path/$tag.err";
				my $tmp_indel = "$output_path/$tag.indel.$outformat";
				
				if (-f $tmp_err){
					`sed 's/^/$tag\t/' $tmp_err >> $file_err`; 
				}
				if (-f $tmp_in){
					`cat $tmp_in >> $file_detail`;
				}
				if (-f $tmp_indel){
					`cat $tmp_indel >> $file_indel`;
				}
				unlink("$tmp_in");
				unlink("$tmp_err");
				unlink("$tmp_indel");
			}
		}
	
		if ((scalar @samples - $num*$num_per_run )>0){
			my @sample_tag = ();
			for $i ($num*$num_per_run .. $#samples){
				sleep(3);
				my $sample_id = $samples[$i];
				my $tag = (split /\./,((split /\//,$sample_id)[-1]))[0];
				push @sample_tag,$tag;
			
       			my $pid = $pm -> start and next;
       			my $sample_file = "$output_path/$sample_id";
				`awk '{FS="\t";OFS="\t"} \$6 == "$sample_id"' $infile > $sample_file`;
        		&main($sample_file);
        		unlink($sample_file);
        		$SIG{INT} = sub {kill 9, $pid;};
				$pm -> finish;
			}
			$pm -> wait_all_children;

			foreach my $tag(@sample_tag){
				my $tmp_in = "$output_path/$tag.result.$outformat";
				my $tmp_err = "$output_path/$tag.err";
				my $tmp_indel = "$output_path/$tag.indel.$outformat";

				if (-f $tmp_err){
					`sed 's/^/$tag\t/' $tmp_err >> $file_err`; 
				}
				if (-f $tmp_in){
					`cat $tmp_in >> $file_detail`;
				}
				
				if (-f $tmp_indel){
					`cat $tmp_indel >> $file_indel`;
				}
				unlink("$tmp_in");
				unlink("$tmp_err");
				unlink("$tmp_indel");

			}
			
		}
	}
	Funseq_SNV::recur($file_detail,$file_recur,$file_driver,$outformat,$cancer_dir,$cancer_type,$score_cut,$user_anno_dir,$weight_file,$recur_db_use);   # last input is the non-coding score to define drivers
		
	if (-s $file_indel){
		Funseq_Indel::combine($file_indel,$outformat);
	}

}else{
	if (scalar @input ==1){
		&main($input[0]);
		my $tag = (split /\./,((split /\//,$input[0])[-1]))[0];
		my $tmp_err = "$output_path/$tag.err";
		if (-f $tmp_err){
			`mv $tmp_err $file_err`; 
		}
		my $tmp_in = "$output_path/$tag.result.$outformat";
		if (-f $tmp_in){
			`mv $tmp_in $file_detail`;
		}
		unlink($tmp_err);
		
	}else{
	
    	my $pm = new Parallel::ForkManager($num_per_run);
    	my $data;
    	my $i;
    	my $j;
    	my $num = int(scalar @input / $num_per_run);

		for $j (1 .. $num){
			my @sample_tag = ();
			for $i (($j-1)*$num_per_run  .. $j*$num_per_run - 1){
       			sleep(3);     			# to avoid memory overuse 
       			
       			my $tag = (split /\./,((split /\//,$input[$i])[-1]))[0];
				push @sample_tag,$tag;
       			
       			my $pid = $pm -> start and next;
				&main($input[$i]);
				$SIG{INT} = sub {kill 9, $pid;};
				$pm -> finish;
			}
			$pm -> wait_all_children;
			
			foreach my $tag(@sample_tag){
				my $tmp_in = "$output_path/$tag.result.$outformat";
				my $tmp_err = "$output_path/$tag.err";
				if (-f $tmp_err){
					`sed 's/^/$tag\t/' $tmp_err >> $file_err`; 
				}
				if (-f $tmp_in){
					`cat $tmp_in >> $file_detail`;
				}
				unlink("$tmp_in");
				unlink("$tmp_err");
			}
		}
	
		if ((scalar @input - $num*$num_per_run )>0){
			my @sample_tag=();
			
			for $i ($num*$num_per_run .. $#input){
				sleep(3);
				my $tag = (split /\./,((split /\//,$input[$i])[-1]))[0];
				push @sample_tag,$tag;
				
       			my $pid = $pm -> start and next;
        		&main($input[$i]);
        		$SIG{INT} = sub {kill 9, $pid;};
				$pm -> finish;
			}
			$pm -> wait_all_children;
			foreach my $tag(@sample_tag){
				my $tmp_in = "$output_path/$tag.result.$outformat";
				my $tmp_err = "$output_path/$tag.err";
				if (-f $tmp_err){
					`sed 's/^/$tag\t/' $tmp_err >> $file_err`; 
				}
				if (-f $tmp_in){
					`cat $tmp_in >> $file_detail`;
				}
				unlink("$tmp_in");
				unlink("$tmp_err");
			}
			
		}
	}
	Funseq_SNV::recur($file_detail,$file_recur,$file_driver,$outformat,$cancer_dir,$cancer_type,$score_cut,$user_anno_dir,$weight_file,$recur_db_use);    # last input is the non-coding score to define drivers
}





##   Pipeline ## 
sub main{
	my ($infile) = @_;
	die "Error: input file not found or empty...\n" unless (-f $infile && -s $infile);
	my $tag = (split /\./,(split /\//,$infile)[-1])[0];
	if (glob("$output_path/$tag.*")){
		`rm $output_path/$tag.*`;
	}
	
	# Record Error. 
	my $err = "$output_path/$tag.err";
	open(ERR,">$err")||die;
	STDERR->fdopen( \*ERR,  'w' ) or die $!;
	
	# Intermedia files. Delete after finishing ...
	my $out_snv_filter= "$output_path/$tag.out.snv.filter";  				# variants passed 1kg filter
	my $out_nc = "$output_path/$tag.out.nc"; 								# non-coding variants
	my $out_cds = "$output_path/$tag.out.cds"; 							  	# coding variants
	my $out_vat = "$output_path/$tag.out.vat"; 								# vat analysis 
	my $out_motif = "$output_path/$tag.out.motif";                          # motif gain analysis
	my $out_indel = "$output_path/$tag.out.indel";  						# indels 
	
	
	# Output File 
	my $run_out = "$output_path/$tag.result.$outformat";
 	my $status;
	
	# Get variants associated with certain genes... 
	if ($gene_list ne ""){
		my $out_gene_sel = "$output_path/$tag.out.gene.sel";
		`awk '{print "\t"\$_"\$"}' $gene_list| grep -f - $cds | cut -f 1,2,3 > $out_gene_sel`;
		`awk '{print "\t"\$_"\$"}' $gene_list| grep -f - $promoter | cut -f 1,2,3 >> $out_gene_sel`;
		`awk '{print "\t"\$_"\$"}' $gene_list| grep -f - $intron | cut -f 1,2,3 >> $out_gene_sel`;
		`awk '{print "\t"\$_"\$"}' $gene_list| grep -f - $utr | cut -f 1,2,3 >> $out_gene_sel`;
		`awk '{print "\t"\$_"\$"}' $gene_list | grep -f - $enhancer | cut -f 1,2,3 >> $out_gene_sel`;
		`intersectBed -u -a $infile -b $out_gene_sel > $output_path/$tag.infile`;
		$infile = "$output_path/$tag.infile";
		die "Error: no variants occurred in / associated with requested genes...\n" unless (-f $infile && -s $infile);
	}
 	
 	
 	### ------------------------   Pipeline  ------------------------------------ ###

 	print "------------- Running: $tag starts (0%)---------------\n";
	
	my $de_data = "$output_path/DE.gene.txt";

 	# check for input file format 
    print "... Input format check : $informat ...\n";
    $status = Funseq_SNV::format_check($infile,$informat,\*ERR);
	if ($status ==0){
		print "... Format ok ...\n";
	}else{
		exit;
	}
	
	
	# * Step1: 1000 genomes SNV filtering 
	print "... Start filtering SNVs with minor allele frequency = $maf ...\n";
	
	my $data = new Funseq_SNV;
	$data -> snv_filter($infile,$informat,$tgp_snp,$maf,$out_snv_filter,$out_indel,$sv_length_cut);     # filter variants
	
	if (-e "$output_path/$tag.infile"){
		unlink ("$output_path/$tag.infile");
	}
	
	if (-z $out_snv_filter){
		print "Warning: sample $tag - no SNVs left after filtering against natrual variations ...\n";
		print ERR "Warning: no SNVs left after filtering against natrual variations ...\n";
	}else{
	
		$data -> gerp_score($gerp_file,$out_snv_filter);     # Obtain GERP++ scores
	
		# * Step2: Get Annotated (Coding & Non-coding) SNVs
		
    	print "... Annotate variants (30%) ...\n";
    	`intersectBed -u -a $out_snv_filter -b $cds > $out_cds`; # Obtain coding variants 
    	if (-s $out_cds){
        	    `intersectBed -v -a $out_snv_filter -b $out_cds > $out_nc`;
    	}else{
        	    `cp $out_snv_filter $out_nc`;
    	}
    	
        $data -> read_encode($out_snv_filter,$encode_annotation);   # ENCODE annotations for non-coding variants
		
		$data -> user_annotation($out_snv_filter,$user_anno_dir);
	
	
		# * Step 3: Non-coding 
		print "... Start Non-coding SNVs analysis (50%)...\n";
	
		$data -> conserved($out_snv_filter,$conserved);        # Conserved regions
		$data -> sensitive($out_snv_filter,$sensitive);	           # Sensitive regions
		$data -> hot_region($out_snv_filter,$hot_file);                           # non-coding variants in highly occupied regions or not 
		$data -> motif_break($out_snv_filter, $ancestral_file, $bound_motif, $genome_mode, $motif_pfm);      # Motif disruption analysis
		$data -> gene_link($out_snv_filter,$promoter,$enhancer,$intron,$utr,$network_dir);         # Associate non-coding variants to genes & network analysis 
			# Motif gain 
		$data -> motif_gain($motif_pfm,$reference_file,$score_file,$motif_p_value_cut,$out_motif);         # Gain-of-motif analysis. This step should be done after gene_link
	
	

		# * Step 4: Coding 
		print "... Start Coding SNVs analysis (80%)...\n";
	
		$status = Funseq_SNV::coding($out_snv_filter,$coding_interval,$coding_fasta,$network_dir,$selection,$out_vat,$output_path,$nc_mode,$cds);   # Coding analysis
		
		if ($status == 1){
			print "Warning: coding variants not analyzed by VAT ... \n";
			print ERR "Warning: coding variants not analyzed by VAT ... \n";
		}
	
	
	
		# * Step 5: Scoring SNVs 
		$data -> intergrate($outformat,$tag,$out_nc,$out_vat,$gene_info_dir,$reg_net,$run_out,$de_data,$weight_mode,$weight_file);
	}
	
	# * Indel analysis : if any
    if (-s $out_indel){
        print "... Start Indels analysis ...\n";
        
    	my $indel = new Funseq_Indel;
    	my $run_indel_out = "$output_path/$tag.indel.$outformat";
    	my $out_indel_motif = "$output_path/$tag.out.indel.motif";
    
		$indel -> gerp_score($gerp_file,$out_indel);
		$indel -> annotations($out_indel,$conserved,$hot_file,$sensitive,$encode_annotation,$user_anno_dir,$bound_motif);
		$indel -> gene_link($out_indel,$promoter,$enhancer,$intron,$utr,$network_dir,$cds,$selection);
		$indel -> coding($out_indel,$coding_interval,$coding_fasta,$nc_mode);
		$indel -> motif_gain($motif_pfm,$reference_file,$score_file,$motif_p_value_cut,$out_indel_motif); 
		$indel -> intergrate($outformat,$tag,$gene_info_dir,$reg_net,$run_indel_out,$de_data);
	}
    
    print "------------- $tag Done (100%)---------------\n";
    

	# delete all intermedia results 
	my $choice = 1;   # 1 will delete all intermedia results; 0 is used to test scripts. 
	if ($choice == 1){
		`rm $output_path/$tag.out*`;
	}
	
	close ERR;
	
}   

          



