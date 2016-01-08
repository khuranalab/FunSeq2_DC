#!/bin/bash

user_input=
maf=1
genome_mode=1
input_format=
output_format=vcf
gene_list=
expression=
class=
nc_mode=0
weight_mode=1
exp_format=
parallel=5
out_path=out
cancer_type=all
score_cut=1.5
user_anno=/home/fuyao/FunSeq_PCAWG/data_context/user_annotations
recurdb_use=0
sv_length_cut=20

function usage 
{
	echo " 
	
	* Usage : ./run.sh -f file -maf MAF -m <1/2> -len length_cut -inf <bed/vcf> -outf <bed/vcf> -nc -o path -g file -exp file -cls file -exf <rpkm/raw> -p int -cancer cancer_type -s score -uw -ua user_annotations_directory -db
	
	Options : 
		-f		[Required] User Input SNVs File
		-inf	 	[Required] Input format - BED or VCF
		-len            [Optional] Maximum length cutoff for Indel analysis, default = 20. Set to 'inf', if no filter.  
		-maf 		[Optional] Minor Allele Frequency Threshold to filter 1KG SNVs,default = 1 
		-m		[Optional] 1 - Somatic Genome (default); 2 - Germline or Personal Genome
		-outf	 	[Optional] Output format - BED or VCF,default is VCF
		-nc		[Optional] Only do non-coding analysis, no need of VAT (variant annotation tool)
		-o		[Optional] Output path, default is the directory 'out'
		-g		[Optional] gene list, only output variants associated with selected genes. 	
		-exp		[Optional] gene expression matrix
		-cls		[Optional] class file for samples in gene expression matrix
		-exf		[Optional] gene expression format - rpkm / raw
		-p		[Optional] Number of genomes to parallel, default = 5
		-cancer		[Optional] cancer type from recurrence database, default is all of the cancer type
		-uw		[Optional] Use unweighted scoring scheme, defalut is weighted
		-s		[Optional] Score threshold to call non-coding candidates, default = 1.5 for weighted scoring & default = 5 for unweighted scoring
		-ua             [Optional] The directory for user-specific annotations, default will be read from directory 'data/user_annotations'
		-db		[Optional] Use the recurrence database to score variants. Recurrence gets a additional score. 


	* Multiple Genomes with Recurrent Output
		
		Option 1: Separate multiple files by ','
		Example: ./run.sh -f file1,file2,file3,... -maf MAF -m <1/2> -inf <bed/vcf> -outf <bed/vcf> ...
		
		Option 2: Use the 6th column of BED file to specify samples
		Example: ./run.sh -f file -maf MAF -m <1/2> -inf bed -outf <bed/vcf> ...
	
	NOTE: Please make sure you have sufficient memory, at least 3G. 

"
}


## Get inputs
if [ -e $1 ];then
usage
exit
fi


while [ "$1" != "" ]; do
    case $1 in
        -f | --file )           shift
                                user_input=$1
                                ;;
	-maf)			shift
				maf=$1
				;;
	-len)			shift
				sv_length_cut=$1
				;;
	-m | --mode)		shift
				genome_mode=$1
				;;
	-inf)			shift
				input_format=$1
				;;
	-outf)			shift
				output_format=$1
				;;
	-nc)			nc_mode=1
				;;
	-o)			shift
				out_path=$1
				;;
	-g)			shift
				gene_list=$1
				;;
       	-exp)			shift
				expression=$1
				;;
	-cls)			shift
				class=$1
				;;
	-exf)			shift
				exp_format=$1
				;;
	-p)			shift
				parallel=$1
				;;
	-cancer)		shift
				cancer_type=$1
				;;
	-uw)			weight_mode=0
				score_cut=5
				;;
	-s)         		shift
				score_cut=$1
				;;
	-ua)			shift
				user_anno=$1
				;;
	-db)			recurdb_use=1
				;;
	-h | --help )		usage
				exit
				
	esac
        shift
done

## check commands ...
if [[ $nc_mode == 0 ]]
then    
	NEEDED_COMMANDS="bedtools tabix perl snpMapper TFMpvalue-sc2pv awk sed"
else
	NEEDED_COMMANDS="bedtools tabix perl TFMpvalue-sc2pv awk sed"
fi

for cmd in ${NEEDED_COMMANDS} ; do
    if ! command -v ${cmd} &> /dev/null ; then
        echo Please install ${cmd}!
        exit -1
    fi
done

## run programs...
if [[ $user_input != "" &&  $maf != ""  &&  $genome_mode != ""  &&  $input_format != ""  &&  $output_format  != "" && $out_path != "" ]]
then
	if [[ $expression != "" || $class != "" || $exp_format != "" ]]
	then
		if [[ $expression != "" && $class != "" && $exp_format != "" && $gene_list != "" ]]
		then
			perl code/funseq2.pl $user_input $maf $genome_mode $input_format $output_format $nc_mode $out_path $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut $gene_list $expression $class $exp_format	
		elif [[ $expression != "" && $class != "" && $exp_format != "" && $gene_list == "" ]]
		then
			perl code/funseq2.pl $user_input $maf $genome_mode $input_format $output_format $nc_mode $out_path $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut $expression $class $exp_format
		else	
			echo "Please input both expression , class label and expression format data"
	
		fi
	else
		if [[ $gene_list != "" ]]
	 	then
			perl code/funseq2.pl $user_input $maf $genome_mode $input_format $output_format $nc_mode $out_path $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut $gene_list
		else
			perl code/funseq2.pl $user_input $maf $genome_mode $input_format $output_format $nc_mode $out_path $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut
		fi
	fi
else
	usage
fi 
