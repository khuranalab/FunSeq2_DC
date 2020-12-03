#!/bin/bash

#PLEASE ADD YOUR PATH FOR DATA CONTEXT
data_context=

#PLEASE ADD YOUR PATH FOR YOUR INPUT
user_input=

# PARAMETERS CAN BE MODIFIED BY THE USER.
maf=0
genome_mode=1
input_format=vcf
output_format=vcf
gene_list=
expression=
class=
nc_mode=0
weight_mode=1
exp_format=
parallel=5
out_path=output
cancer_type=all
score_cut=1.5

#PLEASE ADD YOUR PATH FOR ANNOTATION
user_anno=data_context/user_annotations
recurdb_use=0
sv_length_cut=20

function usage
{
	echo "

	* Usage : ./funseq2-docker.sh -f file -maf MAF -m <1/2> -len length_cut -inf <bed/vcf> -outf <bed/vcf> -nc -o path -g file -exp file -cls file -exf <rpkm/raw> -p int -cancer cancer_type -s score -uw -ua user_annotations_directory -db

	Options :
		-f		[Required] User Input SNVs File
		-inf	 	[Required] Input format - BED or VCF
		-len            [Optional] Maximum length cutoff for Indel analysis, default = 20. Set to 'inf', if no filter.
		-maf 		[Optional] Minor Allele Frequency Threshold to filter 1KG SNVs,default = 0
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

if [ "$1" == "-h" ];then
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

input_rel_dir=$(dirname $user_input)
input_dir=`cd "$input_rel_dir"; pwd`
input_name=$(basename $user_input)
container_input="/root/Funseq2_docker_latest/input/$input_name"
data_context_dir=`cd $data_context; pwd`
out_dir_name=$(dirname "$out_path")
out_dir=`cd $out_dir_name; pwd`
out_dir="$out_dir/$(basename $out_path)"

## run programs...
if [[ $user_input != "" &&  $maf != ""  &&  $genome_mode != ""  &&  $input_format != ""  &&  $output_format  != "" && $out_path != "" && $data_context != "" ]]
then
	if [[ $expression != "" || $class != "" || $exp_format != "" ]]
	then
		if [[ $expression != "" && $class != "" && $exp_format != "" && $gene_list != "" ]]
		then
			docker run -v $data_context_dir:/root/Funseq2_docker_latest/data_context/ -v $input_dir:/root/Funseq2_docker_latest/input/ -v $out_dir:/root/Funseq2_docker_latest/output/ funseq2-docker-latest perl /root/Funseq2_docker_latest/code/funseq2.pl $container_input $maf $genome_mode $input_format $output_format $nc_mode /root/Funseq2_docker_latest/output/ $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut $gene_list $expression $class $exp_format
		elif [[ $expression != "" && $class != "" && $exp_format != "" && $gene_list == "" ]]
		then
			docker run -v $data_context_dir:/root/Funseq2_docker_latest/data_context/ -v $input_dir:/root/Funseq2_docker_latest/input/ -v $out_dir:/root/Funseq2_docker_latest/output/ funseq2-docker-latest perl /root/Funseq2_docker_latest/code/funseq2.pl $container_input $maf $genome_mode $input_format $output_format $nc_mode /root/Funseq2_docker_latest/output/ $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut $expression $class $exp_format
		else
			echo "Please input both expression , class label and expression format data"

		fi
	else
		if [[ $gene_list != "" ]]
	 	then
			docker run -v $data_context_dir:/root/Funseq2_docker_latest/data_context/ -v $input_dir:/root/Funseq2_docker_latest/input/ -v $out_dir:/root/Funseq2_docker_latest/output/ funseq2-docker-latest perl /root/Funseq2_docker_latest/code/funseq2.pl $container_input $maf $genome_mode $input_format $output_format $nc_mode /root/Funseq2_docker_latest/output/ $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut $gene_list
		else
			docker run -v $data_context_dir:/root/Funseq2_docker_latest/data_context/ -v $input_dir:/root/Funseq2_docker_latest/input/ -v $out_dir:/root/Funseq2_docker_latest/output/ funseq2-docker-latest perl /root/Funseq2_docker_latest/code/funseq2.pl $container_input $maf $genome_mode $input_format $output_format $nc_mode /root/Funseq2_docker_latest/output/ $parallel $cancer_type $score_cut $weight_mode $user_anno $recurdb_use $sv_length_cut
		fi
	fi
else
	usage
fi
