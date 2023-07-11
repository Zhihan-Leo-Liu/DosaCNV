#!/bin/bash

# Initialize variables.
inputfile=""
max_n_gene=100

# Parse command line options.
while [[ $# -gt 0 ]]
do
    key="$1"
    case $key in
        -i|--input)
        inputfile="$2"
        shift # past argument
        shift # past value
        ;;
        -m|--max_n_gene)
        max_n_gene="$2"
        shift # past argument
        shift # past value
        ;;
        *)
        # unknown option
        echo "Error: Unknown option $1"
        exit 1
        ;;
    esac
done

base_filename=$(basename "$inputfile" .bed)
output_filename="${base_filename}_index.txt"

# Processing
if awk 'NR==1 {exit !(NF==5)}' $inputfile
then
    bedtools intersect -wa -wb -a $inputfile -b ./data/CDS_hg19_basic.bed|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' > test_del_gene.bed
    echo 'variant_id,gene_id,label' > intermediate.csv
    sort -k4,4 -k9,9 test_del_gene.bed | bedtools groupby -g 4,9 -c 5 -o first | awk 'BEGIN{OFS=","}{print $1,$2,$3}' >> intermediate.csv
    awk -F, -v limit=$max_n_gene 'NR==1{print;next}{lines[NR]=$0;variant[$1]++;line_to_variant[NR]=$1} END{for(i in lines){if(variant[line_to_variant[i]]<=limit) print lines[i]}}' intermediate.csv > intermediate1.csv
    awk 'BEGIN{FS=",";OFS=","}{if(NR==1) print "variant_id,gene_id,label,type";else print $1,$2,$3,"DEL"}' intermediate1.csv > $output_filename
else
    bedtools intersect -wa -wb -a $inputfile -b ./data/CDS_hg19_basic.bed|awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8}' > test_del_gene.bed
    echo 'variant_id,gene_id' > intermediate.csv
    sort -k4,4 -k8,8 test_del_gene.bed | bedtools groupby -g 4,8 -c 8 -o first | awk 'BEGIN{OFS=","}{print $1,$2}' >> intermediate.csv
    awk -F, -v limit=$max_n_gene 'NR==1{print;next}{lines[NR]=$0;variant[$1]++;line_to_variant[NR]=$1} END{for(i in lines){if(variant[line_to_variant[i]]<=limit) print lines[i]}}' intermediate.csv > intermediate1.csv
    awk 'BEGIN{FS=",";OFS=","}{if(NR==1) print "variant_id,gene_id,type";else print $1,$2,"DEL"}' intermediate1.csv > $output_filename
fi

rm intermediate.csv
rm test_del_gene.bed
rm intermediate1.csv
