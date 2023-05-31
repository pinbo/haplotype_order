#!/bin/bash

echo "The script you are running has basename `basename "$0"`, dirname `dirname "$0"`"
echo "The present working directory is `pwd`"


vcf_file=$1 # vcf
#meta_data=$5 # metadata for lines to merge to SNP table
maxmissl=$2
maxmissm=$3

cp $vcf_file input.vcf

script_path="`dirname "$0"`/bin" # scripts folder
echo $script_path

#step 2: convert from GT to SNP alleles
cmd2="$script_path/convert_vcf_calls_to_SNP.py input.vcf out.txt"
echo "Step 2: Convert vcf to SNP table: $cmd2"
eval $cmd2
# Step 3: sort lines by similarity for haplotype analysis
cmd3="$script_path/haplotype_analysis_no_metadata.R out.txt $maxmissl $maxmissm"
echo "Step 3: Get the sorted SNP table for haplotype analysis: $cmd3"
eval $cmd3

echo "Finish successfully! Output data is 'haplotype_ordered_genotable.txt' "