#!/bin/bash

echo "The script you are running has basename `basename "$0"`, dirname `dirname "$0"`"
echo "The present working directory is `pwd`"

chrom=$1
from_bp=$2
to_bp=$3
vcf_file=$4 # vcf.gz
meta_data=$5 # metadata for lines to merge to SNP table
maxmiss=$6

script_path="`dirname "$0"`/bin" # scripts folder
echo $script_path

# step 1: extract SNPs
#cmd1="zcat < $vcf_file | gawk '\$1 == \"$chrom\" &&  \$2 > $to_bp {exit} /^#CHROM/ || (\$1 == \"$chrom\" && \$2 >= $from_bp && \$2 <= $to_bp)' > test.txt"
cmd1="tabix -h $vcf_file $chrom:$from_bp-$to_bp > subset.vcf"
echo "Step 1: Extract SNPs: $cmd1"
eval $cmd1
#step 2: convert from GT to SNP alleles
cmd2="$script_path/convert_vcf_calls_to_SNP.py subset.vcf out.txt"
echo "Step 2: Convert vcf to SNP table: $cmd2"
eval $cmd2
# Step 3: sort lines by similarity for haplotype analysis
cmd3="$script_path/haplotype_analysis.R out.txt $meta_data $maxmiss"
echo "Step 3: Get the sorted SNP table for haplotype analysis: $cmd3"
eval $cmd3

echo "Finish successfully! Output data is 'haplotype_ordered_genotable.txt' "
