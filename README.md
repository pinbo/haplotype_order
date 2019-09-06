# Haplotype order
Extract variations in a chromosome region and order the lines for haplotype analysis. I use it to extract SNPs and do haplotype analysis from the [1000 wheat exomes project](http://wheatgenomics.plantpath.ksu.edu/1000EC/). If you only need to extract SNPs, you can also download them form [T3 website](https://triticeaetoolbox.org/wheat/genotyping/download-vcf.php).

**Haplotype_order includes three parts:**

1. Extract a chromosome region from an indexed vcf file with `tabix`.

2. Convert the subsetted vcf file to a hapmap-like genotyping table (ATGC).

3. Do a cluster analysis with euclidean distance in R to group the individuals as haplotypes. So you do not need to cut and paste by eyes.

## Input parameters

The input parameters are position-based:

1. chromosome to subset
1. start position
1. stop position
1. vcf to subset. Need to run 'tabix -C xxxx.vcf.gz' to index the vcf file.
1. the metadata of line information for the output SNP table
1. maximum missing proportion for lines (0-1). I use 0.5 to filter out lines with missing values >50%.

## Output in Galaxy

Four files are output from the analysis:

1. A pdf file with dendrogram graph of lines
1. A raw table with line information and A/T/G/C coded genotyping results
1. A haplotype ordered table with line information and A/T/G/C coded genotyping results
1. A haplotype ordered table with line information and 0/1/2 coded genotyping results: 0 is homozygous reference allele, 2 is homozygous alternative allele, and 1 is heterozygous.
