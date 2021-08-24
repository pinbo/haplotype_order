#!/usr/bin/env Rscript
# download the packages here
# https://github.com/tpook92/HaploBlocker
# install.packages("~/Downloads/RandomFieldsUtils_0.6.6.tar.gz",type="source", repo=NULL)
# install.packages("~/Downloads/HaploBlocker_1.6.06.tar.gz",type="source", repo=NULL)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
  stop("Please give the vcf file name, min_majorblock, window_size, and other options for block_calculation function\n")
}

vcf_file = args[1]
min_majorblock = as.numeric(args[2])
window_size = as.numeric(args[3])
maxmissl = as.numeric(args[4]) # maximum missing values to keep a line
maxmissm = as.numeric(args[5]) # maximum missing values to keep a marker
other_options = args[6]
suppressMessages(library(HaploBlocker))
library(RColorBrewer)
# hap1b = read.delim("haplotype_ordered_genotable.txt", na.strings = "-")
# hap1b = hap1b[-nrow(hap1b),]
# dim(hap1b)
data_file <- data_import(vcf_file, inbred=T)
dhm <- data_file[[1]]
bp_map <- data_file[[2]]
rm(data_file)
# get line names
allLines = readLines(vcf_file, n=100) # read 100 lines to get the headers of vcf
hh = grep("^#CHROM",allLines,value = T) # headers
hh = unlist(strsplit(hh,"\t")) # string to vector
rm(allLines)
# add col and row names for dhm
colnames(dhm)=hh[-c(1:9)]
rownames(dhm)=paste0("bp",bp_map)
# filter lines and SNPs
min.marker = (1 - maxmissm) * nrow(dhm)
dhm = dhm[,colSums(dhm==".")<=min.marker]
dim(dhm)
min.line = (1 - maxmissl) * ncol(dhm)
dhm = dhm[rowSums(dhm==".")<=min.line,]
cat("dhm after filtering\n")
dim(dhm)
# blocklist = block_calculation(dhm=dhm,bp_map=bp_map,min_majorblock=min_majorblock,window_size = window_size)
blocklist <- eval(parse(text = paste("block_calculation(dhm=dhm,bp_map=bp_map,min_majorblock=min_majorblock,window_size = window_size,", other_options, ")")))

# blocklist = block_calculation(dhm=vcf_file,
#   window_size = 5,
#   prefilter=T,
#   inbred = T,
#   min_majorblock=min_majorblock
# )


# start and end of each block
blockinfo = blocklist_startend(blocklist)
# Calculate a block-dataset according to the block library 
newdata = block_matrix_construction(blocklist)
cat("block matrix dimention\n")
dim(newdata)
colnames(newdata) = colnames(dhm)[as.numeric(substring(colnames(newdata),6))]
newdata2 = cbind(Block=rownames(newdata), blockinfo, newdata)
write.table(newdata2, "haplo_block_matrix.txt", sep="\t",row.names = F, quote = F)

pdf(file = "haplo-block-graphs.pdf", width = 10, height = 10)
heatmap(t(newdata),Colv = NA, scale = "none")
plot_block(blocklist)
plot_block(blocklist,type = "bp",xlim=range(bp_map))
blocklist_plot(blocklist)
blocklist_plot_xsize(blocklist)
dev.off()


