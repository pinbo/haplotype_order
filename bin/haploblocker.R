#!/usr/bin/env Rscript
# download the packages here
# https://github.com/tpook92/HaploBlocker
# install.packages("~/Downloads/RandomFieldsUtils_0.6.6.tar.gz",type="source", repo=NULL)
# install.packages("~/Downloads/HaploBlocker_1.6.06.tar.gz",type="source", repo=NULL)

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 1) {
  min_majorblock = 5000
} else {
  min_majorblock = as.numeric(args[1])
}

suppressMessages(library(HaploBlocker))
library(RColorBrewer)
hap1b = read.delim("haplotype_ordered_genotable.txt", na.strings = "-")
hap1b = hap1b[-nrow(hap1b),]
dim(hap1b)
blocklist = block_calculation(dhm=hap1b[,-c(1:4)],bp_map=hap1b[,2],min_majorblock=min_majorblock,window_size = 5,prefilter=T)
# blocklist = block_calculation(dhm="hap1B-filtered-lines.vcf")#,
                                 # window_size = 5,
                                 # prefilter=T
                                 # )
pdf(file = "haplo-block-graphs.pdf")
plot_block(blocklist)
plot_block(blocklist,type = "bp",xlim=range(hap1b[,2]))
blocklist_plot(blocklist)
blocklist_plot_xsize(blocklist)
dev.off()

# start and end of each block
blockinfo = blocklist_startend(blocklist)
# Calculate a block-dataset according to the block library 
newdata = block_matrix_construction(blocklist)
newdata2 = cbind(Block=rownames(newdata), blockinfo, newdata)
write.table(newdata2, "haplo_block_matrix.txt", sep="\t",row.names = F, quote = F)




