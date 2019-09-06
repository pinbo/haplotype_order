#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please give the marker data and metadata file\n")
}

snpfile = args[1]
metadata = args[2]
maxmiss = as.numeric(args[3]) # maximum missing values to keep a line

dd= read.delim(snpfile, check.names = F, as.is = T, na.strings = "N")
md = read.delim(metadata)
#md[1:5,1:4]
rownames(md) = md[,1]
dim(dd)

nc = ncol(dd)
nr = nrow(dd)

snp = t(dd[,6:nc])
#snp[1:5,1:5]
snpname = apply( dd[ , 1:2 ] , 1 , paste , collapse = "-" )

dd2 = sapply(1:nr, function(x) {
  tt = factor(dd[x, 6:nc], levels = c(dd[x,4], "H", dd[x,5]))
  tt2 = as.numeric(tt)
})

#dd2[1:20,1:6]

dimnames(dd2) = list(colnames(dd)[6:nc], snpname)
dimnames(snp) = list(colnames(dd)[6:nc], snpname)

dd3 = dd2 - 1

# to see the maximum missing values I can keep.
for (i in rev(seq(0.1, 0.9, 0.05))){
  dd4 = dd3[rowSums(is.na(dd3)) < ncol(dd3)*i, ]
  snp.dist <- dist(dd4)
  if (sum(is.na(snp.dist)) == 0) break
}
i # the final missing proportion
maxmiss = ifelse(maxmiss < i, maxmiss, i)
dd4 = dd3[rowSums(is.na(dd3)) < ncol(dd3)*maxmiss, ]
dim(dd4)

##
snp.dist <- dist(dd4)
fit <- hclust(snp.dist,method="ward.D2")
str(fit)
pdf(file="dendrogram_of_lines.pdf", width = 10)
plot(fit, hang=-1, labels=FALSE)
# rect.hclust(fit, h = 10, border = "red")
dev.off()
ordered_line_names = fit$labels[fit$order]
snp2 = snp[ordered_line_names,]
md2 = md[ordered_line_names,]
dd5 = dd4[ordered_line_names,]
all(rownames(snp2) == rownames(md2))
outdata = cbind(md2, snp2)
outdata2 = cbind(md2, dd5)

outdata3 = cbind(md[rownames(snp),], snp)
#outdata[1:4,1:30]

### write all output data
write.table(outdata, "haplotype_ordered_genotable.txt", row.names = F, sep="\t", quote=F, na = "-")

write.table(outdata2, "haplotype_ordered_genotable_numbers.txt", row.names = F, sep="\t", quote=F, na = "-")

write.table(outdata3, "raw_extracted_genotable.txt", row.names = F, sep="\t", quote=F, na = "-")
