#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 4) {
  stop("Please give the marker data, metadata file, maximum missing data for lines, and maximum missing data for markers\n")
}

snpfile = args[1]
metadata = args[2]
maxmissl = as.numeric(args[3]) # maximum missing values to keep a line
maxmissm = as.numeric(args[4]) # maximum missing values to keep a line

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


dimnames(dd2) = list(colnames(dd)[6:nc], snpname)
dimnames(snp) = list(colnames(dd)[6:nc], snpname)

dd2.1 = dd2 - 1

# remove markers with too many missing values
dd3 = dd2.1[ , colSums(is.na(dd2.1)) < nrow(dd2.1)*maxmissm ]

# to see the maximum missing values I can keep.
for (i in rev(seq(0.1, 0.9, 0.05))){
  dd4 = dd3[rowSums(is.na(dd3)) < ncol(dd3)*i, ]
  snp.dist <- dist(dd4)
  if (sum(is.na(snp.dist)) == 0) break
}
i # the final missing proportion
maxmissl = ifelse(maxmissl < i, maxmissl, i)
dd4 = dd3[rowSums(is.na(dd3)) < ncol(dd3)*maxmissl, ]
dim(dd4)

# remove monomorphic markers
tokeep = apply(dd4, 2, function(x){
  n1 = sum(x == 1, na.rm = T)
  n2 = sum(x == 2, na.rm = T)
  min(n1, n2) > 1 # now keep markers with at least 2 minor alleles
})

dd4 = dd4[,tokeep]

##
snp.dist <- dist(dd4)
fit <- hclust(snp.dist,method="ward.D2")
#str(fit)
pdf(file="dendrogram_of_lines.pdf", width = 14)
plot(fit, hang=-1, labels=FALSE)
# rect.hclust(fit, h = 10, border = "red")
dev.off()

# get big group number
maxheight = max(fit$height)
ng = cutree(fit, h = maxheight/4) # number of group, use 1/4 of the max height as cut threshold
table(ng)
dd4.1 = cbind(dd4, ng)
snp1.1 = cbind(snp, ng)

ordered_line_names = fit$labels[fit$order]
dd5 = dd4.1[ordered_line_names,]
snp2 = snp1.1[ordered_line_names, colnames(dd4.1)]
md2 = md[ordered_line_names,]

all(rownames(snp2) == rownames(md2))
outdata = cbind(md2, snp2) # with A, T, G, C alleles
outdata2 = cbind(md2, dd5) # with numbers: 0 is ref, 2 is alt, 1 is het.
outdata3 = cbind(md[rownames(snp),], snp) # non-filtered raw data
#outdata[1:4,1:30]

### write all output data
write.table(outdata, "haplotype_ordered_genotable.txt", row.names = F, sep="\t", quote=F, na = "-")
write.table(outdata2, "haplotype_ordered_genotable_numbers.txt", row.names = F, sep="\t", quote=F, na = "-")
write.table(outdata3, "raw_extracted_genotable.txt", row.names = F, sep="\t", quote=F, na = "-")
