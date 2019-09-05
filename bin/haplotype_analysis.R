#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("Please give the marker data and metadata file\n")
}

snpfile = args[1]
metadata = args[2]

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

##
snp.dist <- dist(dd3)
fit <- hclust(snp.dist,method="ward.D2")
str(fit)
pdf(file="dendrogram_of_lines.pdf", width = 10)
plot(fit, hang=-1, labels=FALSE)
# rect.hclust(fit, h = 10, border = "red")
dev.off()

snp2 = snp[fit$order,]
md2 = md[rownames(snp2),]
all(rownames(snp2) == rownames(md2))
outdata = cbind(md2, snp2)
#outdata[1:4,1:30]
write.table(outdata, "haplotype_ordered_genotable.txt", row.names = F, sep="\t", quote=F)

