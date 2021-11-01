#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args) < 3) {
  stop("Please give the marker data, maximum missing data for lines, and maximum missing data for markers\n")
}

snpfile = args[1]
maxmissl = as.numeric(args[2]) # maximum missing values to keep a line
maxmissm = as.numeric(args[3]) # maximum missing values to keep a marker

dd= read.delim(snpfile, check.names = F, as.is = T, na.strings = "N")
dim(dd)
dd = as.matrix(dd)
nc = ncol(dd)
nr = nrow(dd)

snp = t(dd[,6:nc])
#snp[1:5,1:5]
snpname = apply( dd[ , c(1,2,4,5) ] , 1 , paste , collapse = "-" )
#snpname
dd2 = sapply(1:nr, function(x) {
  #tt = factor(dd[x, 6:nc], levels = c(dd[x,4], "H", dd[x,5]))
  tt = factor(dd[x, 6:nc])
  tt2 = as.numeric(tt)
})


dimnames(dd2) = list(colnames(dd)[6:nc], snpname)
dimnames(snp) = list(colnames(dd)[6:nc], snpname)

dd2[1:4,1:4]
snp[1:4,1:4]

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
  n1 = sum(x == 0, na.rm = T)
  n2 = sum(x == 1, na.rm = T)
  min(n1, n2) > 0 # now keep markers with at least 1 minor alleles
})
sum(tokeep)

dd4 = dd4[,tokeep]
dim(dd4)


##
snp.dist <- dist(dd4)
fit <- hclust(snp.dist,method="ward.D2")
ht = max(fit$height) # max height
#str(fit)
pdf(file="dendrogram_of_lines.pdf", width = 10+round(nrow(dd4)/35) )
plot(fit, hang=-1, cex=1/(1+nrow(dd4)/250))
rect.hclust(fit, h = ht/5, border = "blue")
dev.off()

# get big group number
#maxheight = max(fit$height)
ng = cutree(fit, h = ht/5) # number of group, use 5 as cut threshold
table(ng)

ordered_line_names = fit$labels[fit$order]
ng = ng[ordered_line_names]
dd5 = dd4[ordered_line_names,]
snp2 = snp[ordered_line_names, colnames(dd4)]
snp3 = t(snp2)
cn = colnames(snp2)
snp.infor = data.frame(do.call('rbind', strsplit(cn,'-',fixed=TRUE)))
colnames(snp.infor) = c("Chrom", "Pos", "Ref", "Alt")
snp4 = cbind(snp.infor, snp3)
outdata = rbind(as.matrix(snp4), c("Group", rep("-", 3), ng)) # with A, T, G, C alleles


### write all output data
write.table(outdata, "haplotype_ordered_genotable.txt", row.names = F, sep="\t", quote=F, na = "-")

