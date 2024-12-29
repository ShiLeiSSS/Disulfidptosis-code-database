#if (!require("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install("maftools")


library(maftools)      
setwd("C:\\Users\\lexb\\Desktop\\disuloOmics\\12.maftools")      

geneRT=read.table("gene.txt", header=T, sep="\t", check.names=F, row.names=1)
gene=row.names(geneRT)

pdf(file="oncoplot.pdf", width=7, height=6.5)
maf=read.maf(maf="input.maf")
oncoplot(maf=maf, genes=gene, fontSize=0.8, draw_titv=T)
dev.off()

