#setwd("C:/Dev/Projects/an_cryptic")
library(DESeq)

# List and load all count files.
raw_counts_files = list(CK2.rep.1  = "output-new/raw_counts/CK2.readcounts.csv",
                        CK2.rep.2  = "output-new/raw_counts/CK2-2.readcounts.csv",
                        Spt6.rep.1 = "output-new/raw_counts/Spt6-SA.readcounts.csv",
                        Spt6.rep.2 = "output-new/raw_counts/Spt6-SA-2.readcounts.csv",
                        WT.rep.1   = "output-new/raw_counts/WT.readcounts.csv",
                        WT.rep.2   = "output-new/raw_counts/WT-2.readcounts.csv")
raw_counts_tables = lapply(raw_counts_files, read.table, sep="\t")

# Counts are on the second column.
raw_counts = as.data.frame(lapply(raw_counts_tables, "[[", 2))

# Transcript ID is on the first column.
rownames(raw_counts) <- raw_counts_tables[[1]][,1]

# Make sure all row names are identical
all.row.names = as.data.frame(lapply(raw_counts_tables, "[[", 1))
stopifnot(!any(!apply(all.row.names, 1, function(x){length(unique(x))==1})))

# Remvoe the meta-rows.
raw_counts = raw_counts[!grepl("^__", rownames(raw_counts)),]

# Get into a CountDataSet object and normalize
conditions = gsub("\\..*", "", colnames(raw_counts))
cds <- newCountDataSet(raw_counts, conditions)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds, method="pooled")

ck2.deseq <- nbinomTest(cds, "WT", "CK2")
spt6.deseq <- nbinomTest(cds, "WT", "Spt6")

# Build DGE list
library(edgeR)

libSize <- colSums(raw_counts)
dge.list=DGEList(counts=as.matrix(raw_counts), group=conditions , lib.size=libSize)

# Normalize counts
dTMM <- calcNormFactors(dge.list)
dTMM <- estimateCommonDisp(dTMM)
dTMM <- estimateTagwiseDisp(dTMM)
ck2.de <- exactTest(dTMM, c("WT", "CK2"))
spt6.de  <- exactTest(dTMM, c("WT", "Spt6"))

both.dge = list(CK2=list(DESeq=ck2.deseq, edgeR=ck2.de),
                Spt6=list(DESeq=spt6.deseq, edgeR=spt6.de))

to.df <- function(dge.results, label, p.threshold=0.05, fc.threshold=log2(1.5)) {
    res = data.frame(edgeR.logFC=dge.results$edgeR$table$logFC,
                     DESeq.logFC=dge.results$DESeq$log2FoldChange,
                     edgeR.pval=dge.results$edgeR$table$PValue,
                     DESeq.pval=dge.results$DESeq$pval)
    res$edgeR.adjpval = p.adjust(res$edgeR.pval, method="fdr", n=nrow(res)*2)
    res$DESeq.adjpval = p.adjust(res$DESeq.pval, method="fdr", n=nrow(res)*2)
    
    res$significant = res$edgeR.adjpval <= 0.05 & res$DESeq.adjpval <= 0.05 &
                      abs(res$edgeR.logFC) >= fc.threshold & abs(res$DESeq.logFC) >= fc.threshold
    
    res$direction = ifelse(!res$significant, "None", ifelse(res$edgeR.logFC > 0, "Up", "Down"))
    colnames(res) = paste0(label, ".", colnames(res))
    return(res)
}
     
     
combined = cbind(Gene=ck2.deseq$id, to.df(both.dge$CK2, "CK2"), to.df(both.dge$Spt6, "Spt6"))

combined$CK2.direction = factor(combined$CK2.direction)
combined$Spt6.direction = factor(combined$Spt6.direction)

write.table(combined, "output-new/analysis/DGE/DGE results.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)



# Correlation across all samples.
normalized.counts <- t( t(counts(cds)) / sizeFactors(cds) )

conditions = colnames(normalized.counts)
expr.df = data.frame()
cor.df = data.frame() 
for(i in 1:(length(conditions))) {
  cond1 = conditions[i]
  for(j in (i):length(conditions)) {
    cond2 = conditions[j]
    count1 = normalized.counts[,cond1]
    count2 = normalized.counts[,cond2]
    cor.val = cor(count1, count2)
    cor.label = sprintf("R=%.2f", cor.val)
    
    expr.df = rbind(expr.df, data.frame(Cond1=cond1, Cond2=cond2, Count1=count1, Count2=count2))
    cor.df = rbind(cor.df, data.frame(Cond1=cond1, Cond2=cond2, Cor=cor(count1, count2), CorLabel=cor.label))
  }
}
cor.df$Count1 = min(expr.df$Count1)
cor.df$Count2=max(expr.df$Count2)

library(ggplot2)
ggplot(expr.df, aes(x=log2(Count1+1), y=log2(Count2+1))) +
    geom_point() +
    geom_label(data=cor.df, mapping=aes(label=CorLabel), hjust=-0.1, vjust=1) +
    facet_grid(Cond1~Cond2)
ggsave("output-new/analysis/DGE/Correlation.pdf", width=14, height=14)

# Compare with known transcription levels
young.levels = read.table("Input/orf_transcriptome.txt", header=TRUE, sep="\t", comment="", stringsAsFactors=FALSE)
for(column.name in colnames(young.levels)[-1]) {
    young.levels[[column.name]][grepl("#N/A", young.levels[[column.name]])] = NA
    young.levels[[column.name]] <- as.numeric(young.levels[[column.name]])
}

young.match = match(rownames(normalized.counts), young.levels$ORF)
comp.matrix = cbind(normalized.counts, young.levels[young.match,-1])
comp.matrix = comp.matrix[!apply(is.na(comp.matrix), 1, any),]
comp.rank = apply(comp.matrix, 2, rank)
res = data.frame()
for(i in c(1:6, 7:12, 19:24)) {
    res = rbind(res, data.frame(Sample=colnames(comp.rank)[i],
                                ExpressionLevel=comp.rank[,"ExpressionLevel"],
                                RNALevel=comp.rank[,i]))
}

ggplot(res, aes(x=ExpressionLevel, y=RNALevel)) + geom_point(alpha=0.1) + facet_wrap(~Sample, nrow=3)


comp.quantile = floor(comp.rank / nrow(comp.rank) * 5)
apply(apply(comp.quantile[,1:(ncol(comp.quantile)-1)], 2, "==", comp.quantile[,ncol(comp.quantile)]), 2, sum)/nrow(comp.quantile)
apply(apply(comp.quantile[,1:(ncol(comp.quantile)-1)], 2, "==", comp.quantile[,7]), 2, sum)/nrow(comp.quantile)

hailu.expr = read.table("input/GSE85460_edgeR_Original_Unprocessed_Expression_data.txt", header=TRUE)
hailu.match = match(rownames(normalized.counts), hailu.expr$Gene)
comp.matrix = cbind(normalized.counts, hailu.expr[hailu.match,-1], young.levels[young.match,2])
comp.matrix = comp.matrix[!apply(is.na(comp.matrix), 1, any),]
comp.rank = apply(comp.matrix, 2, rank)

pang.expr = read.table("Input/GSE80357_counts_of_mapped_reads.csv", sep=",", header=TRUE)
pang.match = match(rownames(normalized.counts), pang.expr$oln_id)
comp.matrix = cbind(normalized.counts, hailu.expr[hailu.match,-1], pang.expr[pang.match, -1], young.levels[young.match,2, drop=FALSE])
comp.matrix = comp.matrix[!apply(is.na(comp.matrix), 1, any),]
comp.rank = apply(comp.matrix, 2, rank)
comp.quantile = floor(comp.rank / nrow(comp.rank) * 5)
apply(apply(comp.quantile[,1:(ncol(comp.quantile)-1)], 2, "==", comp.quantile[,ncol(comp.quantile)]), 2, sum)/nrow(comp.quantile)
