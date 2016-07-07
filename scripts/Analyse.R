library(ggplot2)
library(reshape2)

# Set work directory.
setwd("C:/Users/Eric/Desktop/Amine")

#whichAligner <- "star"
whichAligner <- "bwa"

outputDir <- file.path("Output", whichAligner)
dir.create(outputDir, recursive=TRUE)

# Read in result files.
WT = read.table(file.path("Input", whichAligner, "results.WT"), sep="\t", quote="", header=TRUE)
Spt6 = read.table(file.path("Input", whichAligner, "results.Spt6-SA"), sep="\t", quote="", header=TRUE)
CK2 = read.table(file.path("Input", whichAligner, "results.CK2"), sep="\t", quote="", header=TRUE)

# Remove genes overlapping another gene in the opposite sense.
ambiguous.genes = scan("Input/ambiguous.genes", what=character())

# Remove multi-exon genes
multiexons = read.table("Input/Multiexon.bed", sep="\t", stringsAsFactors=FALSE)


# Subset to keep only exons with enough reads that the ratios won't fluctuate wildly.
# keep = (WT$Sense.All > 100) & (Spt6$Sense.All > 100) & (CK2$Sense.All > 100) & !(WT$Gene %in% ambiguous.genes) & !(WT$Gene %in% multiexons$V4)

#
min.sense.coverage = (WT$Sense.All > 100) & (Spt6$Sense.All > 100) & (CK2$Sense.All > 100)
single.exon.non.ambiguous = !(WT$Gene %in% ambiguous.genes) & !(WT$Gene %in% multiexons$V4)

# Minimume sense, anti-sense coverage at 5' and 3' ends
min.5.prime.coverage = WT$Sense.5 >= 20 & Spt6$Sense.5 >= 20 & CK2$Sense.5 >= 20
min.3.prime.coverage = WT$Sense.3 >= 10 & Spt6$Sense.3 >= 10 & CK2$Sense.3 >= 10

keep = min.sense.coverage & single.exon.non.ambiguous & min.5.prime.coverage & min.3.prime.coverage

#keep = !(WT$Gene %in% ambiguous.genes) & !(WT$Gene %in% multiexons$V4)
WT <- WT[keep,]
Spt6 <- Spt6[keep,]
CK2 <- CK2[keep,]

# Concatenate all results into a single matrix.
concatResults <- cbind(WT, Spt6[,6:11], CK2[,6:11])
colnames(concatResults) <- c(colnames(WT)[1:5], 
                             paste("WT", colnames(WT)[6:11], sep="."),
                             paste("Spt6", colnames(WT)[6:11], sep="."),
                             paste("CK2", colnames(WT)[6:11], sep="."))

# Add symbols to genes
library(org.Sc.sgd.db)
allGeneNames <- as.data.frame(org.Sc.sgdGENENAME)
concatResults$Symbol <- allGeneNames$gene_name[match(concatResults$Gene, allGeneNames$systematic_name)]                             
                             
# Calculate antisense percentage
for(lib in c("WT", "Spt6", "CK2")) {
    senseColumn <- paste(lib, ".Sense.All", sep="")
    antisenseColumn <- paste(lib, ".Antisense.All", sep="")
    ratioColumn <- paste(lib, ".Antisense.Percentage", sep="")
    concatResults[, ratioColumn] <- concatResults[,antisenseColumn] / (concatResults[,senseColumn] + concatResults[,antisenseColumn])
}

# Calculate antisense enrichment
concatResults$CK2.Antisense.enrichment <- log2(concatResults$"CK2.Antisense.Percentage" / concatResults$"WT.Antisense.Percentage")
concatResults$Spt6.Antisense.enrichment <- log2(concatResults$"Spt6.Antisense.Percentage" / concatResults$"WT.Antisense.Percentage")
                             
# Calculate sense 3'/5' ratio for all libraries.
concatResults$WT.3.on.5 <- (concatResults$WT.Sense.3 / concatResults$WT.Sense.5)
concatResults$Spt6.3.on.5 <- (concatResults$Spt6.Sense.3 / concatResults$Spt6.Sense.5)
concatResults$CK2.3.on.5 <- (concatResults$CK2.Sense.3 / concatResults$CK2.Sense.5)

# Comparison of 3' to 5' ratio between Spt6/CK2 and WT.
concatResults$Spt6.WT.3.on.5 <- log2(concatResults$Spt6.3.on.5 / concatResults$WT.3.on.5)
concatResults$CK2.WT.3.on.5 <- log2(concatResults$CK2.3.on.5 / concatResults$WT.3.on.5 )

# Rank various metrics.
concatResults$CK2.Antisense.rank <- order(concatResults$CK2.Antisense.enrichment, decreasing=TRUE)
concatResults$CK2.Antisense.ecdf <-  ecdf(concatResults$CK2.Antisense.enrichment)(concatResults$CK2.Antisense.enrichment)
concatResults$Spt6.Antisense.rank <- order(concatResults$Spt6.Antisense.enrichment, decreasing=TRUE)
concatResults$Spt6.Antisense.ecdf <-  ecdf(concatResults$Spt6.Antisense.enrichment)(concatResults$Spt6.Antisense.enrichment)

concatResults$CK2.3Enrichment.rank <- order(concatResults$CK2.WT.3.on.5, decreasing=TRUE)
concatResults$CK2.3Enrichment.ecdf <- ecdf(concatResults$CK2.3Enrichment.rank)(concatResults$CK2.3Enrichment.rank)
concatResults$Spt6.3Enrichment.rank <- order(concatResults$Spt6.WT.3.on.5, decreasing=TRUE)
concatResults$Spt6.3Enrichment.ecdf <- ecdf(concatResults$Spt6.3Enrichment.rank)(concatResults$Spt6.3Enrichment.rank)

# Write the table to disk.
write.table(concatResults, file.path(outputDir, "All Results.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

# Plot antisense percentage for all libraries.
antisenseDF <- data.frame(AntisensePerc = c(concatResults$WT.Antisense.Percentage, concatResults$Spt6.Antisense.Percentage, concatResults$CK2.Antisense.Percentage),
                          Library = rep(c("WT", "Spt6", "CK2"), each=nrow(concatResults)))
ggplot(data=antisenseDF) +
     geom_line(mapping=aes(x=AntisensePerc, color=Library), size=1, stat="density") +
     xlab("Percentage of antisense reads")
ggsave(file.path(outputDir, "Density plot of percentage of antisense reads.pdf"), width=7, height=7)

# Plot antisense enrichments for Spt6 and CK2
antisenseEnrichmentDF <- data.frame(Enrichment=c(concatResults$CK2.Antisense.enrichment, concatResults$Spt6.Antisense.enrichment),
                                    Contrast=rep(c("CK2 vs WT", "Spt6 vs WT")), each=nrow(concatResults))
ggplot(data=antisenseEnrichmentDF) +
     geom_line(mapping=aes(x=Enrichment, color=Contrast), size=1, stat="density") +
     xlab("Antisense enrichment relative to WT")
ggsave(file.path(outputDir, "Density plot of antisense enrichment.pdf"), width=7, height=7)                                    
                                    
# Get antisense confidence interval
sink(file=file.path(outputDir, "Statistical results.txt"))
cat("Is CK2 antisense enrichment different from 0?\n")
no.infinity <- concatResults$CK2.Antisense.enrichment
no.infinity[is.infinite(no.infinity)] <- NA
t.test(no.infinity)
cat("Is Spt6 antisense enrichment different from 0?\n")
no.infinity <- concatResults$Spt6.Antisense.enrichment
no.infinity[is.infinite(no.infinity)] <- NA
t.test(no.infinity)
sink(NULL)

# Plot distribution of 3' to 5' ratios.
threePrimeRatioDF <- data.frame(Ratio=log2(c(concatResults$WT.3.on.5, concatResults$Spt6.3.on.5, concatResults$CK2.3.on.5)),
                                Library = rep(c("WT", "Spt6", "CK2"), each=nrow(concatResults)))
ggplot(data=threePrimeRatioDF) +    
    geom_line(mapping=aes(x=Ratio, color=Library), size=1, stat="density") +
    xlab("log2(Ratio of the number of 3' reads over the number of 5' reads)")
ggsave(file.path(outputDir, "Density plot of 3 prime on 5 prime reads ratio.pdf"), width=7, height=7)    

# Plot distribution of 3' to 5' enrichments.
threePrimeEnrichmentDF <- data.frame(Enrichment=c(concatResults$CK2.WT.3.on.5, concatResults$Spt6.WT.3.on.5),
                                     Contrast=rep(c("CK2 vs WT", "Spt6 vs WT")), each=nrow(concatResults))
ggplot(data=antisenseEnrichmentDF) +
     geom_line(mapping=aes(x=Enrichment, color=Contrast), size=1, stat="density") +
     xlab("3' enrichment relative to WT")
ggsave(file.path(outputDir, "Density plot of three prime enrichment.pdf"), width=7, height=7)                                    

sink(file=file.path(outputDir, "Statistical results.txt"), append=TRUE)
cat("Is CK2 3' enrichment different from 0?\n")
t.test(concatResults$CK2.WT.3.on.5)
cat("Is Spt6 3' enrichment different from 0?\n")
t.test(concatResults$Spt6.WT.3.on.5)
sink(NULL)

gene.length.below.5K <- (concatResults$End - concatResults$Start) < 5000

# Gene length effect
ggplot(data=concatResults[gene.length.below.5K,], mapping=aes(y=CK2.WT.3.on.5, x=End-Start)) +
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    ylab("3' enrichment") + xlab("Gene length")
ggsave(file.path(outputDir, "Gene length effect on 3 prime enrichment (CK2).pdf"), width=7, height=7)
    
ggplot(data=concatResults[gene.length.below.5K,], mapping=aes(y=Spt6.WT.3.on.5, x=End-Start)) +
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    ylab("3' enrichment") + xlab("Gene length")
ggsave(file.path(outputDir, "Gene length effect on 3 prime enrichment (Spt6).pdf"), width=7, height=7)

# Read PolII enrichment data.
rnap2.enrichment = read.table("input/ck2.wt.rnap2.list.xls", sep="\t", header=TRUE)
rnap2.enrichment$Symbol = gsub(":.*$", "", rnap2.enrichment$X.1)
rnap2.enrichment$Gene = gsub("^.*:", "", rnap2.enrichment$X.1)

concatResults$RNAII.enrichment = rnap2.enrichment[match(concatResults$Gene, rnap2.enrichment$Gene), 1]

# Plot PolII enrichment vs antisense anrichment
ggplot(data=concatResults) + geom_point(mapping=aes(x=RNAII.enrichment, y=CK2.Antisense.enrichment))
ggsave("PolII vs Antisense enrichment.pdf")

ggplot(data=concatResults) + geom_point(mapping=aes(x=order(RNAII.enrichment), y=order(CK2.Antisense.enrichment)))
ggsave("PolII vs Antisense enrichment (Rank).pdf")

# Plot PolII enrichment vs 3' anrichment
ggplot(data=concatResults) + geom_point(mapping=aes(x=RNAII.enrichment, y=CK2.WT.3.on.5))
ggsave("PolII vs 3 prime enrichment.pdf")

ggplot(data=concatResults) + geom_point(mapping=aes(x=order(RNAII.enrichment), y=order(CK2.WT.3.on.5)))
ggsave("PolII vs 3 prime enrichment (Rank).pdf")

# Do metagene type plot.
# Build a data.frame describing the matrix files.
matrix.files = list.files("input/deeptools", pattern=".matrix$", full.names=TRUE)
input.files = data.frame(Filename  = matrix.files,
                         Name      = gsub("input/deeptools/", "", gsub("_sorted.matrix", "", matrix.files)),
                         stringsAsFactors=FALSE)

input.files$Strain    = factor(gsub("_.*", "", input.files$Name), levels=c("WT", "Spt6-SA", "CK2"))
input.files$Direction = factor(gsub(".*_", "", input.files$Name), levels=c("sense", "antisense"))

# Load all matrices and normalize relative to gene length.
raw.matrices = list()
matrices.rpk <- list()
for(i in 1:nrow(input.files)) {
    # Read the data from disk.
    raw.matrix = read.table(input.files$Filename[i], sep="\t", na.strings="nan", quote="", skip=1)
    raw.matrix[is.na(raw.matrix)] <- 0
    raw.matrices[[ input.files$Name[i] ]] <- raw.matrix
    
    # Drop annotation columns and normalize for gene lengths.
    gene.lengths = ((raw.matrix[,3] - raw.matrix[,2])) / 100
    matrices.rpk[[ input.files$Name[i] ]] <- apply(raw.matrix[,7:106] + 1, 2, "/", gene.lengths)
}

# Match coordinates matrix rows to the concatResults data.frame.
# Use coordinates, since this is the only information present in both.
concat.coords = paste(concatResults$Chr, concatResults$Start, concatResults$End)
matrix.coords = paste(raw.matrix[,1], raw.matrix[,2], raw.matrix[,3])

# Only keep rows from the same subset as the one used for the previous analyses.
gene.subset.all <- match(concat.coords, matrix.coords)

# Produces graph for all genes, but also for specific genes.
genesOfInterest <- c("All", "YHL035C", "YIL156W", "YHR186C", "YML123C", "YGL092W","YDR227W")
for(gene in genesOfInterest) {
    
    # Subset genes if we want to graph a specific one.
    if(gene != "All") {
        gene.subset <- gene.subset.all[which(as.character(concatResults$Gene) == gene)]
    } else {
        gene.subset <- gene.subset.all
    }

    # Produce a final matrix giving the read counts for all 100 bins for all read categories.
    final.matrix <- matrix(0, nrow=length(matrix.files), ncol=100)
    
    # Loop over mutants calculating coverage across genes.
    for(strain in unique(input.files$Strain)) {
        this.strain = input.files$Strain == strain
        sense.direction = input.files$Direction=="sense"
        
        sense.index = which(this.strain & sense.direction)
        antisense.index = which(this.strain & !sense.direction)
        
        # Calculate the total number of reads for this strain for normalization.
        reads.sense     = sum(matrices.rpk[[sense.index]][gene.subset,], na.rm=TRUE)
        reads.antisense = sum(matrices.rpk[[antisense.index]][gene.subset,], na.rm=TRUE)
        totalReads = reads.sense + reads.antisense
        
        final.matrix[sense.index,]     <- apply(matrices.rpk[[ sense.index     ]][gene.subset,, drop=FALSE], 2, mean, na.rm=TRUE) / totalReads * 1000000
        final.matrix[antisense.index,] <- apply(matrices.rpk[[ antisense.index ]][gene.subset,, drop=FALSE], 2, mean, na.rm=TRUE) / totalReads * 1000000
    }

    # Give the final matrix dimensions' names.
    rownames(final.matrix) <- names(matrices.rpk)
    colnames(final.matrix) <- 1:100    

    # Turn the matrix into a data frame to plot it with ggplot.
    finalDF <- melt(log2(final.matrix))
    colnames(finalDF) <- c("Library", "TSS.distance", "RPM")
    
    # Define strain/direction as factors.
    finalDF$Strain = factor(gsub("_.*", "", finalDF$Library), levels=c("WT", "Spt6-SA", "CK2"))
    finalDF$Direction = factor(gsub(".*_", "", finalDF$Library), levels=c("sense", "antisense"))
    
    # Generate the plot.
    ggplot(data=finalDF) +
        geom_line(mapping=aes(x=TSS.distance, y=RPM, colour=Strain, linetype=Direction), size=2) +
        labs(y="log2(RPKM)")
        
    # Save it!
    ggsave(paste0("Metagene for ", gene, ".pdf"))
}


# Compute correlation between antisense enrichments.
no.infinity = subset(concatResults, !is.infinite(CK2.Antisense.enrichment) & !is.nan(CK2.Antisense.enrichment) & !is.infinite(Spt6.Antisense.enrichment) & !is.nan(Spt6.Antisense.enrichment))

cor.results = with(no.infinity, cor.test(CK2.Antisense.enrichment, Spt6.Antisense.enrichment))
cor.interval = sprintf("R^2 95%% conf. interval: %0.2f-%0.2f", cor.results$conf.int[1], cor.results$conf.int[2])

ggplot(no.infinity, aes(x=CK2.Antisense.enrichment, y=Spt6.Antisense.enrichment)) +
    geom_point(mapping=) +
    geom_smooth(method=lm, se=FALSE) +
    annotate("text", x=max(no.infinity$CK2.Antisense.enrichment), y=max(no.infinity$Spt6.Antisense.enrichment), label=cor.interval, hjust=1, size=5)

ggsave("Antisense enrichment correlation.pdf")    

# Do the same for 3' enrichments.
no.infinity = subset(concatResults, !is.infinite(CK2.WT.3.on.5) & !is.nan(CK2.WT.3.on.5) & !is.infinite(Spt6.WT.3.on.5) & !is.nan(Spt6.WT.3.on.5))

# Remove outliers
no.infinity = subset(no.infinity, CK2.WT.3.on.5 < quantile(CK2.WT.3.on.5, 0.995) &
                                  CK2.WT.3.on.5 > quantile(CK2.WT.3.on.5, 0.005) & 
                                  Spt6.WT.3.on.5 < quantile(Spt6.WT.3.on.5, 0.995) & 
                                  Spt6.WT.3.on.5 > quantile(Spt6.WT.3.on.5, 0.005))

cor.results = with(no.infinity, cor.test(CK2.WT.3.on.5, Spt6.WT.3.on.5))
cor.interval = sprintf("R^2 95%% conf. interval: %0.2f-%0.2f", cor.results$conf.int[1], cor.results$conf.int[2])

ggplot(no.infinity, aes(x=CK2.WT.3.on.5, y=Spt6.WT.3.on.5)) +
    geom_point(mapping=) +
    geom_smooth(method=lm, se=FALSE) +
    annotate("text", x=max(no.infinity$CK2.WT.3.on.5), y=max(no.infinity$Spt6.WT.3.on.5), label=cor.interval, hjust=1, size=5)

ggsave("Three prime enrichment correlation.pdf")    

