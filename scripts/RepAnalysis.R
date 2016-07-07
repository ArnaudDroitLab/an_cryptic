library(ggplot2)
library(reshape2)
library(eric.utils)
library(org.Sc.sgd.db)
library(plyr)

# Utility function to work with data-frames: subset columns
# and only keep those that fit a regular expression.
grep.col <- function(df, regex) {
    return(df[,grepl(regex, colnames(df))])
}

# Get output directory, which in our case is also the input directory.
output.dir <- file.path("output-new")


# Read in result files.
count.path = file.path(output.dir, "5-3Counts")
count.files = list.files(count.path, "*.results")
all.counts = read.identical(file.path(count.path, count.files), 1:5, 6:11, gsub(".results", "", count.files), quote="")

# Define useful column subsets
sense.all = grep.col(all.counts, "Sense.All")
antisense.all = grep.col(all.counts, "Antisense.All")
sense.5 = grep.col(all.counts, "Sense.5")
sense.3 = grep.col(all.counts, "Sense.3")

# Calculate new columns



# Perform gene filtering:
# Subset on minimum number of reads so we don't deal with widly fluctuating ratios.
sense.all = grep.col(all.counts, "Sense.All")
min.sense.overall = apply(sense.all, 1, mean) > 100 & apply(sense.all > 50, 1, all)
antisense.all = grep.col(all.counts, "Antisense.All")
min.antisense.overall = apply(antisense.all > 1, 1, all)
sense.5 = grep.col(all.counts, "Sense.5")
min.sense.5 = apply(sense.5 >= 20, 1, all)
sense.3 = grep.col(all.counts, "Sense.3")
min.sense.3 = apply(sense.3 >= 10, 1, all)

# Remove genes overlapping another gene in the opposite sense.
ambiguous.genes = scan("Input/ambiguous.genes", what=character())
non.ambiguous = !(all.counts$Gene %in% ambiguous.genes)

# Remove multi-exon genes
multiexons = read.table("Input/Multiexon.bed", sep="\t", stringsAsFactors=FALSE)
single.exon = !(all.counts$Gene %in% multiexons$V4)

# Minimume sense, anti-sense coverage at 5' and 3' ends
keep = min.sense.overall & min.antisense.overall & min.sense.5 & min.sense.3 & single.exon & non.ambiguous
all.counts = all.counts[keep, ]

# Add symbols to genes
allGeneNames <- as.data.frame(org.Sc.sgdGENENAME)
all.counts$Symbol <- allGeneNames$gene_name[match(all.counts$Gene, all.counts$systematic_name)]                             

                             
# Calculate antisense percentage
antisense.percentage = grep.col(all.counts, "Antisense.All")  / (grep.col(all.counts, "Sense.All") + grep.col(all.counts, "Antisense.All"))
colnames(antisense.percentage) = gsub("Antisense.All", "Antisense.Percentage", colnames(antisense.percentage))

# Calculate antisense percentage
antisense.ratio = log2(grep.col(all.counts, "Antisense.All")  / grep.col(all.counts, "Sense.All"))
colnames(antisense.ratio) = gsub("Antisense.All", "Antisense.Ratio", colnames(antisense.ratio))


# Calculate antisense enrichment
antisense.enrichment = cbind(log2(grep.col(antisense.percentage, "CK2") / grep.col(antisense.percentage, "WT")),
                             log2(grep.col(antisense.percentage, "Spt6") / grep.col(antisense.percentage, "WT")))
colnames(antisense.enrichment) = gsub("Antisense.Percentage", "Antisense.enrichment", colnames(antisense.enrichment))

# Calculate sense 3'/5' ratio for all libraries.
three.on.five = grep.col(all.counts, "Sense.3") / grep.col(all.counts, "Sense.5")
colnames(three.on.five) = gsub("Sense.3", "Three.On.Five", colnames(three.on.five))

three.on.five.ratio = log2(grep.col(all.counts, "Sense.3") / grep.col(all.counts, "Sense.5"))
colnames(three.on.five) = gsub("Sense.3", "Three.On.Five", colnames(three.on.five))


three.on.five.enrichment = cbind(log2(grep.col(three.on.five, "CK2") / grep.col(three.on.five, "WT")),
                                 log2(grep.col(three.on.five, "Spt6") / grep.col(three.on.five, "WT")))
colnames(three.on.five.enrichment) = gsub("Three.On.Five", "Three.On.Five.enrichment", colnames(three.on.five.enrichment))        
        
        
row.t.test <- function(x) {
    CK2 = t.test(x[grepl("CK2", names(x))], x[grepl("WT", names(x))])$p.value
    Spt6 = t.test(x[grepl("Spt6", names(x))], x[grepl("WT", names(x))])$p.value
    return(c(CK2=CK2, Spt6=Spt6))
}        
        
# T-tests
antisense.enrichment.pval = t(apply(antisense.percentage, 1, row.t.test))
antisense.ratio.pval = t(apply(antisense.ratio, 1, row.t.test))
three.on.five.enrichment.pval = t(apply(three.on.five, 1, row.t.test))
three.on.five.ratio.pval = t(apply(three.on.five.ratio, 1, row.t.test))        


test = cbind(all.counts, 
             antisense.percentage, 
             antisense.ratio, 
             antisense.enrichment, 
             three.on.five.ratio, 
             three.on.five.enrichment, 
             antisense.enrichment.pval, 
             antisense.ratio.pval, 
             three.on.five.enrichment.pval, three.on.five.ratio.pval)


antisense.ratio.df = melt(antisense.ratio)
antisense.ratio.df$Mutant = gsub("(-2)*.Antisense.Ratio", "", antisense.ratio.df$variable)
antisense.ratio.df$Rep = ifelse(grepl("-2", antisense.ratio.df$variable), 2, 1)
means.df = ddply(antisense.ratio.df, Mutant~Rep, summarize, mean=mean(value))

ggplot(antisense.ratio.df) + 
    geom_density(mapping=aes(x=value, fill=Mutant), alpha=0.6) + 
    geom_vline(data=means.df, mapping=aes(xintercept=mean), linetype=2) + 
    xlab("log2(Antisense reads / Sense reads)") + 
    facet_grid(Mutant~Rep)
    
ggsave("output/Antisense ratios.pdf")    

# Replicate correlation
rep.correlation = rbind(setNames(grep.col(antisense.ratio, "CK2"), c("Rep1","Rep2")),
                        setNames(grep.col(antisense.ratio, "Spt6"), c("Rep1","Rep2")),
                        setNames(grep.col(antisense.ratio, "WT"), c("Rep1","Rep2")))
rep.correlation$Mutant = rep(c("CK2", "Spt6", "WT"), each=nrow(antisense.ratio))
rep.correlation.melt = melt(rep.correlation, "Mutant", variable.name = "Replicate")

pairwise.correlations <- function(input.df, colname, label) {
    # Plot antisense ratios correlations
    full.df <- NULL
    cor.df <- NULL
    max.limit = max(input.df)
    for(i in c("CK2", "Spt6-SA", "WT")) {
        for(j in c("CK2", "Spt6-SA", "WT")) {
            x <- input.df[,paste0(i, colname)]
            y <- input.df[,paste0(j, "-2", colname)]
            new.df = data.frame(Rep1=x,
                                Rep2=y,
                                Mutant1=i,
                                Mutant2=j,
                                stringsAsFactors=FALSE)
            new.cor = data.frame(Cor=cor(x, y),
                                 x=max.limit,
                                 y=max.limit,
                                 Mutant1=i,
                                 Mutant2=j)
                                
            if(is.null(full.df)) {
                full.df = new.df
                cor.df = new.cor
            } else {
                full.df = rbind(full.df, new.df)
                cor.df = rbind(cor.df, new.cor)
            }
        }
    }
    
    ggplot(full.df, aes(x=Rep1, y=Rep2)) + 
        geom_point() + 
        geom_smooth(method=lm, se=FALSE) + 
        geom_label(data=cor.df, mapping=aes(label=sprintf("R=%.2f", Cor), x=x, y=y), hjust=1, vjust=1) +    
        facet_grid(Mutant1~Mutant2)
    
    ggsave(paste0("output/Pairwise correlations of ", label, ".pdf"))    
}

pairwise.correlations(antisense.ratio, ".Antisense.Ratio", "antisense ratios")
pairwise.correlations(three.on.five.ratio, ".Sense.3.", "three prime to Five prime ratios")

wt.enrichment <- function(input.df, colname, label, ymax) {
    # Plot antisense enrichments distributions
    full.df <- NULL
    cont.df <- NULL
    #max.limit = max(antisense.ratio)
    for(i in c("CK2", "Spt6-SA")) {
        for(j in 1:2) {
            rep.str = ifelse(j==2, "-2", "")
            mutant <- antisense.ratio[,paste0(i, rep.str, ".Antisense.Ratio")]
            wt <- antisense.ratio[,paste0("WT", rep.str, ".Antisense.Ratio")]
            enrichment = mutant - wt
            
            stat.test = t.test(enrichment)
            
            new.df = data.frame(Enrichment=enrichment,
                                Rep=j,
                                Mutant=i,
                                stringsAsFactors=FALSE)
                                
            new.conf = data.frame(Mean=mean(enrichment),
                                  low.conf=stat.test$conf.int[1],
                                  high.conf=stat.test$conf.int[2],
                                  p="p<0.001",
                                  x=max(enrichment),
                                  y=0.4,
                                  Rep=j,
                                  Mutant=i,
                                  stringsAsFactors=FALSE)
                               
            if(is.null(full.df)) {
                full.df = new.df
                conf.df = new.conf
            } else {
                full.df = rbind(full.df, new.df)
                conf.df = rbind(conf.df, new.conf)
            }
        }
    }
    
    ggplot(full.df, aes(x=Enrichment, fill=Mutant)) + 
        geom_density(alpha=0.6) + 
        geom_vline(data=conf.df, mapping=aes(xintercept=Mean), linetype=2) +
        geom_label(data=conf.df, mapping=aes(x=x, y=y, label=p), hjust=1, vjust=1) +
        facet_grid(Rep~Mutant)
    
    ggsave("output/Enrichment of .pdf")   
    
}

wt.enrichment(antisense.ratio, ".Antisense.Ratio", "antisense ratios", 0.4)
wt.enrichment(three.on.five.ratio, ".Sense.3.", "three prime to Five prime ratios", 1.4)


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
write.table(concatResults, file.path(output.dir, "All Results.txt"), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
                                    
# Get antisense confidence interval
sink(file=file.path(output.dir, "Statistical results.txt"))
cat("Is CK2 antisense enrichment different from 0?\n")
no.infinity <- concatResults$CK2.Antisense.enrichment
no.infinity[is.infinite(no.infinity)] <- NA
t.test(no.infinity)
cat("Is Spt6 antisense enrichment different from 0?\n")
no.infinity <- concatResults$Spt6.Antisense.enrichment
no.infinity[is.infinite(no.infinity)] <- NA
t.test(no.infinity)
sink(NULL)

sink(file=file.path(output.dir, "Statistical results.txt"), append=TRUE)
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
ggsave(file.path(output.dir, "Gene length effect on 3 prime enrichment (CK2).pdf"), width=7, height=7)
    
ggplot(data=concatResults[gene.length.below.5K,], mapping=aes(y=Spt6.WT.3.on.5, x=End-Start)) +
    geom_point() +
    geom_smooth(method=lm, se=FALSE) +
    ylab("3' enrichment") + xlab("Gene length")
ggsave(file.path(output.dir, "Gene length effect on 3 prime enrichment (Spt6).pdf"), width=7, height=7)

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

matrix.path = file.path(output.dir, "Metagene")
matrix.files = list.files(matrix.path, pattern=".matrix$")
input.files = data.frame(Filename  = file.path(matrix.path, matrix.files),
                         Name      = gsub(".sorted.matrix", "", matrix.files),
                         stringsAsFactors=FALSE)

input.files$Strain    = factor(gsub("(-2)*\\..*", "", input.files$Name), levels=c("WT", "Spt6-SA", "CK2"))
input.files$Direction = factor(gsub(".*\\.", "", input.files$Name), levels=c("sense", "antisense"))
input.files$Replicate = factor(ifelse(grepl("-2", input.files$Name), 2, 1), levels=c("1", "2"))

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
concat.coords = paste(all.counts$Chr, all.counts$Start, all.counts$End)
matrix.coords = paste(raw.matrix[,1], raw.matrix[,2], raw.matrix[,3])

# Only keep rows from the same subset as the one used for the previous analyses.
gene.subset.all <- match(concat.coords, matrix.coords)

# Produces graph for all genes, but also for specific genes.
genesOfInterest <- c("All", "YHL035C", "YIL156W", "YHR186C", "YML123C", "YGL092W","YDR227W")
for(gene in genesOfInterest) {
    
    # Subset genes if we want to graph a specific one.
    if(gene != "All") {
        gene.subset <- gene.subset.all[which(as.character(all.counts$Gene) == gene)]
    } else {
        gene.subset <- gene.subset.all
    }

    # Produce a final matrix giving the read counts for all 100 bins for all read categories.
    final.matrix <- matrix(0, nrow=length(matrix.files), ncol=100)
    
    # Loop over mutants calculating coverage across genes.
    for(strain in unique(input.files$Strain)) {
        for(replicate.id in 1:2) {
            this.strain = input.files$Strain == strain & input.files$Replicate==replicate.id
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
    }

    # Give the final matrix dimensions' names.
    rownames(final.matrix) <- names(matrices.rpk)
    colnames(final.matrix) <- 1:100    

    # Turn the matrix into a data frame to plot it with ggplot.
    finalDF <- melt(log2(final.matrix))
    colnames(finalDF) <- c("Library", "TSS.distance", "RPM")
    
    # Define strain/direction as factors.
    #finalDF$Strain = factor(gsub("_.*", "", finalDF$Library), levels=c("WT", "Spt6-SA", "CK2"))
    #finalDF$Direction = factor(gsub(".*_", "", finalDF$Library), levels=c("sense", "antisense"))
    #finalDF$Replicate = factor(gsub(".*_", "", finalDF$Library), levels=c("sense", "antisense"))
    
    finalDF$Strain = input.files$Strain
    finalDF$Direction = input.files$Direction
    finalDF$Replicate = input.files$Replicate
        
    # Generate the plot.
    ggplot(data=finalDF) +
        geom_line(mapping=aes(x=TSS.distance, y=RPM, colour=Strain, linetype=Direction), size=2) +
        labs(y="log2(RPKM)") +
        facet_grid(~Replicate)
        
    # Save it!
    ggsave(file.path(output.dir, paste0("Metagene for ", gene, ".pdf")))
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

