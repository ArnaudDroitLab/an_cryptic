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
base.output.dir <- file.path("output-new")
output.dir = file.path(base.output.dir, "analysis")
dir.create(output.dir, recursive=TRUE, showWarnings=FALSE)


# Read in result files.
count.path = file.path(base.output.dir, "5-3Counts")
count.files = list.files(count.path, "*.results")
all.counts = read.identical(file.path(count.path, count.files), 1:5, 6:11, gsub(".results", "", count.files), quote="")

# Define useful column subsets
sense.all = grep.col(all.counts, "Sense.All")
antisense.all = grep.col(all.counts, "Antisense.All")
sense.5 = grep.col(all.counts, "Sense.5")
sense.3 = grep.col(all.counts, "Sense.3")

# Calculate new columns
calculate.enrichment <- function(input.df, old.suffix, new.suffix) {
    result = cbind(grep.col(input.df, "CK2") - grep.col(input.df, "WT"),
                   grep.col(input.df, "Spt6") - grep.col(input.df, "WT"))
    colnames(result) = gsub(old.suffix, new.suffix, colnames(result))

    return(result)
}

calculate.ratio <- function(col1, col2, new.suffix) {
    result = log2(grep.col(all.counts, col1) / grep.col(all.counts, col2))
    colnames(result) = gsub(col1, new.suffix, colnames(result))
    return(result)    
}

# Calculate antisense ratio
antisense.ratio = calculate.ratio("Antisense.All", "Sense.All", "Antisense.Ratio")
antisense.enrichment = calculate.enrichment(antisense.ratio, "Antisense.Ratio", "Antisense.Enrichment")

# Calculate sense 3'/5' ratio for all libraries.
three.on.five.ratio = calculate.ratio("Sense.3.", "Sense.5.", "Three.On.Five.Ratio")
three.on.five.enrichment = calculate.enrichment(three.on.five.ratio, "Three.On.Five.Ratio", "Three.On.Five.Enrichment")

#all.counts = cbind(all.counts, antisense.ratio, antisense.enrichment, three.on.five.ratio, three.on.five.enrichment)   
all.counts = cbind(all.counts, 
                   antisense.ratio,
                   antisense.enrichment,
                   three.on.five.ratio,
                   three.on.five.enrichment)


# Perform gene filtering:
# Subset on minimum number of reads so we don't deal with widly fluctuating ratios.
min.sense.overall = apply(sense.all, 1, mean) > 100 & apply(sense.all > 50, 1, all)
min.antisense.overall = apply(antisense.all > 1, 1, all)
min.sense.5 = apply(sense.5 >= 20, 1, all)
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
all.counts$Symbol <- mapIds(org.Sc.sgd.db, all.counts$Gene, "GENENAME", "ENSEMBL", "first")
        
row.t.test <- function(x) {
    CK2 = t.test(x[grepl("CK2", names(x))], x[grepl("WT", names(x))])$p.value
    Spt6 = t.test(x[grepl("Spt6", names(x))], x[grepl("WT", names(x))])$p.value
    return(c(CK2=CK2, Spt6=Spt6))
}        

all.t.test <- function(input.df, colname) {
    results = t(apply(input.df, 1, row.t.test))
    colnames(results) <- paste(colnames(results), colname, "Pval", sep=".")
    return(results)
}
        
# T-tests
# antisense.enrichment.pval = t(apply(antisense.percentage, 1, row.t.test))
# antisense.ratio.pval = t(apply(antisense.ratio, 1, row.t.test))
# three.on.five.enrichment.pval = t(apply(three.on.five, 1, row.t.test))
# three.on.five.ratio.pval = t(apply(three.on.five.ratio, 1, row.t.test))        

all.counts = cbind(all.counts,
                   all.t.test(antisense.ratio[keep, ], "Antisense.Ratio"),
                   all.t.test(three.on.five.ratio[keep, ], "Three.On.Five.Ratio"))

write.table(all.counts, file=file.path(output.dir, "Genewise data.txt"), sep="\t", col.names=TRUE, row.names=FALSE)                   
                   
plot.raw.ratios <- function(input.df, colname, x.label, file.label) {
    melt.df = melt(input.df)
    melt.df$Mutant = gsub(paste0("(-2)*.", colname), "", melt.df$variable)
    melt.df$Rep = ifelse(grepl("-2", melt.df$variable), 2, 1)
    means.df = ddply(melt.df, Mutant~Rep, summarize, mean=mean(value, na.rm=TRUE))
    
    ggplot(melt.df) + 
        geom_density(mapping=aes(x=value, fill=Mutant), alpha=0.6) + 
        geom_vline(data=means.df, mapping=aes(xintercept=mean), linetype=2) + 
        xlab(x.label) + 
        facet_grid(Mutant~Rep)
        
    ggsave(file.path(output.dir, paste0(file.label, ".pdf")))
}

plot.raw.ratios(antisense.ratio[keep,], "Antisense.Ratio", "log2(Antisense reads / Sense reads)", "Antisense ratios")
plot.raw.ratios(three.on.five.ratio[keep,], "Three.On.Five.Ratio", "log2(3' reads / 5' reads)", "Three prime to five prime ratios")

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
    
    ggsave(file.path(output.dir, paste0("Pairwise correlations of ", label, ".pdf")))    
}

pairwise.correlations(antisense.ratio[keep,], ".Antisense.Ratio", "antisense ratios")
pairwise.correlations(three.on.five.ratio[keep,], ".Three.On.Five.Ratio", "three prime to five prime ratios")

wt.enrichment <- function(input.df, colname, label, ymax) {
    # Plot antisense enrichments distributions
    full.df <- NULL
    cont.df <- NULL
    #max.limit = max(antisense.ratio)
    for(i in c("CK2", "Spt6-SA")) {
        for(j in 1:2) {
            rep.str = ifelse(j==2, "-2", "")
            mutant <- input.df[,paste0(i, rep.str, colname)]
            wt <- input.df[,paste0("WT", rep.str, colname)]
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
    
    ggsave(file.path(output.dir, paste0("Enrichment of ", label, ".pdf")))
    
}

wt.enrichment(antisense.ratio[keep,], ".Antisense.Ratio", "antisense ratios", 0.4)
wt.enrichment(three.on.five.ratio[keep,], ".Three.On.Five.Ratio", "three prime to Five prime ratios", 1.4)

stats.df = NULL
for(metric in c("Antisense.Enrichment", "Three.On.Five.Enrichment")) {
    for(mutant in c("CK2", "Spt6-SA")) {
        for(replicate.id in c(1, 2)) {
            stats = t.test(all.counts[, paste0(mutant, ifelse(replicate.id==2, "-2", ""), ".", metric)])
            new.stats = data.frame(Metric=metric,
                                  Mutant=mutant,
                                  Replicate=replicate.id,
                                  Mean=stats$estimate,
                                  CI.low=stats$conf.int[1],
                                  CI.high=stats$conf.int[2],
                                  Mean.Absolute=2^stats$estimate,
                                  CI.low.Absolute=2^stats$conf.int[1],
                                  CI.high.Absolute=2^stats$conf.int[2],
                                  stringsAsFactors=FALSE)
                                  
            if(is.null(stats.df)) {
                stats.df = new.stats
            } else {
                stats.df = rbind(stats.df, new.stats)
            }
        }
    }
}
write.table(stats.df, file=file.path(output.dir, "Statistics.txt"), row.names=FALSE, col.names=TRUE, sep="\t")



# Do metagene type plot.
# Build a data.frame describing the matrix files.
matrix.path = file.path(base.output.dir, "Metagene")
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

