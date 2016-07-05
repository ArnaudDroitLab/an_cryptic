library(rtracklayer)
library(metagene)


bam.files = list.files("output/bwa", pattern=".*sense_sorted.bam$", full.names=TRUE)
regions = import("input/exons.cut.bed", format="BED")

# Filter regions < 100 nt
regions.above.100 = regions[as.data.frame(regions)$width > 100]


mg <- metagene$new(regions = regions.above.100, bam_files = bam.files)
save(mg, "mg.RData")