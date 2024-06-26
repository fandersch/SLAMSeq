#!/usr/bin/env Rscript

library(getopt)
library(data.table)
library(tidyverse)

spec = matrix(c(
	'counts',           'c', 1, "character",   "Table of sample counts.",
	'factors',          'f', 1, "character",   "Table of sample scaling factors,",
	'nid',              'i', 1, "integer",     "Number of row-ID columns at the left of the table (0). All of the will be kept.",
	'outpref',          'o', 1, "character",   "Destination file."
), byrow=TRUE, ncol=5)

opt <- getopt(spec)
# opt <- list(counts='/Volumes/groups/zuber/USERS/florian.andersch/data/workspace/SLAMseq/2023_06_13_Martina/slamdunk/counts_collapsed/counts.txt', factors='/Volumes/groups/zuber/USERS/florian.andersch/data/workspace/SLAMseq/2023_06_13_Martina/slamdunk/spike_summary_counts.factors.txt', nid=3, outpref='/Volumes/groups/zuber/USERS/florian.andersch/data/workspace/SLAMseq/2023_06_13_Martina/slamdunk/counts_collapsed/counts')


if (is.null(opt$nid)) {
	opt$nid <- 0
}

outscaled <- paste0(opt$outpref, '_scaled.txt')
outrest <- paste0(opt$outpref, '_rest.txt')


# Input
DT <- fread(opt$counts)
factors <- fread(opt$factors)
setnames(factors, c('row_ID', 'factor'))

#only consider sample-id
colnames_original <- colnames(DT)
colnames(DT) <- colnames_original %>% str_split("_") %>% lapply('[[', 1) %>% unlist
factors$row_ID <- factors$row_ID %>% str_split("_") %>% lapply('[[', 1) %>% unlist

# Save columns for which no factor is provided
fwrite(DT[, c(1:opt$nid, which(!names(DT) %in% factors$row_ID)), with=FALSE], file=outrest, quote=FALSE, sep="\t")

# Select count columns for which a factor is provided.
DT <- DT[, c(1:opt$nid, which(names(DT) %in% factors$row_ID)), with=FALSE]

# Apply scaling
for (id in factors[[1]]) {
	set(DT, i=NULL, j=id, DT[[id]] * factors[row_ID==id, factor])
}

#restore original column header
colnames(DT) <- colnames_original

fwrite(DT, file=outscaled, quote=FALSE, sep="\t")


