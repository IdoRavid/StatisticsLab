if (!require('data.table')) { install.packages('data.table'); library('data.table') }
# Ex2 function to count reads
count_reads_table <- function(reads, beg_region, end_region) {
  locs <- reads$Loc[reads$Loc >= beg_region & reads$Loc <= end_region]
  counts <- tabulate(locs - beg_region + 1, nbins = end_region - beg_region + 1)
  counts
}
