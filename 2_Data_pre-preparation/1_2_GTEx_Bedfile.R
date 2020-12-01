# Save the GTEx data as bedfile

library(readr)
dat <- read_csv("../data/GTExV8_specificity.csv")

write_delim(dat[, c("chromosome_name", "start_position", "end_position", "gene_id", "gene_name")], "../data/GTExV8.bed", delim="\t")
