# Extract the exons only from the annotation file

# Read the hg38 gtf file we downloaded in the previous step
hg38  <- read.table("../GRCh38/gencode.v26.annotation.gtf",sep="\t")
temp1 <- stringr::str_split(hg38[,9],pattern=";")

geneout <- stringr::str_replace(sapply(X = temp1, FUN = stringr::str_subset,
                                       pattern = "gene_id"),
                            pattern = "gene_id ", replacement = "")

genetype <- stringr::str_replace(sapply(X = temp1, FUN = stringr::str_subset,
                                        pattern = "gene_type"),
                                 pattern = " gene_type ", replacement = "")

genename <- stringr::str_replace(sapply(X = temp1, FUN = stringr::str_subset,
                                       pattern = "gene_name"),
                            pattern = "gene_name", replacement = "")

exon2    <- cbind(hg38[, 1:8],  gene_name = genename, gene_biotype = genetype, gene_id = geneout)

colnames(exon2)[1:8] <- c("chromsome_name", "annotation_source", "feature_type", "gene_start",
                         "gene_end", "score", "strand", "phase")

exon2            <-exon2[,c(1,4,5,7,2,3,9:11)]
gene_information <- dplyr::filter(exon2, feature_type=="gene")

readr::write_delim(gene_information, "../GRCh38/EXON_ANNOTATION", delim="\t")
