# Custom function from yarn package to download GTExV8

library(yarn)

# Load all the dependencies (of note, these all need to be installed)
library ("repmis")
library ("Rgraphviz")
library ("ggplot2")
library ("gridExtra")
library ("reshape")
library ("org.Hs.eg.db")
library ("NMF")
library ("GenomicFeatures")
library ("doBy")
library ("plyr")
library ("WriteXLS")
library ("ontologyIndex")
library ("PRROC")
library ("yarn")
library ("readr")
library ("data.table")
library ("randomForest")
library ("parallel")
library ("precrec")
library ("dplyr")
library ("caret")
library ("CAGEr")
library ("e1071")
library ("ranger")
library ("MLmetrics")
library ("OOBCurve")
library ("genefilter")
library ("GeneOverlap")
library ("ontologyPlot")
library("biomaRt")
library ("BiocGenerics")

# We need to set the downloadGTEx.fixed function
downloadGTExV8 <- function (type = "genes", file = NULL, ...)
{
    phenoFile  <- "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt"
    pheno2File <- "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt"
    geneFile   <- "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz"
    message("Downloading and reading files")
    pdFile <- tempfile("phenodat", fileext = ".txt")
    download.file(phenoFile, destfile = pdFile)
    pd <- read_tsv(pdFile)
    pd <- as.matrix(pd)
    rownames(pd) <- pd[, "SAMPID"]
    ids <- sapply(strsplit(pd[, "SAMPID"], "-"), function(i) paste(i[1:2],
        collapse = "-"))
    pd2File <- tempfile("phenodat2", fileext = ".txt")
    download.file(pheno2File, destfile = pd2File)
    pd2 <- read_tsv(pd2File)
    pd2 <- as.matrix(pd2)
    rownames(pd2) <- pd2[, "SUBJID"]
    pd2 <- pd2[which(rownames(pd2) %in% unique(ids)), ]
    pd2 <- pd2[match(ids, rownames(pd2)), ]
    rownames(pd2) <- colnames(counts)
    pdfinal <- AnnotatedDataFrame(data.frame(cbind(pd, pd2)))
    if (type == "genes") {
        countsFile <- tempfile("counts", fileext = ".gz")
        download.file(geneFile, destfile = countsFile)
        cnts <- suppressWarnings(read_tsv(countsFile, skip = 2))
        genes <- unlist(cnts[, 1])
        geneNames <- unlist(cnts[, 2])
        counts <- cnts[, -c(1:2)]
        counts <- as.matrix(counts)
        rownames(counts) <- genes
        for (i in 1:nrow(problems(cnts))) {
            counts[problems(cnts)$row[i], problems(cnts)$col[i]] <- 1e+05
        }
        throwAway <- which(rowSums(counts) == 0)
        counts <- counts[-throwAway, ]
        genes <- sub("\\..*", "", rownames(counts))
        #host <- "dec2013.archive.ensembl.org"
        #biomart <- "ENSEMBL_MART_ENSEMBL"
        #dataset <- "hsapiens_gene_ensembl"
        attributes <- c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
            "start_position", "end_position","strand", "gene_biotype")
    }
    message("Creating ExpressionSet")
    pdfinal <- pdfinal[match(colnames(counts), rownames(pdfinal)),
        ]
    es <- ExpressionSet(as.matrix(counts))
    phenoData(es) <- pdfinal
    #pData(es)["GTEX-YF7O-2326-101833-SM-5CVN9", "SMTS"] <- "Skin"
    #pData(es)["GTEX-YEC3-1426-101806-SM-5PNXX", "SMTS"] <- "Stomach"

# This step is not annotating everything correctly, so later, we need to annotate the data with
# the gencode annotation file (GRCh38/hg19, v26)
    message("Annotating from biomaRt")
    es <- annotateFromBiomart(obj = es, genes = genes, attributes = attributes)
    message("Cleaning up files")
    unlink(pdFile)
    unlink(pd2File)
    unlink(countsFile)
    if (!is.null(file))
        saveRDS(es, file = file)
    return(es)
}


#===========================================================================================================
# Actually downloading GTExV8
#===========================================================================================================
gtex8 <- downloadGTExV8(type="genes", file="gtex8.rds")
