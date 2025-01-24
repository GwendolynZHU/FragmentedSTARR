#!/programs/R-4.2.1-r9/bin/Rscript

# Filename: 2_DESeq_STARR.r
# Description: This script is used to perform DESeq2 analysis on STARR-seq data
# Author: Yutong Zhu (yz2676@cornell.edu)
# Date: 2025-01-13
# Dependencies: DESeq2, dplyr, edgeR, argparse

#################################################
suppressPackageStartupMessages({
    library(argparse)
    library(DESeq2)
    library(dplyr)
    library(edgeR)
    library("R.utils")
});
#################################################
##### Constants #####
defaults <- list(
    "WHG-STARR" = list(
        DNA_rep = c(1),
        RNA_rep = c(1,2,3),
        ctrl_file = "data/normalization/srt_WHG_exon_ctrl.txt"
    ),
    "deep-ATAC-STARR" = list(
        DNA_rep = c(1,2,3,4,5,6),
        RNA_rep = c(1,2,3,4),
        ctrl_file = "data/normalization/srt_deep-ATAC_exon_ctrl.txt"
    )
)


#### Helper functions ####
# Create parser object and add parser arguments
setup_parser <- function() {
    parser <- ArgumentParser()

    parser$add_argument("--output", 
                        required=TRUE, 
                        help="Output file directory, \
                        should be the same directory as the previous script")

    parser$add_argument("--design", 
                        nargs="+", 
                        required=TRUE, 
                        help="Name of the design(s): \
                        They should be analyzed in the same script \
                        (generate a single DESeq table)")

    parser$add_argument("--starr", 
                        choices=c("deep-ATAC-STARR", "WHG-STARR"), 
                        default="deep-ATAC-STARR", 
                        help="File source of STARR-seq data")

    parser$add_argument("--resolution", 
                        type="integer",
                        default=5,
                        help="Resolution of the data (same as previous script)")

    parser$add_argument("--DNA_rep", 
                        type="integer", 
                        nargs="+", 
                        default=c(1,2,3,4,5,6), 
                        help="Specify which replicates to use for DNA, \
                        the max rep for deep-ATAC-STARR is 6, \
                        must set to 1 for WHG-STARR")

    parser$add_argument("--RNA_rep", 
                        type="integer", 
                        nargs="+", 
                        default=c(1,2,3,4), 
                        help="Specify which replicates to use for RNA, \
                        the max replicate for deep-ATAC-STARR is 4, \
                        must set to at most 3 for WHG-STARR")

    parser$add_argument("--cutoff", 
                        type="integer", 
                        default=5, 
                        help="cutoff value")

    parser$add_argument("--cpm_filter", 
                        default=FALSE, 
                        help="whether to use cpm_filter (stringent)")

    args <- parser$parse_args()
    return(args)
}


# Helper function to read files
read_file <- function(design) {
    forward_file <- read.table(file.path(getwd(), args$output, paste(args$starr, str(args$resolution), sep = "_"), 
                                        design, "combined_counts_f.txt"))
    reverse_file <- read.table(file.path(getwd(), args$output, paste(args$starr, str(args$resolution), sep = "_"), 
                                        design, "combined_counts_r.txt"))

    rownames(forward_file) <- paste0(design, ":forward:", forward_file$V4)
    rownames(reverse_file) <- paste0(design, ":reverse:", reverse_file$V4)
    return (rbind(forward_file, reverse_file))
}


# Helper function to calculate CPM cutoff
calculate_cpm_cutoff <- function(data, DNA_columns, cutoff) {
    cpm <- cpm(counts(dds), log=FALSE)
    lib_size <- mean(colSums(counts(dds))[,DNA_columns])
    message(paste("Library size (mean):", lib_size))
    cpm_cutoff <- args$cutoff/(lib_size/10^6)
    message(paste("CPM filter cutoff:", cpm_cutoff))

    return(list(cpm = cpm, cpm_cutoff = cpm_cutoff))
}


# Helper function to apply CPM filter
apply_cpm_filter <- function(cpm, DNA_columns, cpm_cutoff, min_libs = 4) {
    if(length(DNA_columns) == 1){
        keep <- cpm[, DNA_columns] >= cpm_cutoff
    } else {
        keep <- rowSums(cpm[, DNA_columns] >= cpm_cutoff) >= min_libs
    }
    return(keep)
}


# Helper function to apply raw count filter
apply_raw_count_filter <- function(dds, DNA_columns, cutoff) {
    if (length(DNA_columns) == 1){
        keep <- (rowSums(counts(dds)) >= cutoff) & (counts(dds)[, DNA_columns] >= 1)
    } else {
        keep <- (rowSums(counts(dds)) >= cutoff) & rowSums(counts(dds)[, DNA_columns] >= 1) > 0
    }
    return(keep)
}


# Helper function to save the z-score file
save_zscore_file <- function(result, control_indexes, deseq_directory) {
    mean_neg_ctrl <- mean(result[control_indexes,]$log2FoldChange, na.rm=TRUE)
    sd_neg_ctrl <- sd(result[control_indexes,]$log2FoldChange, na.rm=TRUE)
    result$z_score <- (result$log2FoldChange - mean_neg_ctrl) / sd_neg_ctrl
    logFC_cutoff <- 1.64 * sd_neg_ctrl + mean_neg_ctrl

    lfc_file <- file.path(deseq_directory, "LogFC_cutoff.txt")
    message(paste("LogFC cutoff:", logFC_cutoff))
    writeLines(paste(logFC_cutoff), lfc_file)

    return(list(result = result, logFC_cutoff = logFC_cutoff))
}


#### Processing functions ####
# DESeq2 with negative control provided
DESeq_with_ctrl <- function(data, args, deseq_directory) {
    # Append negative control
    neg_ctrl <- read.table(args$ctrl_file)
    rownames(neg_ctrl) <- paste0("neg_ctrl_", 1:nrow(neg_ctrl))
    all <- rbind(data, neg_ctrl)

    # Remove unnecessary columns
    count_table <- all[, -c(1,2,3,4)]

    # Choose the replicates
    colnames(count_table) <- c(
        paste0("DNA", args$original_dna_rep),
        paste0("RNA", args$original_rna_rep)
    )
    colnames(neg_ctrl) <- c(
        paste0("DNA", args$original_dna_rep),
        paste0("RNA", args$original_rna_rep)
    )
    DNA_columns <- paste0("DNA", args$DNA_rep)
    RNA_columns <- paste0("RNA", args$RNA_rep)

    all_data <- dplyr::select(count_table, all_of(DNA_columns), all_of(RNA_columns))
    ctrl_data <- dplyr::select(neg_ctrl, all_of(DNA_columns), all_of(RNA_columns))

    full_data <- all_data[1:(nrow(all_data)-nrow(ctrl_data)),]

    # Prepare count matrix (cts) and coldata
    cts <- all_data

    coldata <- data.frame(row.names=colnames(cts))
    coldata$condition <- c(rep("DNA", length(args$DNA_rep)), rep("RNA", length(args$RNA_rep)))
    coldata$condition <- factor(coldata$condition)

    stopifnot(all(rownames(coldata) %in% colnames(cts)))
    stopifnot(all(rownames(coldata) == colnames(cts)))

    dds <- DESeqDataSetFromMatrix(countData = cts, 
                                  colData = coldata, 
                                  design= ~ condition)
    
    # Filter data
    if (args$cpm_filter) {
        cpm_info <- calculate_cpm_cutoff(dds, DNA_columns, args$cutoff)
        keep <- apply_cpm_filter(cpm_info$cpm, DNA_columns, cpm_info$cpm_cutoff)
    } else {
        keep <- apply_raw_count_filter(dds, DNA_columns, args$cutoff)
    }

    dds_filtered <- dds[keep,]

    # Perform DESeq2 analysis
    control_genes <- rownames(dds_filtered) %in% rownames(ctrl_data)
    control_indexes <- which(control_genes)

    filtered_full_data <- rownames(dds_filtered) %in% rownames(full_data)
    filtered_full_indexes <- length(which(filtered_full_data))

    dds_filtered <- estimateSizeFactors(dds_filtered, controlGenes=control_indexes)#, type="poscounts")
    dds_filtered <- estimateDispersions(dds_filtered)
    dds_filtered <- nbinomWaldTest(dds_filtered)

    resLFC <- results(dds_filtered, contrast=c("condition", "RNA", "DNA"))

    if (args$starr == "deep-ATAC-STARR") {
        resLFC <- lfcShrink(dds_filtered, coef="condition_RNA_vs_DNA", type="apeglm")
    }

    # Save the results
    resLFC_results <- save_zscore_file(resLFC, control_indexes, deseq_directory)
    output <- resLFC_results$result[1:filtered_full_indexes,]
    write.table(output, 
                file.path(deseq_directory, "DE_results_nctrl.txt"), 
                sep="\t", row.names=TRUE, col.names=TRUE)
}


##### Main function #####
main <- function(args) {
    deseq_directory <- file.path(args$output, paste(args$starr, str(args$resolution), sep = "_"), "DESeq2")
    if (!dir.exists(deseq_directory)) {
        dir.create(deseq_directory)
    }

    # Set the parameters
    if (args$starr %in% names(defaults)) {
        args$original_dna_rep <- defaults[[args$starr]]$DNA_rep
        args$original_rna_rep <- defaults[[args$starr]]$RNA_rep
        args$ctrl_file <- defaults[[args$starr]]$ctrl_file
    } else {
        stop("Invalid STARR-seq data source")
    }

    if (args$starr == "WHG-STARR") {
        args$DNA_rep <- c(1)
        args$RNA_rep <- c(1,2,3)
    }

    # Generate DESeq2 table (including negative control)
    if (file.exists(args$ctrl_file)) {
        message("Generating DESeq2 table ... ")
    } else if (args$starr == "WHG-STARR") {
        message("Control file not found. Decompressing required files ...")
        gunzip(file.path("data", "normalization", "srt_WHG_exon_ctrl.txt.gz"), remove=FALSE)
    } else {
        stop("Check to make sure the negative control files exists at \
            data/normalization/")
    }
    # Read in the data
    design_count_table <- lapply(args$design, function(design) {
        read_file(design)
    })
    # print(dim(design_count_table[[1]]))
    # print(design_count_table[[1]][1:5,])
    combined_count_table <- do.call(rbind, design_count_table)

    # Perform DESeq2 analysis
    DESeq_with_ctrl(combined_count_table, args, deseq_directory)
}


# Execute the script
args <- setup_parser()
main(args)