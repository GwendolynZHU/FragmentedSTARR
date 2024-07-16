#!/usr/bin/env /programs/R-4.2.1-r9/bin/Rscript

# Sample command to run the script: 
# /programs/R-4.2.1-r9/bin/Rscript 2_DESeq_STARR.r -o /path_to_file/EnhancerNet --name 5p, 3p
# /programs/R-4.2.1-r9/bin/Rscript 2_DESeq_STARR.r -o ../new_data --name TSS_p TSS_n TSS_b --RNA_rep 1 2 3 4 --cutoff 10 --starr deep_ATAC_STARR --either TRUE

write_params <- function(args, file) {
  conn <- file(file, "w")
  on.exit(close(conn))
  
  lapply(names(args), function(arg) {
    write(paste(arg, args[[arg]]), conn)
  })
}


def_colnames <- function(dataframe, name, suffixes){
  base_suffixes <- c("chr","start","end")
  dynamic_names <- sapply(suffixes, function(suffix) paste0(name, "_", suffix))
  new_colnames <- c(base_suffixes, dynamic_names)
  
  colnames(dataframe) = new_colnames
  return(dataframe)
}


# Function to bind each pair with control_gene
bind_each_pair <- function(name_data) {
  # Binding the forward, reverse data frames
  bind_rows(name_data$forward, name_data$reverse)
}


DESeq_regular <- function(ctrl, DNA_rep, RNA_rep, output_file){
  ctrl <- ctrl[which(rowSums(ctrl)>10),]
  
  cts <- ctrl
  rownames(cts) <- rownames(ctrl)
  
  coldata <- data.frame(colnames(cts))
  coldata$condition <- c(rep("DNA", length(DNA_rep)),rep("RNA", length(RNA_rep)))
  
  rownames(coldata) <- coldata$colnames.cts.
  coldata[,1] <- NULL
  coldata$condition <- factor(coldata$condition)
  
  all(rownames(coldata) %in% colnames(cts))
  all(rownames(coldata) == colnames(cts))
  
  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)
  dds <- DESeq(dds)
  res <- results(dds, contrast=c("condition", "RNA", "DNA"))
  resLFC <- lfcShrink(dds, coef="condition_RNA_vs_DNA", type="apeglm")
  
  write.table(res, file=output_file,
              append = FALSE, sep = "\t",
              row.names = TRUE, col.names = TRUE)
}


DESeq_with_ctrl <- function(all, ctrl, outdir, DNA_rep, RNA_rep, cutoff, cpm_filter, starr){
    all <- all[, -(1:3)]
    ctrl <- ctrl[, -(1:3)]
    if (starr == "deep_ATAC_STARR") {
        ori_DNA_columns <- paste0("DNA_", c(1,2,3,4,5,6))
        ori_RNA_columns <- paste0("RNA_", c(1,2,3,4,5,6,7))
    }
    else {
        ori_DNA_columns <- paste0("DNA_", c(1))
        ori_RNA_columns <- paste0("RNA_", c(1,2,3))
    }
    colnames(all) <- c(ori_DNA_columns, ori_RNA_columns)
    colnames(ctrl) <- c(ori_DNA_columns, ori_RNA_columns)
    

    ### add in flexibility to allow a subset of replicates
    DNA_columns <- paste0("DNA_", DNA_rep)
    RNA_columns <- paste0("RNA_", RNA_rep)
    all <- dplyr::select(all, all_of(DNA_columns), all_of(RNA_columns))
    ctrl <- dplyr::select(ctrl, all_of(DNA_columns), all_of(RNA_columns))
    full <- all[1:(nrow(all)-nrow(ctrl)),]
    
    cts <- all
    rownames(cts) <- rownames(all)
    
    coldata <- data.frame(colnames(cts))
    coldata$condition <- c(rep("DNA", length(DNA_rep)),rep("RNA", length(RNA_rep)))
    
    rownames(coldata) <- coldata$colnames.cts.
    coldata[,1] <- NULL
    coldata$condition <- factor(coldata$condition)
    
    ## reset cts row index
    #row.names(cts) <- 1:nrow(cts)
    
    all(rownames(coldata) %in% colnames(cts))
    all(rownames(coldata) == colnames(cts))
    
    dds <- DESeqDataSetFromMatrix(countData = cts,
                                    colData = coldata,
                                    design = ~ condition)
    # dds

    if (cpm_filter) { # at least 4 DNA libs >= cpm cutoff
        cpm = cpm(counts(dds), log=FALSE)
        lib_size = mean(colSums(counts(dds))[1:length(DNA_rep)])
        message(paste("mean_lib_size", lib_size))
        cpm_cutoff = cutoff/(lib_size/10^6)
        message(paste("CPM_filter_cutoff", cpm_cutoff))
        if(length(DNA_rep) == 1){
            keep = cpm[, DNA_columns] >= cpm_cutoff
        } else {
            keep = rowSums(cpm[, DNA_columns] >= cpm_cutoff) >= 4
        }
    } else { #raw counts filter
        if (length(DNA_rep) == 1){
            # print(head(counts(dds)[, DNA_columns] >= 1))
            keep <- (rowSums(counts(dds)) >= cutoff) & (counts(dds)[, DNA_columns] >= 1)
        }
        else {keep <- (rowSums(counts(dds)) >= cutoff) & rowSums(counts(dds)[, DNA_columns] >= 1) > 0}
    }

    dds_filtered <- dds[keep,]
    control_genes <- rownames(dds_filtered) %in% rownames(ctrl)
    control_indexes <- which(control_genes)

    filtered_full <- rownames(dds_filtered) %in% rownames(full)
    filtered_full_indexes <- length(which(filtered_full))

    dds_filtered <- estimateSizeFactors(dds_filtered, controlGenes=control_indexes)#, type="poscounts")
    dds_filtered <- estimateDispersions(dds_filtered)
    dds_filtered <- nbinomWaldTest(dds_filtered)

    # Check if the directory exists
    qc_directory <- file.path(outdir, "plots", "dispersion")
    if (!file.exists(qc_directory)) {
    dir.create(qc_directory, recursive = TRUE)
    }

    png(filename = file.path(qc_directory, paste0("nctrl_", cutoff, "_dispersion_plot.png")))
    plotDispEsts(dds_filtered)
    dev.off()

    #   resultsNames(dds_filtered)
    resLFC <- results(dds_filtered, contrast=c("condition", "RNA", "DNA"))
    
    if (starr == "deep_ATAC_STARR") {resLFC <- lfcShrink(dds_filtered, coef="condition_RNA_vs_DNA", type="apeglm")}
    # resLFC
  
    ## first option: LFC>=1 & padj<0.1
    # output <- resLFC[1:filtered_full_indexes, ][resLFC[1:filtered_full_indexes, ]$padj < 0.1 & 
    #                                             resLFC[1:filtered_full_indexes, ]$log2FoldChange >= 1, ]
    ## second option: z-score > 1.64 - 95th percentile for a one-sided Gaussian distribution
    mean_neg_ctrl <- mean(resLFC[control_indexes,]$log2FoldChange)
    std_neg_ctrl = sd(resLFC[control_indexes,]$log2FoldChange)

    resLFC$"z-score" = (resLFC$log2FoldChange-mean_neg_ctrl)/std_neg_ctrl
    logFC_cutoff = 1.64*std_neg_ctrl+mean_neg_ctrl
    lfc_file = file.path(outdir, "LogFC_cutoff.txt")
    if (!file.exists(lfc_file)) {
        message("LFC: ", logFC_cutoff)
        conn <- file(lfc_file, "w")
        write(logFC_cutoff, conn)
        close(conn)  
        }

    output <- resLFC[1:filtered_full_indexes,]
    write.table(output, file=file.path(outdir, paste0("DE_results_nctrl.txt")),
                append = FALSE, sep = "\t", 
                row.names = TRUE, col.names = TRUE)
    message("DESeq results saved.")
    
    filtered_full <- resLFC[grepl("^(full|either)", rownames(resLFC)), ] 
    if (starr == "deep_ATAC_STARR") {
        full_output <- filtered_full[filtered_full$padj < 0.1 & filtered_full$log2FoldChange >= logFC_cutoff, ]
        } else {
            full_output <- filtered_full[filtered_full$log2FoldChange >= logFC_cutoff, ]
        }
    write.table(full_output, file=file.path(outdir, paste0("DE_results_full.txt")),
                append = FALSE, sep = "\t", 
                row.names = TRUE, col.names = TRUE)
}


save_enhancers <- function(deseq, outdir, design, starr, either) {
    ### find the row index for full_f elements
    f_idx <- length(which(substring(rownames(deseq), 1, nchar(design)+2) == paste0(design,"_f")))
    r_idx <- length(which(substring(rownames(deseq), 1, nchar(design)+2)== paste0(design,"_r")))

    ### save the elements for each subgroup
    full_f = substring(rownames(deseq)[1:f_idx], first=nchar(design)+3)
    full_r = substring(rownames(deseq)[(f_idx+1):(f_idx+r_idx)], first=nchar(design)+3)

    common = intersect(full_f, full_r)
    # print(common)
    
    ### need to find the original element
    ori_f <- read.table(file.path(outdir, starr, design, paste0("srt_", design, "_f.bed")))
    ori_r <- read.table(file.path(outdir, starr, design, paste0("srt_", design, "_r.bed")))
    
    rownames(ori_f) <- paste0(design, "_f", 1:nrow(ori_f))
    rownames(ori_r) <- paste0(design, "_r", 1:nrow(ori_r))
    
    names_f <- paste0(design, "_f", common)
    names_r <- paste0(design, "_r", common)
    
    # Filter dataframe based on the numbers in names
    filtered_rows_f <- ori_f[rownames(ori_f) %in% names_f, ]
    filtered_rows_r <- ori_r[rownames(ori_r) %in% names_r, ]
    
    # re-indexing the dataframes so that they can be added together
    rownames(filtered_rows_r) <- rownames(filtered_rows_f)
    
    result_values <- filtered_rows_f[,4:ncol(filtered_rows_f)] + filtered_rows_r[,4:ncol(filtered_rows_r)]
    result <- cbind(filtered_rows_f[,1:3], result_values)
  
    ## add in either-orientation elements
    if (either == "either_ori"){
        either_f_idx <- length(which(substring(rownames(deseq), 1, nchar("either")+2) == "either_f"))
        either_f = substring(rownames(deseq)[(f_idx+r_idx+1):(f_idx+r_idx+either_f_idx)], first=nchar("either")+3)
        either_r = substring(rownames(deseq)[(f_idx+r_idx+either_f_idx+1):nrow(deseq)], first=nchar("either")+3)

        ori_either_f <- read.table(file.path(outdir, starr, design, "unpaired", "srt_full_either_f.bed"))
        ori_either_r <- read.table(file.path(outdir, starr, design, "unpaired", "srt_full_either_r.bed"))
        
        filtered_either_rows_f <- ori_either_f[rownames(ori_either_f) %in% either_f, ]
        filtered_either_rows_r <- ori_either_r[rownames(ori_either_r) %in% either_r, ]

        result <- bind_rows(result, filtered_either_rows_f, filtered_either_rows_r)
    }
    
    write.table(result, file=file.path(outdir, starr, design, paste0("srt_", design, "_", either, "_e.bed")),
                append = FALSE, sep = "\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)

    
    if (nrow(result)>=1) {
        message("Enhancer elements successfully identified and saved.")
    }
    else {
        message("No enhancer elements captured.")
    }
}


read_in_bedfiles <- function(outdir, starr, either, name, orientation="f") {
    if (name == "full" && either){
        file <- read.table(file.path(outdir, starr, name, "unpaired", paste0("srt_", name, "_either_", orientation, ".bed")))
        rownames(file) <- paste0("either_", orientation, 1:nrow(file))
    } else {
        file <- read.table(file.path(outdir, starr, name, paste0("srt_", name, "_", orientation, ".bed")))
        rownames(file) <- paste0(name, "_", orientation, 1:nrow(file))
    }
    return(file)
}
################################################################################
suppressPackageStartupMessages({
  library(argparse)
  library(DESeq2)
  library(dplyr)
  library(edgeR)
});

# create parser object and add parser arguments
parser <- ArgumentParser()

parser$add_argument("-o", "--output", required=T, help="Output file directory")
parser$add_argument("--DNA_rep", type="integer", nargs="+", default=c(1,2,3,4,5,6), help="number of replicate")
parser$add_argument("--RNA_rep", type="integer", nargs="+", default=c(1,2,3,4,5,6,7), help="number of replicate")
parser$add_argument("-c", "--cutoff", type="integer", default=10, help="cutoff value")
parser$add_argument("--cpm_filter", type="logical", default=FALSE, help="whether to use cpm_filter (stringent)")
parser$add_argument("--name", nargs="+", required=T, help="design of element, e.g. TSS_b, TSS_p, TSS_n, pause_site_b, ...")
parser$add_argument("--starr", default="deep_ATAC_STARR", help="deep_ATAC_STARR or WHG_STARR")
parser$add_argument("--either", default=FALSE, help="Whether to include either-orientation enhancers")
parser$add_argument("--neg_ctrl", default=TRUE, help="Whether to make a DE file containing all elements")

args <- parser$parse_args()

################################################################################
# Define the directory path
orientation_dir <- ifelse(args$either, "either_ori", "ori_ind")
deseq_directory <- file.path(args$output, "DESeq", args$starr, orientation_dir)

# Check if the directory exists
if (!file.exists(deseq_directory)) {
dir.create(deseq_directory, recursive = TRUE)
}

### Global files & parameters
if (args$starr == "WHG_STARR") {
    if (length(args$DNA_rep) == 6) {
            args$DNA_rep <- c(1)
    }
    if (length(args$RNA_rep) == 7){
            args$RNA_rep <- c(1,2,3)
    }
}
# Save parameters:
write_params(args, file.path(deseq_directory, "params.txt"))

ctrl_file <- ifelse(args$starr == "deep_ATAC_STARR", "srt_deep_ATAC_exon_ctrl.bed", "srt_WHG_exon_ctrl.txt")
control_gene <- read.table(file.path("/fs","cbsuhy01","storage","yz2676","data","STARR-seq","normalization", ctrl_file))
neg_ctrl_file <- ifelse(args$starr == "deep_ATAC_STARR", "DESeq_ATAC_ctrl.txt", "DESeq_WHG_ctrl.txt")
neg_ctrl <- file.path("/fs","cbsuhy01","storage","yz2676","data","STARR-seq","normalization", neg_ctrl_file)

# Validate that the negative ctrl is normally distributed
if (!file.exists(neg_ctrl)){
  message("Testing negative controls ... ")
  ctrl <- control_gene
  ctrl <- ctrl[, -(1:4)]
  DESeq_regular(ctrl, args$DNA_rep, args$RNA_rep, neg_ctrl)
  message("Negative ctrl results saved.")
}

ctrl <- control_gene[, -(4)]
rownames(ctrl) <- paste0("ctrl", 1:nrow(ctrl))

# Generate a combined DE results - including controls
if (args$neg_ctrl){
    message("Generating negative control DE results ... ")

    partial_list <- list()
    for (n in args$name) {
        partial_f <- read_in_bedfiles(args$output, args$starr, args$either, n, "f")
        partial_r <- read_in_bedfiles(args$output, args$starr, args$either, n, "r")
        partial_list[[paste0(n, "_f")]] <- partial_f
        partial_list[[paste0(n, "_r")]] <- partial_r
    }

    full_f <- read_in_bedfiles(args$output, args$starr, FALSE, "full", "f")
    full_r <- read_in_bedfiles(args$output, args$starr, FALSE, "full", "r")
    colnames(ctrl) <- colnames(full_f)

    if (args$either) {
        either_f <- read_in_bedfiles(args$output, args$starr, TRUE, "full", "f")
        either_r <- read_in_bedfiles(args$output, args$starr, TRUE, "full", "r")
        all <- bind_rows(c(list(full_f), list(full_r), list(either_f), list(either_r), partial_list, list(ctrl)))
    } else {
        all <- bind_rows(c(list(full_f), list(full_r), partial_list, list(ctrl)))
    }
    DESeq_with_ctrl(all, ctrl, deseq_directory, args$DNA_rep, args$RNA_rep, args$cutoff, args$cpm_filter, args$starr)
}


# Output full enhancer elements with orientation-independence
# if (!file.exists(file.path(args$output, args$starr, "full", paste0("srt_full_",orientation_dir,"_e.bed")))){
#   message("Started searching for enhancers ... ")
  
#   full <- read.table(file.path(deseq_directory, "DE_results_full.txt"))
#   save_enhancers(full, args$output, "full", args$starr, orientation_dir)
# }