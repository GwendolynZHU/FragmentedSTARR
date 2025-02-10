#!/programs/R-4.2.1-r9/bin/Rscript

# Filename: 4_visualization.r
# Description: This script is used to visualize the results of DESeq2 analysis on STARR-seq data
# Author: Yutong Zhu (yz2676@cornell.edu)
# Date: 2025-02-05
# Dependencies: ggplot2, dplyr, argparse, wesanderson, systemfonts, showtext, ggtext

################################################################################
suppressPackageStartupMessages({
    if (!requireNamespace("systemfonts", quietly = TRUE)) {
        install.packages("systemfonts")
        install.packages("showtext")
        install.packages("ggtext")
    }
    
    library(argparse)
    library(ggplot2)
    library(dplyr)
    library(wesanderson)
    library(systemfonts)
    library(showtext)
    library(ggtext)

    font_add("HelveticaNeueFont", "/fs/cbsuhy01/storage/yz2676/package/fonts/HelveticaNeue/HelveticaNeue.ttf")
    showtext_auto()
});
################################################################################
##### Helper functions #####
# Create parser object and add parser arguments
setup_parser <- function() {
    parser <- ArgumentParser()

    parser$add_argument("--output", 
                        required=TRUE, 
                        help="Output file directory")

    parser$add_argument("--baseline",
                        required=TRUE,
                        help="The baseline design for comparison")

    parser$add_argument("--comparisons",
                        required=TRUE,
                        nargs="+",
                        help="The groups to compare")

    parser$add_argument("--double_sided",
                        action="store_true",
                        help="Whether to compare within groups or between groups.\
                        Default is between groups. If specified, the comparison \
                        will be within groups. (i.e. tata.down, tata.up)")

    parser$add_argument("--starr", 
                        choices=c("deep-ATAC-STARR", "WHG-STARR"), 
                        default="deep-ATAC-STARR", 
                        help="File source of STARR-seq data")

    parser$add_argument("--resolution", 
                        type="integer",
                        default=5,
                        help="Resolution of the data (same as previous scripts)")

    args <- parser$parse_args()
    return(args)
}


load_lfc <- function(deseq_dir, group) {
    lfc_file <- read.table(file.path(deseq_dir, paste0(group, "_lfc.txt")), 
                           header=TRUE, sep="\t")
    rownames(lfc_file) <- lfc_file$name
    lfc_file$name <- NULL
    return(lfc_file)
}


get_effect_size <- function(df, index) {
    rn <- rownames(df)
    es <- ifelse(index %in% rn, df[index, "log2FC"], NA)
}


# Plot the comparison within groups
plot_two_sides <- function(baseline, group1, group2, deseq_dir, visualization_dir) {
    baseline_file <- load_lfc(deseq_dir, baseline)
    group1_file <- load_lfc(deseq_dir, group1)
    group2_file <- load_lfc(deseq_dir, group2)

    full <- filter(baseline_file, enhancer=="True")
    full <- full[rownames(full) %in% union(rownames(group1_file), rownames(group2_file)), ]
    filtered_group1 <- group1_file[rownames(group1_file) %in% rownames(full), ]
    filtered_group2 <- group2_file[rownames(group2_file) %in% rownames(full), ]
    
    max_group1 <- filter(filtered_group1, prominent_tss=="maximum")
    min_group1 <- filter(filtered_group1, prominent_tss=="minimum")
    max_group2 <- filter(filtered_group2, prominent_tss=="maximum")
    min_group2 <- filter(filtered_group2, prominent_tss=="minimum")
    max_tss_delta <- rbind(max_group1, max_group2)
    min_tss_delta <- rbind(min_group1, min_group2)

    # Identify the common elements
    common12 <- intersect(rownames(full), rownames(min_tss_delta))
    common13 <- intersect(rownames(full), rownames(max_tss_delta))
    merge <- union(common12, common13)
    print(merge)

    # Compute the effect sizes for statistical tests
    effect_size12 <- get_effect_size(full, common12)
    effect_size21 <- get_effect_size(min_tss_delta, common12)
    effect_size13 <- get_effect_size(full, common13)
    effect_size31 <- get_effect_size(max_tss_delta, common13)

    # Compute effect sizes for plotting
    es_0 <- get_effect_size(full, merge)
    es_1 <- get_effect_size(min_tss_delta, merge)
    es_2 <- get_effect_size(max_tss_delta, merge)

    design_labels <- list(
        "tata" = list(name = "Entire Core Promoter", 
                      max_label = "\u0394Maximum_TSS", 
                      min_label = "\u0394Minimum_TSS"),
        "pause" = list(name = "Pause Region", 
                       max_label = "Maximum_TSS\u0394pause", 
                       min_label = "Minimum_TSS\u0394pause")
    )
    labels <- if (grepl("tata", group1, ignore.case = TRUE)
            ) design_labels[["tata"]] else design_labels[["pause"]]
    design <- labels$name

    # Prepare dataframe for plotting
    df <- data.frame("None" = es_0, es_1, es_2, index = seq_along(es_1))
    colnames(df)[2:3] <- c(labels$min_label, labels$max_label)
    df_long <- tidyr::pivot_longer(df, cols = c("None", labels$min_label, labels$max_label), names_to = "group", values_to = "values")
    df_long$group <- factor(df_long$group, levels = c(labels$min_label, "None", labels$max_label))

    # Perform paired t-tests
    result_1 <- t.test(effect_size12, effect_size21, alternative = "greater", paired = TRUE)
    result_2 <- t.test(effect_size13, effect_size31, alternative = "greater", paired = TRUE)

    p <- ggplot(df_long, aes(x = group, y = values, color = group, group = index)) +
        geom_point(size = 5) +
        geom_line(color = "grey") +
        xlab("") +
        ylab("Activity (log2)") +
        ggtitle(paste0(design, " Loss")) +
        scale_y_continuous(limits = c(-1, 7), breaks = seq(0, 6, by = 2)) +
        scale_color_manual(values = c("#98C1D6", "black", "#D89A88")) +
        annotate("text", x = 1.1, y = 6.4, label = paste("p =", format(result_1$p.value, digits = 2)), color = "black", size = 6, hjust = 0) +
        annotate("text", x = 2.1, y = 6.4, label = paste("p =", format(result_2$p.value, digits = 2)), color = "black", size = 6, hjust = 0) +
        annotate("text", x = 1.25, y = 6, label = paste("n =", length(effect_size12)), color = "black", size = 6, hjust = 0) +
        annotate("text", x = 2.25, y = 6, label = paste("n =", length(effect_size13)), color = "black", size = 6, hjust = 0) +
        geom_hline(yintercept = 2.236, col = "red", lwd = 0.5, lty = 2) +
        theme(
            panel.background = element_rect(fill = "white"),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "none",
            axis.line = element_line(color = "black"),
            axis.line.x.top = element_blank(),
            axis.line.y.right = element_blank(), 
            axis.title = element_text(size = 20, color = "black"),
            axis.text = element_text(size = 17, color = "black"),
            plot.title = element_text(hjust = 0.5, size = 20, color = "black")
        )

    ggsave(file.path(visualization_dir, paste0(design, ".pdf")), plot = p, width = 7, height = 7, dpi = 300)
    message("Finished saving the plot for pairwise comparisons for ", design, ".")
}


# Plot the comparison between groups
plot_one_side <- function(baseline, group, deseq_dir, visualization_dir) {
    stop("NotImplementedError: 'plot_one_side' has not been implemented yet.")
}

################################################################################
##### Main function #####
main <- function(args) {
    visualization_dir <- file.path(args$output, paste0(args$starr, "_", 
                            as.character(args$resolution)), "Visualization")

    if (!dir.exists(visualization_dir)) {
        dir.create(visualization_dir)
    }

    deseq_dir <- file.path(args$output, paste0(args$starr, "_", 
                            as.character(args$resolution)), "DESeq2")

    if (args$double_sided) {
        if (length(args$comparisons) %% 2 != 0) {
            stop("Number of comparison groups must be even for double-sided mode.")
        }
        for (i in seq(1, length(args$comparisons), by=2)) {
            plot_two_sides(args$baseline, args$comparisons[i], args$comparisons[i+1], 
                            deseq_dir, visualization_dir)
        }
    } else {
        for (comparison in args$comparisons) {
            plot_one_side(args$baseline, comparison, 
                            deseq_dir, visualization_dir)
        }
    }
}

################################################################################
# Execute the script
args <- setup_parser()
main(args)