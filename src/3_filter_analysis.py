import os
import argparse

import pandas as pd
import numpy as np

import helpers

class ExtractEnhancer:
    """
    Extract enhancers given the DESeq2 results.
    """
    def __init__(self):
        self.svr_dpr = "src/SVRb_DPR.model"
        self.svr_tata = "src/SVRtata_TATAbox.model"
        self.fasta = "data/reference/hg38.fa"
    
    @staticmethod
    def add_arguments(parser):
        """
        Add arguments to the parser.
        """
        parser.add_argument("--outdir", required=True, 
                            help="Path to the output directory.")
        
        parser.add_argument("--starr", choices=["WHG-STARR", "deep-ATAC-STARR"], 
                            required=True, help="File source")
        
        parser.add_argument("--design", required=True, 
                            nargs="+", help="Name of design(s)")
        
        parser.add_argument("--resolution", default=5, type=int,
                            help="Resolution of the data.")
        
        parser.add_argument("--either", default=False, action="store_true",
                            help="Whether to include either-orientation enhancers.")
        
        parser.add_argument("--padj", default=0.05, type=float,
                            help="Adjusted p-value cutoff.")
        
        parser.add_argument("--logfc", default="avg", choices=["avg", "max"],
                            help="How to combine the activities of different orientations.")
        
        parser.add_argument("--calc_active_rate", default=False, action="store_true",
                            help="Whether to calculate the active rate of enhancers.")

        parser.add_argument("--annotate_extra", default=False, action="store_true",
                            help="Whether to annotate extra elements information.")
        
        parser.add_argument("--original_bed", nargs="+",
                            help="Path to the original bed file(s). \
                            Must be in the same order as the design names. \
                            Required if --annotate_extra or --core_promoter is True.")
        
        parser.add_argument("--annotate_bigwig", default=False, action="store_true",
                            help="Whether to annotate prominent TSSs.")
        
        parser.add_argument("--bigwig", default=None, type=str, nargs="+",
                            help="Path to the bigwig files. \
                            Must be in the same order of plus, minus strand, \
                            Required if --annotate_bigwig is True.")

        parser.add_argument("--core_promoter", default=False, action="store_true",
                            help="Whether to annotate core promoters.")

    def annotate_enhancers(self, args):
        """
        Save the annotated enhancers.
        """
        deseq_dir = f"{args.outdir}/{args.starr}_{args.resolution}/DESeq2"
        deseq_results = pd.read_csv(deseq_dir+"/DE_results_nctrl.txt", 
                                    sep="\t", index_col=0, header=0)
        indexes = deseq_results.index.to_series().str.split(":", expand=True)
        
        deseq_results[["design", "orientation", "name"]] = indexes
        deseq_results.reset_index(drop=True, inplace=True)

        with open(deseq_dir+"/LogFC_cutoff.txt", "r") as f:
            logfc_cutoff = round(float(f.readline().strip()), 3)
        
        for design in args.design:
            design_deseq = deseq_results[deseq_results["design"] == design]
            design_dir = f"{args.outdir}/{args.starr}_{args.resolution}/{design}"
            extract_enhancers = self.extract_enhancers(design_deseq, logfc_cutoff, 
                                                       design_dir, args.either, 
                                                       args.padj, args.logfc)
            if not extract_enhancers.empty:
                if args.original_bed:
                    original_bed_design = args.original_bed[args.design.index(design)]
                    if "full" in design.lower():
                        full_original = original_bed_design

                if args.annotate_extra and not args.original_bed:
                    raise ValueError("Please provide the path to the original bed file(s).")
                elif args.annotate_extra:
                    extract_enhancers = self.annotate_extra_columns(
                                            extract_enhancers, original_bed_design)
                
                if args.calc_active_rate:
                    active_rate_file = f"{deseq_dir}/{design}_active_rate.txt"
                    self.calculate_active_rate(design, extract_enhancers, 
                                               active_rate_file)
                
                if args.annotate_bigwig and not args.bigwig:
                    raise ValueError("Please provide the path to the bigwig files.")
                elif args.annotate_bigwig:
                    extract_enhancers = self.annotate_prominent_tss(
                            design, extract_enhancers, args.bigwig, full_original)
                
                if args.core_promoter and not args.original_bed:
                    raise ValueError("Please provide the path to the original bed file(s).")
                elif args.core_promoter:
                    extract_enhancers = self.annotate_core_promoters(
                        design, extract_enhancers, original_bed_design)

                enhancer_file_path = f"{deseq_dir}/{design}_lfc.txt"
                extract_enhancers.to_csv(enhancer_file_path, sep="\t", 
                                         index=False, header=True)
            else:
                print(f"No enhancers found for {design}.")

    def extract_enhancers(self, deseq_table: pd.DataFrame, logfc_cutoff: float, 
                          design_dir: str, either: bool, padj: float, lfc_method: str
                          ) -> pd.DataFrame:
        """
        Extract enhancers given the DESeq2 results.

        Args:
            deseq_table (pd.DataFrame): DESeq2 results.
            logfc_cutoff (float): LogFC cutoff.
            design_dir (str): Name of the design directory.
            either (bool): Whether to include either-orientation enhancers.
                Both orientation enhancers: logFC > cutoff, padj < padj_cutoff in both forward and reverse.
                Either orientation enhancers: only mapped in one direction 
                        (elements that only have forward or reverse orientation data available)
                                        and logFC > cutoff, padj < padj_cutoff.
            padj (float): Adjusted p-value cutoff.
            lfc_method (str): How to combine the activities of different orientations. 
                Default is taking averages.
        """
        forward_deseq = deseq_table[deseq_table["orientation"] == "forward"].drop(columns=["orientation"])
        reverse_deseq = deseq_table[deseq_table["orientation"] == "reverse"].drop(columns=["orientation"])
        if either:
            either_enhancers = self.filter_either_orientation_enhancers(
                    forward_deseq, reverse_deseq, design_dir, logfc_cutoff, padj)
        else:
            either_enhancers = pd.DataFrame()

        both_enhancers = self.filter_both_orientation_enhancers(
                    forward_deseq, reverse_deseq, logfc_cutoff, padj, lfc_method)

        enhancers = pd.concat([both_enhancers, either_enhancers])
        return enhancers

    def filter_either_orientation_enhancers(self, forward_deseq: pd.DataFrame, 
                                            reverse_deseq: pd.DataFrame,
                                            design_dir: str, logfc_cutoff: float, 
                                            padj: float) -> pd.DataFrame:
        """
        Filter either orientation enhancers.
        """
        forward_specific = forward_deseq[~forward_deseq["name"].isin(reverse_deseq["name"])]
        reverse_specific = reverse_deseq[~reverse_deseq["name"].isin(forward_deseq["name"])]

        forward_tested = pd.read_csv(f"{design_dir}/combined_counts_f.txt", sep="\t", header=None)
        reverse_tested = pd.read_csv(f"{design_dir}/combined_counts_r.txt", sep="\t", header=None)
        forward_tested.rename(columns={3: "name"}, inplace=True)
        reverse_tested.rename(columns={3: "name"}, inplace=True)

        forward_specific = forward_specific[
                ~forward_specific["name"].isin(reverse_tested["name"].values)]
        reverse_specific = reverse_specific[
                ~reverse_specific["name"].isin(forward_tested["name"].values)]

        forward_specific["log2FC"] = forward_specific["log2FoldChange"]
        reverse_specific["log2FC"] = reverse_specific["log2FoldChange"]

        forward_specific["enhancer"] = np.where(
                (forward_specific["padj"] < padj) &
                (forward_specific["log2FoldChange"] > logfc_cutoff),
                True, False
                )
        reverse_specific["enhancer"] = np.where(
                (reverse_specific["padj"] < padj) &
                (reverse_specific["log2FoldChange"] > logfc_cutoff),
                True, False
                )

        return pd.concat([forward_specific[["name", "log2FC", "enhancer"]],
                          reverse_specific[["name", "log2FC", "enhancer"]]])

    def filter_both_orientation_enhancers(self, forward_deseq: pd.DataFrame, 
                                          reverse_deseq: pd.DataFrame,
                                          logfc_cutoff: float, padj: float,
                                          lfc_method: str) -> pd.DataFrame:
        """
        Filter both orientation enhancers.
        """
        both_deseq = forward_deseq.merge(reverse_deseq, on=["name", "design"], how="inner", 
                                        suffixes=("_forward", "_reverse"))
        
        calc_method = "mean" if lfc_method == "avg" else "max"
        both_deseq.loc[:, "log2FC"] = both_deseq[
                    ["log2FoldChange_forward", "log2FoldChange_reverse"]
                ].agg(calc_method, axis=1)
        
        both_deseq["enhancer"] = np.where(
                            (both_deseq["padj_forward"] < padj) & 
                            (both_deseq["padj_reverse"] < padj) & 
                            (both_deseq["log2FoldChange_forward"] > logfc_cutoff) & 
                            (both_deseq["log2FoldChange_reverse"] > logfc_cutoff),
                            True, False
                            )
        
        return both_deseq[["name", "log2FC", "enhancer"]]

    def calculate_active_rate(self, design: str, enhancers: pd.DataFrame, 
                              active_rate_file: str):
        """
        Calculate the active rate of enhancers.

        Active rate is defined as: Number of active enhancers / total number 
                                            of mapped elements (passed count cutoff).
        
        Args:
            design (str): Name of the design.
            enhancers (pd.DataFrame): DataFrame of extracted LFC values for elements.
            active_rate_file (str): Path to the output file.
        """
        print(f"Design: {design}")

        number_of_tested = enhancers.shape[0]
        print(f"Number of tested elements: {number_of_tested}")

        number_of_active = enhancers[enhancers["enhancer"]].shape[0]
        print(f"Number of active elements: {number_of_active}")

        active_rate = number_of_active / number_of_tested
        print(f"Active rate: {active_rate:.3f}\n")

        with open(active_rate_file, "w") as f:
            f.write(f"Number of tested elements: {number_of_tested}\n")
            f.write(f"Number of active elements: {number_of_active}\n")
            f.write(f"Active rate: {active_rate:.3f}\n")
        f.close()

    def annotate_extra_columns(self, enhancers: pd.DataFrame, original_bed: str) -> pd.DataFrame:
        """
        Annotate extra columns to the enhancers DataFrame.
        (chrom, start, end)

        Args:
            enhancers (pd.DataFrame): DataFrame of extracted LFC values for elements.
            original_bed (str): Path to the original bed file.
        """
        original_bed_df = pd.read_csv(original_bed, sep="\t", header=None)
        original_bed_df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "name"}, inplace=True)
        original_bed_df = original_bed_df[["name", "chrom", "start", "end"]]
        enhancers = enhancers.merge(original_bed_df, on="name", how="left")
        return enhancers

    def annotate_prominent_tss(self, design: str, enhancers: pd.DataFrame, 
                               bigwigs: str, full_original: str) -> pd.DataFrame:
        """
        Annotate the enhancer dataframe with the prominent TSS.

        Args:
            design (str): Name of the design.
            enhancers (pd.DataFrame): Enhancer regions.
            bigwigs (str): Path to the bigwig files.
            full_original (str): Path to the full original bed file.
        """
        plus_bigwig, minus_bigwig = bigwigs
        prominent_tss = helpers.get_prominent_tss(plus_bigwig, minus_bigwig, design, 
                                                  enhancers, full_original)
        return prominent_tss
    
    def annotate_core_promoters(self, design: str, enhancers: pd.DataFrame,
                                original_bed: str) -> pd.DataFrame:
        """
        Annotate the enhancer dataframe with the core promoters.

        Args:
            design (str): Name of the design.
            enhancers (pd.DataFrame): Enhancer regions.
            original_bed (str): Path to the original bed file.
        """
        svr_dir = f"{args.outdir}/{args.starr}_{args.resolution}/SVR_predictions"
        os.makedirs(svr_dir, exist_ok=True)
        if not os.path.exists(self.fasta):
            raise FileNotFoundError(f"Please have hg38.fa ready at {self.fasta}.")
        
        original_bed_df = pd.read_csv(original_bed, sep="\t", header=None)
        original_bed_df.rename(columns={3: "name", 6: "thickStart", 7: "thickEnd"}, inplace=True)
        original_bed_df = original_bed_df[["name", "thickStart", "thickEnd"]]
        enhancers = enhancers.merge(original_bed_df, on="name", how="left")

        if "tata" in design.lower():
            core_promoters = helpers.add_svr_predictions(design, enhancers, 
                                                         self.svr_tata, svr_dir, self.fasta)
        elif "pause" in design.lower():
            core_promoters = helpers.add_svr_predictions(design, enhancers, 
                                                         self.svr_dpr, svr_dir, self.fasta)
        else:
            enhancers["Groups"] = "Others"
            core_promoters = enhancers
        
        return core_promoters

def main(args):
    """
    Main function to extract enhancers.

    Args:
        args: Namespace object from the argument parser.
    """
    process = ExtractEnhancer()
    process.annotate_enhancers(args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract enhancers.")
    ExtractEnhancer.add_arguments(parser)
    args = parser.parse_args()
    main(args)