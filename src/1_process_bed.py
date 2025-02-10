#!/usr/bin/env python

import os
import argparse
import pybedtools
import tempfile

import helpers

import pandas as pd

from concurrent.futures import ProcessPoolExecutor


class ProcessSTARR:
    """
    Generate counts table from STARR raw data given a reference BED file.
    """
    def __init__(self, input_path: str, file_source: str, 
                       design: str, outdir= "./output/", 
                       resolution= 5, ncores= 8):
        """
        Initialize ProcessSTARR process.

        Args:
            input_path (str): Input BED file path for reads extraction. 
                    PINTS format (chr, start, end, name, score, strand, tss_start, tss_end). #TODO: Check this.
            file_source (str): Raw STARR-seq data source (WHG-STARR, deep-ATAC-STARR).
            design (str): Name for the design (name for the folder access).
            outdir (str): Output directory path. Default is "./output/".
            resolution (int): Element resolution. Default is 5.
            ncores (int): Number of cores to use. Default is 8.
            
        """
        self._temp_dir = tempfile.TemporaryDirectory()
        self.input_path = input_path
        self.file_source = file_source
        self.design = design
        self.outdir = outdir
        self.resolution = resolution
        self.ncores = ncores

    @staticmethod
    def add_data_specific_args(parser):
        """
        Add data-specific arguments to the parser.

        Args:
            parser: ArgumentParser object
        """
        parser.add_argument("--input_path", type=str, required=True, help="Input BED file path")
        parser.add_argument("--file_source", choices=["WHG-STARR", "deep-ATAC-STARR"], required=True, help="File source")
        parser.add_argument("--design", type=str, required=True, help="Name for the design")
        parser.add_argument("--outdir", type=str, default="./output/", help="Output directory")
        parser.add_argument("--resolution", type=int, default=5, help="Resolution of the elements")
        parser.add_argument("--n_cores", type=int, default=8, help="Number of cores to use")

    def extract_reads(self, args):
        """
        Extract reads from different STARR-seq datasets. 
        This method can be overridden by subclasses.
        """
        raise NotImplementedError("Subclasses must implement this method.")
    
    def process_replicate(self, raw_path: str, elements_df: pd.DataFrame,
                          input_file_df: pd.DataFrame, umi: bool, dir_path: str):
        """
        Process a single replicate.

        Args:
            raw_path (str): Path to raw STARR counts file.
            elements_df (pd.DataFrame): Elements dataframe.
            input_file_df (pd.DataFrame): Input BED region dataframe.
            umi (bool): Whether the raw data deduplicated by UMI.
            dir_path (str): Path to the directory to store the processed counts.
        """
        raw_forward_df, raw_reverse_df = helpers.get_orientation_specific_raw_counts(raw_path, self.ncores)

        tasks = [(raw_forward_df, elements_df), (raw_reverse_df, elements_df)]

        with ProcessPoolExecutor(max_workers=2) as executor:
            results = list(executor.map(helpers.align_wrapper, tasks))
            aligned_forward, aligned_reverse = results

        helpers.count_mapped_fragments(aligned_forward, raw_forward_df, input_file_df, umi, dir_path, "f")
        helpers.count_mapped_fragments(aligned_reverse, raw_reverse_df, input_file_df, umi, dir_path, "r")
    
    def merge_replicates(self, dnas, rnas):
        """
        Generate the orientation-independent counts table across DNA and RNA replicates.
        """
        combined_path = f"{self.outdir}/{self.file_source}_{self.resolution}/{self.design}"

        for orientation in ["f", "r"]:
            DNA_files = [f"{combined_path}/DNA{dna_idx}/aligned_count_{orientation}.bed" 
                                for dna_idx in range(1, dnas+1)]
            RNA_files = [f"{combined_path}/RNA{rna_idx}/aligned_count_{orientation}.bed" 
                                for rna_idx in range(1, rnas+1)]

            ls = DNA_files + RNA_files
            helpers.combine_replicates(ls, combined_path, orientation)
    
    def configure_tempdir(self):
        """
        Configure temporary directory for processing.
        """
        pybedtools.helpers.set_tempdir(self._temp_dir.name)

    def cleanup(self):
        """
        Clean up temporary files.
        """
        pybedtools.cleanup(remove_all=True)
        self._temp_dir.cleanup()


class ProcessSTARRDeepATAC(ProcessSTARR):
    """
    Generate counts table from STARR raw data given a reference BED file.

    Attributes:
        DNA_path (str): path to deep-ATAC-STARR DNA raw counts data
        RNA_path (str): path to deep-ATAC-STARR RNA raw counts data
        new_RNA_path (str): path to updated experiment RNA raw counts data
        dnas (int): number of deep-ATAC-STARR DNA replicates
        rnas (int), number of deep-ATAC-STARR RNA replicates
        dna_umi (bool): Whether original raw DNA data deduplicates reads
        rna_umi (bool): Whether original raw RNA data deduplicates reads
    """
    def __init__(self, args):
        super().__init__(args.input_path, args.file_source, args.design, args.outdir, args.resolution, args.n_cores)
        self.DNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/deep_ATAC_STARR/processing_data_v1/out_DNA_no_UMI/"
        self.RNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/deep_ATAC_STARR/processing_data_v1/out_RNA_with_UMI/"
        self.dnas = 6
        self.rnas = 4
        self.dna_umi = False
        self.rna_umi = True

    def extract_reads(self, args):
        """
        Extract reads from deep-ATAC-STARR dataset. (ref: link TODO)

        Args:
            input_path (str): Input BED file path
            file_source (str): Raw STARR-seq data source (WHG-STARR, deep-ATAC-STARR)
            design (str): Name for the design
            outdir (str): Output directory path
            resolution (int): Element resolution, default is 5
        """
        input_file_df = pybedtools.BedTool(args.input_path).to_dataframe(disable_auto_names=False, header=None)
        elements_df = helpers.generate_region_fragments(input_file_df, args.resolution)

        for dna_idx in range(1, self.dnas + 1):
            print("Mapping to DNA replicate " + str(dna_idx) + "... ")
            dir_path = helpers.make_counts_dir(args.outdir, args.file_source, args.resolution, 
                                               args.design, f"DNA{dna_idx}")

            dna_raw_count_path = self.DNA_path + f"DNA{dna_idx}/all/count.bed.gz"
            super().process_replicate(dna_raw_count_path, elements_df, input_file_df, self.dna_umi, dir_path)
            
        for rna_idx in range(1, self.rnas + 1):
            print("Mapping to RNA replicate " + str(rna_idx) + "... ")
            dir_path = helpers.make_counts_dir(args.outdir, args.file_source, args.resolution,
                                               args.design, f"RNA{rna_idx}")

            rna_paths = {
                1: "RNA1/all/count.bed.gz",
                2: "corrected_bam/RNA1/all/count.bed",
                3: "RNA3/all/count.bed.gz",
                4: "corrected_bam_RNA4/RNA1/all/count.bed"
            }
            rna_raw_count_path = self.RNA_path + rna_paths.get(rna_idx)
            super().process_replicate(rna_raw_count_path, elements_df, input_file_df, self.rna_umi, dir_path)


class ProcessSTARRWHG(ProcessSTARR):
    """
    Generate counts table from STARR raw data given a reference BED file.

    Attributes:
        DNA_path (str): path to WHG-STARR DNA raw counts data
        RNA_path (str): path to WHG-STARR RNA raw counts data
        dnas (int): number of WHG-STARR DNA replicates
        rnas (int), number of WHG-STARR RNA replicates
        dna_umi (bool): Whether original raw DNA data deduplicates reads
        rna_umi (bool): Whether original raw RNA data deduplicates reads
    """
    def __init__(self, args):
        super().__init__(args.input_path, args.file_source, args.design, args.outdir, args.resolution, args.n_cores)
        self.DNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/WHG_STARR_seq_TR/processing_data_v1/out_DNA_no_UMI/"
        self.RNA_path = "/fs/cbsuhy01/storage/jz855/STARR_seq_dataset/WHG_STARR_seq_TR/processing_data_v1/out_RNA_no_UMI/"
        self.dnas = 1
        self.rnas = 3
        self.dna_umi = False
        self.rna_umi = False

    def extract_reads(self, args):
        """
        Extract reads from WHG-STARR dataset. (ref: link TODO)

        Args:
            input_path (str): Input BED file path
            file_source (str): Raw STARR-seq data source (WHG-STARR, deep-ATAC-STARR)
            design (str): Name for the design
            outdir (str): Output directory path
            resolution (int): Element resolution, default is 5
        """
        input_file_df = pybedtools.BedTool(args.input_path).to_dataframe(disable_auto_names=False, header=None)
        elements_df = helpers.generate_region_fragments(input_file_df, args.resolution)

        for dna_idx in range(1, self.dnas + 1):
            print("Mapping to DNA replicate " + str(dna_idx) + "... ")
            dir_path = helpers.make_counts_dir(args.outdir, args.file_source, args.resolution,
                                               args.design, f"DNA{dna_idx}")

            dna_raw_count_path = self.DNA_path + f"DNA{dna_idx}/all/count.bed"
            super().process_replicate(dna_raw_count_path, elements_df, input_file_df, self.dna_umi, dir_path)

        for rna_idx in range(1, self.rnas + 1):
            print("Mapping to RNA replicate " + str(rna_idx) + "... ")
            dir_path = helpers.make_counts_dir(args.outdir, args.file_source, args.resolution,
                                               args.design, f"RNA{rna_idx}")

            rna_raw_count_path = self.RNA_path + f"RNA{rna_idx}/all/count.bed"
            super().process_replicate(rna_raw_count_path, elements_df, input_file_df, self.rna_umi, dir_path)

    
def main(args):
    """
    Main function to process STARR-seq data.

    Args:
        args: Namespace object
    """
    os.makedirs(args.outdir, exist_ok=True)

    if args.file_source == "deep-ATAC-STARR":
        process = ProcessSTARRDeepATAC(args)
    elif args.file_source == "WHG-STARR":
        process = ProcessSTARRWHG(args)
    else:
        raise ValueError("Invalid file source")

    process.configure_tempdir()

    process.extract_reads(args)
    process.merge_replicates(process.dnas, process.rnas)
    process.cleanup()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process STARR-seq data.")
    ProcessSTARR.add_data_specific_args(parser)
    args = parser.parse_args()
    main(args)