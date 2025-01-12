import os
import subprocess

import pybedtools
import pandas as pd

from concurrent.futures import ThreadPoolExecutor


def process_chunk(counts_path: str, start: int, num_lines: int
                  ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ 
    Process a chunk of the raw counts file.

    Args:
        counts_path (str): Path to a specific replicate STARR raw counts file.
        start (int): Starting line number of the chunk.
        num_lines (int): Number of lines to read.
    
    Returns:
        forward, reverse (pd.DataFrame, pd.DataFrame): 
            Forward and reverse orientation separated raw chunk data.
    """
    chunk = pd.read_csv(counts_path, sep="\t", header=None, 
                    skiprows=start, nrows=num_lines,
                    compression="gzip" if counts_path.endswith(".gz") else None)
    
    forward = chunk[chunk[3] == "+"]
    reverse = chunk[chunk[3] == "-"]
    return forward, reverse


def get_orientation_specific_raw_counts(counts_path: str, num_chunks: int
                                        ) -> tuple[pd.DataFrame, pd.DataFrame]:
    """ 
    Return the orientation-separated STARR counts file.

    Args:
        counts_path (str): Path to a specific replicate STARR raw counts file.
        num_chunks (int): Number of cores to use.
    
    Returns:
        forward_df, reverse_df (pd.DataFrame, pd.DataFrame): 
            Forward and reverse orientation separated raw counts data.
    """
    lc_command = "gunzip -c " + counts_path + " | wc -l" if counts_path.endswith(".gz") else "wc -l < " + counts_path
    total_lines = int(subprocess.run(lc_command, 
                                        shell=True, 
                                        stdout=subprocess.PIPE).stdout)

    chunk_size = total_lines // num_chunks
    if total_lines % num_chunks != 0:
        chunk_size += 1
    
    # Read the counts file parallelly
    tasks = [(counts_path, i * chunk_size, min(chunk_size, 
                                        total_lines - i * chunk_size)) for i in range(num_chunks)]
    
    # Process the chunks in parallel
    with ThreadPoolExecutor(max_workers=44) as executor:
        results = list(executor.map(lambda args: process_chunk(*args), tasks))
    
    # Combine all chunks into a single DataFrame
    forward_chunks = [res[0] for res in results]
    reverse_chunks = [res[1] for res in results]
    forward_df = pd.concat(forward_chunks, ignore_index=True)
    reverse_df = pd.concat(reverse_chunks, ignore_index=True)
    return forward_df, reverse_df


def make_counts_dir(outdir: str, file_source: str, design: str, replicate: str) -> str:
    """ 
    Create a directory for storing the processed counts data.

    Args:
        outdir (str): Path to the output directory.
        file_source (str): Source of the STARR-seq data.
        design (str): Name for the design.
        replicate (str): Replicate number.

    Returns:
        dir_path (str): Path to the newly created directory.
    """
    dir_path = f"{outdir}/{file_source}/{design}/{replicate}"
    os.makedirs(dir_path, exist_ok=True)
    return dir_path


def generate_region_fragments(file: pd.DataFrame, gSize= 5) -> pd.DataFrame:
    """
    Generate a sorted dataframe of targeted fragments.
    
    Default resolution is 5bp at the edges: This means that the regions 
    with at most 5bp differences (on either ends) with the original regions 
    are all counted as possible matches.

    Args: 
        file (pd.DataFrame): Input file containing the regions.
        gSize (int): The resolution of the regions. Default is 5bp.

    Returns:
        pd.DataFrame: A DataFrame containing all possible regions.
    """
    print("Generating regions...\n")
    adjustments = [(i, j) for i in range(-gSize, gSize + 1) for j in range(-gSize, gSize + 1)]
    new_rows = []

    for row in file.itertuples(index=False):
        for i_sta, i_end in adjustments:
            new_rows.append((
                row[0], 
                row[1] + i_sta, 
                row[2] + i_end, 
                row[3], row[4], row[5], row[6], row[7]
            ))
    
    new_df = pd.DataFrame(new_rows, columns=file.columns)
    return new_df.sort_values(by=["chrom", "start", "end"])


def align(starr_counts: pd.DataFrame, region_ref: pd.DataFrame) -> pd.DataFrame:
    """
    Align the sorted STARR fragment counts with target regions.
    A reciprocal overlap of 100% is required for a match.

    Args:
        starr_counts (pd.DataFrame): STARR fragment counts.
        region_ref (pd.DataFrame): Target regions.
    
    Returns:
        pd.DataFrame: A DataFrame containing the aligned STARR fragment counts.
        a list that contains 100% sequence coverage file:
    (chr, start, end, counts) #TODO: check if this is the correct format
    """
    starr_counts = pybedtools.BedTool.from_dataframe(starr_counts)
    region_ref = pybedtools.BedTool.from_dataframe(region_ref)

    overlap = region_ref.coverage(starr_counts, sorted=True, counts=True, f=1.0, r=True, n=44)
    overlap = overlap.to_dataframe(disable_auto_names=True, header=None)
    overlap = overlap.rename(
        columns={0: "chrom", 
                 1: "start",
                 2: "end",
                 3: "name",
                 4: "score",
                 5: "strand",
                 6: "thickStart",
                 7: "thickEnd",
                 8: "overlap"})
    
    overlap = overlap[overlap["overlap"]==1].reset_index(drop=True)
    return overlap


def align_wrapper(args):
    """
    Wrapper function to align the STARR counts with target regions.

    Args:
        args (tuple): A tuple containing the STARR counts and target regions.
    
    Returns:
        pd.DataFrame: A DataFrame containing the aligned STARR counts.
    """
    return align(*args)


def count_mapped_fragments(aligned_counts: pd.DataFrame, raw_counts_df: pd.DataFrame,
                           original_ref_df: pd.DataFrame, umi: bool, 
                           out_dir: str, orientation: str):
        """
        Count the mapped fragments and save the data.


        Args:
            aligned_counts (pd.DataFrame): Aligned fragment counts.
            raw_counts_df (pd.DataFrame): Raw STARR counts.
            original_ref_df (pd.DataFrame): Original reference regions.
            umi (bool): Whether the original dataset deduplicated using UMI.
            out_dir (str): Path to the output directory.
        """
        raw_counts_df = raw_counts_df.rename(
                columns={0: "chrom", 
                         1: "start", 
                         2: "end", 
                         3: "strand", 
                         4: "counts"})

        if umi: # already UMI deduplicated - sum the counts in raw data
            aligned_counts["count"] = raw_counts_df.merge(aligned_counts, 
                                                on=["chrom", "start", "end"], how="inner")["counts"]
            aligned_counts = aligned_counts.loc[:, ["name", "count"]]
        else:
            aligned_counts = aligned_counts.loc[:, ["name", "overlap"]]
            aligned_counts = aligned_counts.rename(columns={"overlap": "count"})

        aligned_counts = aligned_counts.groupby("name", as_index=True).agg(
                                                {"count": "sum"})

        # Reassign chr, start, and end given the name index
        if not aligned_counts.empty:
            annotated_aligned_counts = aligned_counts.merge(original_ref_df, 
                                                left_index=True, right_on="name", how="left").reset_index()
            annotated_aligned_counts = annotated_aligned_counts.loc[:, ["chrom", "start", "end", "name", "count"]]
        else:
            annotated_aligned_counts = pd.DataFrame()

        # Save the counts data
        out_path = out_dir + "/aligned_count_" + orientation + ".bed"
        annotated_aligned_counts.to_csv(out_path, sep="\t", header=False, index=False)


def combine_replicates(file_list: list[str], out_dir: str, orientation: str, ):
    """
    Combine the replicate data into a single file.

    Args:
        file_list (list[str]): List of paths to the replicate files.
        out_dir (str): Path to the output directory.
        orientation (str): Orientation of the data.
    """
    combined_df = None

    for idx, filepath in enumerate(file_list):
        file = pybedtools.BedTool(filepath).to_dataframe(disable_auto_names=True, header=None)
        if file.empty:
            print(f"Skipping empty file: {filepath}")
            continue
        
        file = file.rename(columns={4: f"value_{idx}"})
        if combined_df is None:
            combined_df = file
        else:
            combined_df = pd.merge(
                combined_df, file, on=[0, 1, 2, 3], how="outer"
            )
            combined_df = combined_df.fillna(0)

    combined_df["numeric_part"] = combined_df[3].str.extract(r"(\d+)").astype(int)
    combined_df = combined_df.sort_values(by="numeric_part").drop(columns="numeric_part")

    combined_out_path = out_dir + "/combined_counts_" + orientation + ".bed"
    combined_df.to_csv(combined_out_path, sep="\t", header=False, index=False)

    