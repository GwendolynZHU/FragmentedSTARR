import os
import re
import subprocess

import pybedtools
import pyBigWig
import pandas as pd
import numpy as np

from Bio import SeqIO
pd.set_option('future.no_silent_downcasting', True)

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


def make_counts_dir(outdir: str, file_source: str, resolution: int, design: str, replicate: str) -> str:
    """ 
    Create a directory for storing the processed counts data.

    Args:
        outdir (str): Path to the output directory.
        file_source (str): Source of the STARR-seq data.
        resolution (int): Resolution of the data.
        design (str): Name for the design.
        replicate (str): Replicate number.

    Returns:
        dir_path (str): Path to the newly created directory.
    """
    dir_path = f"{outdir}/{file_source}_{resolution}/{design}/{replicate}"
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
    file = file.iloc[:, :4]
    adjustments = [(i, j) for i in range(-gSize, gSize + 1) for j in range(-gSize, gSize + 1)]
    new_rows = []

    for row in file.itertuples(index=False):
        for i_sta, i_end in adjustments:
            new_rows.append((
                row[0], 
                row[1] + i_sta, 
                row[2] + i_end, 
                row[3],
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
                 4: "overlap"})
    
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
            file = pd.DataFrame({0: ["chrN"], 1: [0], 2: [0], 3: [""], 4: [0]})
        
        file = file.rename(columns={4: f"value_{idx}"})
        if combined_df is None:
            combined_df = file
        else:
            combined_df = pd.merge(
                combined_df, file, on=[0, 1, 2, 3], how="outer"
            )
            combined_df = combined_df.fillna(0)
            combined_df = combined_df[combined_df[0] != "chrN"]

    combined_df["numeric_part"] = combined_df[3].str.extract(r"(\d+)").astype(int)
    combined_df = combined_df.sort_values(by="numeric_part").drop(columns="numeric_part")

    combined_out_path = out_dir + "/combined_counts_" + orientation + ".txt"
    combined_df.to_csv(combined_out_path, sep="\t", header=False, index=False)

    
def count_reads_from_bigwig(bw: str, regions: pd.DataFrame) -> pd.DataFrame:
    """
    Count the reads from a bigwig file.

    Args:
        bw (str): Path to the bigwig file.
        regions (pd.DataFrame): Regions to count the reads from.
    
    Returns:
        pd.DataFrame: A DataFrame containing the read counts.
    """
    bw = pyBigWig.open(bw)
    results = []
    for idx, row in regions.iterrows():
        chrom, start, end = row["chrom"], row["originalStart"], row["originalEnd"]
        if pd.isna(chrom) or pd.isna(start) or pd.isna(end):
            print(chrom, start, end)
        sum_expr = abs(np.sum(np.nan_to_num(bw.values(chrom, int(start), int(end)))))
        results.append({"index": idx, "chrom": chrom, "start": start, \
                        "end": end, "expression": sum_expr})
    bw.close()
    return pd.DataFrame(results)


def get_prominent_tss(pl_bw: str, mn_bw: str, design: str, 
                      data: pd.DataFrame, original_bed: str) -> pd.DataFrame:
    """
    Get the prominent TSS regions.

    Args:
        pl_bw (str): Path to the plus strand bigwig file.
        mn_bw (str): Path to the minus strand bigwig file.
        design (str): Name for the design.
        data (pd.DataFrame): Data containing the regions.
        original_bed (str): Path to the original BED file.
    
    Returns:
        pd.DataFrame: A DataFrame with the prominent TSS regions
            annotated as an extra column.
    """
    full_bed = pd.read_table(original_bed, header=None
                            ).rename(columns={1: "originalStart", 2: "originalEnd", 3: "name"})
    full_bed = full_bed.loc[:, ["originalStart", "originalEnd", "name"]]
    data = data.merge(full_bed, on="name", how="left")
    pl_counts = count_reads_from_bigwig(pl_bw, data)
    mn_counts = count_reads_from_bigwig(mn_bw, data)
    data["pl_counts"] = pl_counts["expression"]
    data["mn_counts"] = mn_counts["expression"]

    if "down" in design:
        data["prominent_tss"] = np.where(data["pl_counts"] > data["mn_counts"], 
                                         "maximum", "minimum")
    elif "up" in design:
        data["prominent_tss"] = np.where(data["pl_counts"] > data["mn_counts"], 
                                         "minimum", "maximum")
    else:
        data["prominent_tss"] = ""

    data = data.drop(columns=["pl_counts", "mn_counts", "originalStart", "originalEnd"])
    return data


def save_temp_seq(design: str, data: pd.DataFrame, svr_dir: str, shifts: list[int]):
    """
    Save the temporary sequence data for SVR input.

    Args:
        design (str): Name for the design.
        data (pd.DataFrame): Data containing the regions.
        svr_dir (str): Path to the SVR predictions directory.
        shifts (list[int]): List of shifts to apply.
    """
    config = {
        ("tata", "down"): (-32, -20, "thickEnd"),
        ("tata", "up"): (20, 32, "thickStart"),
        ("pause", "down"): (17, 36, "thickEnd"),
        ("pause", "up"): (-36, -17, "thickStart"),
    }

    parts = re.split(r"\W+", design.lower())
    motif = next((key[0] for key in config if key[0] in parts), None)
    orientation = next((key[1] for key in config if key[1] in parts), None)

    if motif and orientation:
        start_adj, end_adj, column = config[(motif, orientation)]
    else:
        raise ValueError(f"Could not determine config for input: {design}")

    for shift in shifts:
        df = data.loc[:, ["chrom", column]].copy()
        df["motif_start"] = df[column] + start_adj + shift
        df["motif_end"] = df[column] + end_adj + shift
        df = df[["chrom", "motif_start", "motif_end"]].astype(
                        {"motif_start": int, "motif_end": int})
        file_path = os.path.join(svr_dir, f"tmp_{design}_{shift}.bed")
        df.to_csv(file_path, sep="\t", header=None, index=False)


def reverse_complement(sequence: str) -> str:
    """
    Compute the reverse complement of a DNA sequence.

    Args:
        sequence (str): DNA sequence.
    
    Returns:
        str: Reverse complement of the DNA sequence.
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in reversed(sequence))


def fasta_to_svr_input(fasta: str, svr_input: str, design: str):
    """
    Convert the FASTA file to SVR input format.

    Args:
        fasta (str): Path to the sequence FASTA file.
        svr_input (str): Path to the SVR input file.
        design (str): Name for the design.
    """
    with open(svr_input, "w") as f:
        for record in SeqIO.parse(fasta, "fasta"):
            sequence = reverse_complement(record.seq.upper()) if "up" in design else record.seq.upper()
            f.write(f"{sequence}\n")


def assign_predictions(design: str, prediction_file: str) -> pd.DataFrame:
    """
    Assign the SVR predictions based on different designs.

    Args:
        design (str): Name for the design.
        prediction_file (str): Path to the prediction file.
    
    Returns:
        pd.DataFrame: A DataFrame containing group assignments
            based on the SVR predictions.
    """
    pred_df = pd.read_csv(prediction_file, sep="\t", header=0)

    design_conditions = {
        "tata": {
            "conditions": [
                (pred_df["SVR_score"] >= 1),
                (pred_df["SVR_score"] >= 0.25) & (pred_df["SVR_score"] < 1),
                (pred_df["SVR_score"] < 0.25),
            ],
            "choices": ["StrongTATA", "IntermediateTATA", "NoTATA"],
        },
        "pause": {
            "conditions": [
                (pred_df["SVR_score"] >= 2),
                (pred_df["SVR_score"] >= 1.25) & (pred_df["SVR_score"] < 2),
                (pred_df["SVR_score"] >= 0.5) & (pred_df["SVR_score"] < 1.25),
                (pred_df["SVR_score"] < 0.5),
            ],
            "choices": ["StrongDPR", "IntermediateDPR", "WeakDPR", "NoDPR"],
        }
    }
    config = design_conditions["tata"] if "tata" in design.lower() else design_conditions["pause"]
    pred_df["Groups"] = np.select(config["conditions"], config["choices"])
    return pred_df


def svr_predict(svr: str, design: str, svr_dir: str, shifts: list[int], fasta: str
                ) -> pd.DataFrame:
    """
    Predict the motif presence using SVR.

    Args:
        svr (str): Path to the SVR model.
        design (str): Name for the design.
        svr_dir (str): Path to the SVR predictions directory.
        shifts (list[int]): List of shifts to apply.
        fasta (str): Path to the reference genome FASTA file.
    
    Returns:
        pd.DataFrame: A DataFrame containing the SVR predictions.
    """
    predictions = []

    for shift in shifts:
        tmp_seq_file = os.path.join(svr_dir, f"tmp_{design}_{shift}.bed")
        tmp_fasta_file = os.path.join(svr_dir, f"tmp_{design}_{shift}.fa")
        tmp_svr_file = os.path.join(svr_dir, f"tmp_{design}_{shift}.txt")
        tmp_prediction_file = os.path.join(svr_dir, f"tmp_predictions_{design}_{shift}.txt")

        cmd = f"bedtools getfasta -fi {fasta} -bed {tmp_seq_file} -fo {tmp_fasta_file}"
        subprocess.run(cmd, shell=True)

        fasta_to_svr_input(tmp_fasta_file, tmp_svr_file, design)
        cmd = f"/programs/R-4.2.1-r9/bin/Rscript src/SVRpredict.R -i {tmp_svr_file} -m {svr} -o {tmp_prediction_file}"
        subprocess.run(cmd, shell=True)

        prediction_df = assign_predictions(design, tmp_prediction_file)
        predictions.append(prediction_df)
    
    result_df = pd.DataFrame({
        "SVR_score": predictions[-1]["SVR_score"].values, 
        "Groups": predictions[-1]["Groups"].values
    })

    # Vectorized update of SVR_score and Groups
    for pred_df in predictions[:-1]:
        mask = pred_df["SVR_score"].values > result_df["SVR_score"].values
        result_df.loc[mask, "SVR_score"] = pred_df.loc[mask, "SVR_score"].values
        result_df.loc[mask, "Groups"] = pred_df.loc[mask, "Groups"].values

    # print(result_df)
    os.system(f"rm -r {svr_dir}/")
    return result_df



def add_svr_predictions(design: str, data: pd.DataFrame, 
                        svr: str, svr_dir: str, fasta: str) -> pd.DataFrame:
    """
    Add the SVR predictions to the data.

    Args:
        design (str): Name for the design.
        data (pd.DataFrame): Data containing the regions.
        svr (str): Path to the SVR models.
        svr_dir (str): Path to the SVR predictions directory.
        fasta (str): Path to the reference genome FASTA file.
    
    Returns:
        pd.DataFrame: A DataFrame with the SVR predictions
            annotated as an extra column.
    """
    ### input to svr
    ### predictions in a column
    shifts = range(-2,3) # allow for 2bp shifts
    save_temp_seq(design, data, svr_dir, shifts)
    svr_annotation = svr_predict(svr, design, svr_dir, shifts, fasta)
    data["Groups"] = svr_annotation["Groups"]
    return data
