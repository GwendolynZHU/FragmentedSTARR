"""
Downstream analysis.

Sample run:
python3 3_filter_analysis.py -o ../new_data --enh_file ../GRO_cap_PINTS_qpeak_call_117/K562_GROcap_hg38_1.1.7_qpeak_calls_1_divergent_peaks_element_60bp.bed --design TSS pause_site --either True
Author: Yutong Zhu
Date: 2024-4-24
"""
import os
import pandas as pd
import numpy as np
import argparse
import pybedtools
from pybedtools import BedTool
import tempfile
import shutil
from helpers import judge_grocap_groups, pybedtools_read_without_header, process_bed_file
from Bio import SeqIO

temp_dir = tempfile.mkdtemp()
pybedtools.helpers.set_tempdir(temp_dir)


def _save_lfc_avg_partial(deseq_design, design, data_dir, starr, outdir):
    """
    Helper function, Save the annotated deseq results for the design

    Parameters:
    deseq_design: Dataframe, subset of the original deseq results
    design: String, one of ("TSS", "pause_site", "INR") - input from users
    data_dir: String, folder where all the data are saved
    starr: String, deep_ATAC_STARR or WHG_STARR
    outdir: String, folder where all DESeq results are saved
    """
    deseq_ls = []
    for orientation in ("f", "r"):
        deseq_design_ori = deseq_design[deseq_design.index.to_series().apply(lambda x: x.split("_")[-1][0]) == orientation]
        deseq_design_ori["index"] = deseq_design_ori.index.astype(str)
        deseq_design_ori[["chrom", "designStart", "designEnd", "orientation", "name"]] = "nan"
        for strand in ("p", "n", "b"):
            name = design+"_"+strand+"_"+orientation
            design_ref_path = os.path.join(data_dir, "design_ref", "divergent_60bp_without_"+design+"_"+strand+".bed")
            design_ref = pybedtools_read_without_header(design_ref_path)[[0,1,2,3,6,7]]
            
            partial_file_path = os.path.join(data_dir, starr, design+"_"+strand, "srt_"+name+".bed")
            partial = process_bed_file(partial_file_path, name)
            partial = partial.merge(design_ref, on=[0,1,2], how="left")
            partial.rename(columns={0:"chrom", 1:"designStart", 2:"designEnd", 3:"name", 6:"thickStart", 7:"thickEnd"}, inplace=True)
            partial["orientation"] = strand
            deseq_design_ori.loc[
                deseq_design_ori["index"].isin(partial["index"]), 
                ["chrom", "designStart", "designEnd", "orientation", "name", "thickStart", "thickEnd"]
                ] = partial.loc[partial["index"].isin(deseq_design_ori["index"]), 
                                ["chrom", "designStart", "designEnd", "orientation", "name", "thickStart", "thickEnd"]].values
        # deseq_design_ori.drop(columns=["index"], axis=1, inplace=True)
        deseq_ls.append(deseq_design_ori[["log2FoldChange", "chrom", "designStart", "designEnd", "orientation", "name", "thickStart", "thickEnd"]])
        deseq_design_ori.to_csv(os.path.join(outdir, "results_ori", "DE_results_"+design+"_"+orientation+".txt"), sep="\t", header=True, index=False)

    ### save a combined averaged LFC file
    deseq_avg = deseq_ls[0].merge(deseq_ls[1], on=["chrom", "designStart", "designEnd", "orientation", "name", "thickStart", "thickEnd"], how="outer")
    deseq_avg["log2FoldChange"] = deseq_avg.apply(lambda row: np.nanmean([row["log2FoldChange_x"], row["log2FoldChange_y"]]), axis=1)
    deseq_avg = deseq_avg[["log2FoldChange", "chrom", "designStart", "designEnd", "orientation", "name", "thickStart", "thickEnd"]]
    deseq_avg.to_csv(os.path.join(outdir, "DE_results_"+design+".txt"), sep="\t", header=True, index=False)


def _save_lfc_avg_full(outdir):
    """ 
    Helper function, return a list with the annotated deseq results for the design

    Parameters:
    outdir: String, folder where all DESeq files are saved
    """
    deseq_full = pd.read_table(os.path.join(outdir, "DE_results_full.txt"))
    deseq_full["index"] = deseq_full.index.values
    enh_ref = pybedtools_read_without_header(args.enh_file)[[0,1,2,3,6,7]]

    deseq_ls = []
    for orientation in ("f", "r"):
        all_full = process_bed_file(os.path.join(args.outdir, args.starr, "full", "srt_full_"+orientation+".bed"), "full_"+orientation)
        deseq_full_orient = deseq_full.merge(all_full, on="index", how="inner")
        deseq_full_orient["full_index"] = deseq_full_orient["index"].apply(lambda x: x.split("_")[-1][1:]).values
        deseq_full_orient = deseq_full_orient.merge(enh_ref, on=[0,1,2], how="left")
        deseq_full_orient.rename(columns={0:"chrom", 1:"designStart", 2:"designEnd",3:"name", 6:"thickStart", 7:"thickEnd"}, inplace=True)
        deseq_ls.append(deseq_full_orient[["log2FoldChange", "chrom", "designStart", "designEnd", "full_index", "name", "thickStart", "thickEnd"]])

    deseq_avg = deseq_ls[0].merge(deseq_ls[1], on=["full_index", "chrom", "designStart", "designEnd", "name", "thickStart", "thickEnd"], how="inner")
    deseq_avg["log2FoldChange"] = deseq_avg.apply(lambda row: np.nanmean([row["log2FoldChange_x"], row["log2FoldChange_y"]]), axis=1)
    deseq_avg = deseq_avg[["log2FoldChange", "chrom", "designStart", "designEnd", "name", "thickStart", "thickEnd"]]

    if args.either:
        either_full_f = process_bed_file(os.path.join(args.outdir, args.starr, "full", "unpaired", "srt_full_either_f.bed"), "either_f")
        either_full_r = process_bed_file(os.path.join(args.outdir, args.starr, "full", "unpaired", "srt_full_either_r.bed"), "either_r")
        either_full = pd.concat([either_full_f, either_full_r])
        
    deseq_either = deseq_full.merge(either_full, on=["index"], how="inner")
    deseq_either = deseq_either.merge(enh_ref, on=[0,1,2], how="left")
    deseq_either.rename(columns={0:"chrom", 1:"designStart", 2:"designEnd",3:"name",6:"thickStart", 7:"thickEnd"}, inplace=True)
    deseq_either = deseq_either[["log2FoldChange", "chrom", "designStart", "designEnd", "name", "thickStart", "thickEnd"]]

    deseq_output = pd.concat([deseq_avg, deseq_either])
    deseq_output.to_csv(os.path.join(outdir, "DE_results_full_annotated.txt"), sep="\t", header=True, index=False)


def save_lfc_avg(deseq, design, data_dir, starr, outdir):
    """ 
    Save the averaged (lfc averaged, considering elements inserted with different orientations) deseq file.

    Parameter:
    deseq: DataFrame, output with all designs and negative control from DESeq2
    design: String, one of ("TSS", "pause_site", "INR") - input from users
    data_dir: String, folder where all the data are saved
    starr: String, deep_ATAC_STARR or WHG_STARR
    outdir: String, output file folder
    """
    ori_dir = os.path.join(outdir, "results_ori")
    os.makedirs(ori_dir, exist_ok=True)
    deseq_design = deseq[deseq.index.to_series().apply(lambda x: x.split("_")[0]) == design.split("_")[0]]

    if design != "full":
        _save_lfc_avg_partial(deseq_design, design, data_dir, starr, outdir)
    else:
        _save_lfc_avg_full(outdir)


def save_tmp_seq(input, design, orientation, outdir):
    """ 
    """
    config = {
        ("TSS", "p"): (-32, -20, "thickEnd"),
        ("TSS", "n"): (20, 32, "thickStart"),
        ("pause_site", "p"): (17, 36, "thickEnd"),
        ("pause_site", "n"): (-36, -17, "thickStart"),
    }
    start, end, column = config[(design, orientation)]
    df = input.copy()
    df["seq_start"] = df[column] + start
    df["seq_end"] = df[column] + end
    df = df.loc[:,["chrom", "seq_start", "seq_end"]]
    df = df.astype({"seq_start": int, "seq_end": int})
    df.to_csv(os.path.join(outdir, "tmp_"+design+"_"+orientation+".bed"), sep="\t", header=None, index=None)


def reverse_complement(dna_seq):
    """ 
    """
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    reverse_complement_seq = "".join(complement[base] for base in reversed(dna_seq))
    return reverse_complement_seq


def fasta_to_txt(orientation, fasta, output_file):
    """ 
    """
    with open (output_file, "w") as out_handle:
        for record in SeqIO.parse(fasta, "fasta"):
            sequence = record.seq.upper()
            if orientation == "r":
                seq_reversed = reverse_complement(sequence)
                out_handle.write(f"{seq_reversed}\n")
            else:
                out_handle.write(f"{sequence}\n")
    return output_file


def assign_groups(design, prediction_file):
    """ 
    """
    pred_df = pd.read_table(prediction_file)
    if design == "TSS":
        conditions = [(pred_df["SVR_score"]>=1), \
                      (pred_df["SVR_score"]>=0.25)&(pred_df["SVR_score"]<1), (pred_df["SVR_score"]<0.25)]
        choices = ["StrongTATA", "IntermediateTATA", "NoTATA"]
    else:
        conditions = [(pred_df["SVR_score"]>=2), \
                      (pred_df["SVR_score"]>=0.75)&(pred_df["SVR_score"]<2), \
                      (pred_df["SVR_score"]>=0.3)&(pred_df["SVR_score"]<0.75), (pred_df["SVR_score"]<0.3)]
        choices = ["GoodDPR", "IntermediateDPR", "WeakDPR", "NoDPR"]
        
    pred_df["Groups"] = np.select(conditions,choices)
    pred_df.to_csv(prediction_file, sep="\t", header=True, index=None)
    return pred_df


def add_svr_predictions(deseq_df, design, outdir):
    """ 
    """
    svr_dir = os.path.join(outdir, "fetch_seq")
    os.makedirs(svr_dir, exist_ok=True)

    fasta_ref = "/fs/cbsuhy01/storage/yz2676/ref/hg38.fa"
    svr_model = "SVRtata_TATAbox.model" if design == "TSS" else "SVRb_DPR.model"

    for orientation in ("p", "n"):
        df_subselect = deseq_df["orientation"] == orientation
        deseq = deseq_df[df_subselect]
        
        # get the sequence input for svr model
        save_tmp_seq(deseq, design, orientation, svr_dir)
        tmp_seq_file = os.path.join(svr_dir, "tmp_"+design+"_"+orientation+".bed")
        fasta = os.path.join(svr_dir, design+"_"+orientation+".fasta")
        os.system("bedtools getfasta -fi "+fasta_ref+" -bed "+tmp_seq_file+" -fo "+fasta)

        seq_outfile = os.path.join(svr_dir, design+"_"+orientation+".txt")
        seq_file = fasta_to_txt(orientation, fasta, seq_outfile)

        # model prediction
        prediction_file = os.path.join(svr_dir, "prediction_"+design+"_"+orientation+".txt")
        os.system("/programs/R-4.2.1-r9/bin/Rscript SVRpredict.R -i "+seq_file+" -m "+svr_model+" -o "+prediction_file)
        prediction_df = assign_groups(design, prediction_file)

        # update final df with svr_predictions
        if not "Groups" in deseq_df.columns:
            deseq_df["Groups"] = "nan"
        deseq_df.loc[df_subselect, "Groups"] = prediction_df["Groups"].values
        deseq_df.fillna("nan")

    # remove tmp files
    os.system("rm "+svr_dir+"/"+design+"*")
    os.system("rm "+svr_dir+"/tmp*")
    
    return deseq_df
    

def get_visualization_script(deseq_file, outdir, starr, ori_dir):
    """ 
    """
    visual_dir = os.path.join(outdir, "visualization", starr, ori_dir)
    os.makedirs(visual_dir, exist_ok=True)

    full_enh = deseq_file[deseq_file.loc[:,"orientation"] == "nan"][["chrom", "start", "end"]]
    full_enh = full_enh.astype({"start": int, "end": int})

    # print(full_enh)
    full_enh.to_csv(os.path.join(visual_dir, "tmp_full_visual.bed"), sep="\t", header=None, index=None)
    cmds = ["bedtools igv -i "+visual_dir+"/tmp_full_visual.bed -path "+visual_dir+" > "+visual_dir+"/test.bs",
            "rm "+visual_dir+"/tmp_full_visual.bed"]
    for cmd in cmds:
        os.system(cmd)


def parse_args():
    parser = argparse.ArgumentParser(description='Downstream filtering of elements')
    parser.add_argument("-o", '--outdir', required=True, help="Output directory")
    parser.add_argument('--enh_file', required=True, help="Path to original .bed enhancer file, should be in the format of eight or four required fields. (chr, chromStart, chromEnd, name, score, strand, thickStart, thickEnd).")
    parser.add_argument('--design', nargs="+", required=True, help="The names of different designs. e.g. TSS, pause_site")
    parser.add_argument('--starr', default="deep_ATAC_STARR", help="One of the STARR datasets (deep_ATAC_STARR, WHG_STARR)")
    parser.add_argument('--either', default=False, help="Whether to include either-orientation enhancers")
    parser.add_argument('--orientation_independent', default=True, help="Whether to use avg LFC and combine element inserted forwardly/reversely")

    return parser.parse_args()


def main(args):
    """ 
    Filtering.
    Grouping.
    Filter out proximal enhancers, enhancers overlapping with other enhancers, elements with GROcap-reads less than 5
    """
    ### Global variables
    ori_dir = "either_ori" if args.either else "ori_ind"
    outdir = os.path.join(args.outdir, "DESeq", args.starr, ori_dir)
    bw1_file_path = "/fs/cbsuhy01/storage/jz855/Reference/K562/GRO_cap/K562_GROcap_hg38_aligned_pl.bw"
    bw2_file_path = "/fs/cbsuhy01/storage/jz855/Reference/K562/GRO_cap/K562_GROcap_hg38_aligned_mn.bw"

    ref = BedTool(args.enh_file) ## info: chr, +60bp srt & end, TSS srt & end
    # tss = BedTool("/fs/cbsuhy01/storage/yz2676/ref/gencode/gencode.v37.annotation.1kb.TSS.sorted.bed.gz")
    junke_tss = BedTool("/fs/cbsuhy01/storage/jz855/Reference/hg38/proximal_v45/promoter_1kbp_protein_coding_TSS_centered_gencode_v45.bed")

    ### save the avg LFC files for each design
    deseq_file = pd.read_csv(os.path.join(outdir, "DE_results_nctrl.txt"), sep="\t", index_col=0)
    for design in args.design:
        save_lfc_avg(deseq_file, design, args.outdir, args.starr, outdir)
    save_lfc_avg(deseq_file, "full", args.outdir, args.starr, outdir)
    full_enh_file = pd.read_csv(os.path.join(outdir, "DE_results_full_annotated.txt"), sep="\t")

    ### step 1: filter out proximal enhancers and save only those partial that have a corresponding full element
    full = BedTool.from_dataframe(pybedtools_read_without_header(os.path.join(args.outdir, args.starr, "full", "srt_full_"+ori_dir+"_e.bed"))[[0,1,2]])
    overlap = full.coverage(junke_tss, counts=True, f=0.9)
    
    ### step 2: filter out enhancers overlapping with other enhancers
    non_overlap = overlap.intersect(ref, C=True).to_dataframe(disable_auto_names=True, header=None)
    non_overlap = non_overlap[(non_overlap[4] == 1) & (non_overlap[3] == 0)].drop(columns=[3,4])
    non_overlap.rename(columns={0:"chrom", 1:"designStart", 2:"designEnd"}, inplace="True")

    ### whole element:
    non_overlap = non_overlap.merge(full_enh_file, on=["chrom", "designStart", "designEnd"], how="left")
    non_overlap.rename(columns={"designStart": "start", "designEnd": "end"}, inplace=True)

    ref = ref.to_dataframe()
    ### select those partial elements that have corresponding full elements
    for design in args.design:
        deseq_partial = pd.read_csv(os.path.join(outdir, "DE_results_"+design+".txt"), sep="\t")
        deseq_merged_partial = deseq_partial.merge(non_overlap.drop(columns="log2FoldChange"), on=["name", "chrom", "thickStart", "thickEnd"], how="inner")
        ### check GROcap signal within the peak and define maximum TSS and minimum TSS
        deseq_merged_partial = judge_grocap_groups(deseq_merged_partial, bw1_file_path, bw2_file_path)
        deseq_merged_full = non_overlap[non_overlap["name"].isin(deseq_merged_partial["name"])]
        deseq_merged_full[["designStart", "designEnd"]] = deseq_merged_full[["start", "end"]].values
        deseq_merged_full[["orientation", "GROcap_signal"]] = "nan"
        
        deseq_output = pd.concat([deseq_merged_full, deseq_merged_partial])

        ## get the TATA, DPR predictions
        if design in ("TSS", "pause_site"):
            deseq_output_annotated = add_svr_predictions(deseq_output, design, outdir)
        else:
            deseq_output_annotated = deseq_output
        deseq_output_annotated.to_csv(os.path.join(outdir, "DE_results_"+design+"_annotated.txt"), sep="\t", header=True, index=None)
    print("Finished filtering the pairwise elements.")

    visual_deseq_ls = []
    for design in args.design:
        deseq_anno_file = pd.read_table(os.path.join(outdir, "DE_results_"+design+"_annotated.txt")).fillna("nan")
        visual_deseq_ls.append(deseq_anno_file)
    visual_deseq = pd.concat(visual_deseq_ls, ignore_index=True)
    
    get_visualization_script(visual_deseq, args.outdir, args.starr, ori_dir)

    pybedtools.cleanup(remove_all=True)
    shutil.rmtree(temp_dir)


if __name__ == '__main__':
    args = parse_args()
    main(args)