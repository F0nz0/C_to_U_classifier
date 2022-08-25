# import basic modules
import pandas as pd
from datetime import datetime
import pysam
import argparse
import os, sys
from utils import predict_editing_custom_thrs, df_aggr_apobec1_annot
from tqdm import tqdm

### --- utils functions --- ###
def retrieve_depth_stranded(bam_file, region, pos1based, strand):
    pos0based=pos1based-1
    for pileupcolumn in bam_file.pileup(region, pos0based, pos0based+1, truncate=True, max_depth=1000000, min_base_quality=0):
        column = pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True)
        if strand == "+":
            depth_stranded = sum(1 for b in column if b.isupper())
        elif strand == "-":
            depth_stranded = sum(1 for b in column if b.islower())
    return depth_stranded
        


def aggregate_results_genomic_dimension(df_CT_predicted_filepath, bam_filepath, Tproba_thr=0.8, ref_filepath=None, model_type=None, organism=None):
    '''
    Function to aggregates CT reads basing on prediction of iForest or CNN-WaveNet model, on genomic dimension. 
    It also produce from a bam file a column with stranded depth for each site with CT reads detected 
    and T base frequencies before and after the model correction (CT reads filtering via thresholding on 
    model probabilities output).
    '''
    # define output filepath
    output_filepath = os.path.join(os.path.dirname(df_CT_predicted_filepath), "df_CT_predicted_aggregated.tsv")

    print(f"[{datetime.now()}] Aggregating results of reads from df_CT_predicted table: {df_CT_predicted_filepath}", flush=True)
    df_CT = pd.read_table(df_CT_predicted_filepath)

    # aggregate results on native ONT data
    df_CT_aggregated_native = df_CT.groupby(["region", "position", "strand"]).size().to_frame(name = 'T_native').reset_index()

    # apply correction via a threshold on Tproba columns and aggregate results
    if "Tproba" in df_CT.columns:
        df_CT_aggregated_corrected = df_CT.query(f"Tproba > {Tproba_thr}").groupby(["region", "position", "strand"]).size().to_frame(name = 'T_corrected').reset_index()
    else:
        df_CT_aggregated_corrected = df_CT.query(f"pred > {Tproba_thr}").groupby(["region", "position", "strand"]).size().to_frame(name = 'T_corrected').reset_index()

    # merge dataframes
    df_CT_aggregated = pd.merge(df_CT_aggregated_native, 
                                df_CT_aggregated_corrected, 
                                how="left", 
                                left_on=["region", "position", "strand"],
                                right_on=["region", "position", "strand"])

    # fill with 0 sites that after correction results with nan due to no T bases on the site after model filtering
    df_CT_aggregated.fillna(0, inplace=True)

    ### retrieve depth for each position
    # open bam file with pysam
    print(f"[{datetime.now()}] Retrieving stranded depths for each modified position and calculating T base frequecy before and after model correction from bam file: {bam_filepath}", flush=True)
    bam = pysam.AlignmentFile(bam_filepath) # opening bamfile via pysam
    depths_stranded = []
    with tqdm(total=df_CT_aggregated.shape[0], file=sys.stdout) as pbar:
        for s in df_CT_aggregated.itertuples():
            d = retrieve_depth_stranded(bam, s.region, s.position, s.strand)
            depths_stranded.append(d)
            pbar.update(1)
    bam.close()
    # append depths_stranded columns to final dataframe of aggregated data
    df_CT_aggregated["depth_stranded"] = depths_stranded

    # calculate T frequencies before and after correction
    df_CT_aggregated["Tfreq_native"] = df_CT_aggregated["T_native"] / df_CT_aggregated["depth_stranded"]
    df_CT_aggregated["Tfreq_corrected"] = df_CT_aggregated["T_corrected"] / df_CT_aggregated["depth_stranded"]

    # change dtypes to int at T_corrected columns
    df_CT_aggregated["T_corrected"] = df_CT_aggregated["T_corrected"].astype("int64")
    # sort by region and position
    df_CT_aggregated.sort_values(by=["region", "position"])
    

    # if both reference and model type are provided, produce genome-space level prediction of CtoU substituition label.
    if ref_filepath != None and model_type != None:
        if model_type in ("iforest", "CNN"):
            # load custom thresholds computed on cc1 and cc2 samples
            cc1_cc2_merged_freq_thresholds_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)), 
                                                                                   "cc1_cc2_freqs_thrs", 
                                                                                   f"cc1_cc2_merged.{model_type}.freq_thresholds.tsv")
            # define basal cut-off value based on 95th (for CNN-WaveNet) percentile of CT residual errors on synthetic runs after
            # algorithms corrections
            if model_type == "iforest":
                freq_thrs = 0.01
                print(f"[{datetime.now()}] Predicting CtoU substitutions events using {model_type}.")
            elif model_type == "CNN":
                freq_thrs = 0.01
                print(f"[{datetime.now()}] Predicting CtoU substitutions events using {model_type} custom 5-mer specific cut-off values.")
            # produce label of CtoU substitutions for each sites in aggregated genome space data.
            
            df_CT_aggregated = predict_editing_custom_thrs(df_CT_aggregated,
                                                           cc1_cc2_merged_freq_thresholds_filepath,
                                                           ref_filepath,
                                                           freq_thrs)
        else:
            print(f"[{datetime.now()}] ERROR! The provided <{model_type}> is wrong to predict CtoU substitution on genome space aggregated data.", flush=True)

    # remove rows with nan on y_hat
    df_CT_aggregated.dropna(axis=0, inplace=True)

    # annotate with the apobec signature
    if ref_filepath != None and (organism in ("human", "murine")):
        df_CT_aggregated = df_aggr_apobec1_annot(df_CT_aggregated, ref_filepath, organism=organism)

    # save to csv on disk at the same path of the starting df_CT_predicted.tsv dataframe
    df_CT_aggregated.to_csv(output_filepath, sep="\t", index=False)




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Function to aggregates CT reads basing on prediction of CNN-WaveNet model, on genomic dimension. It also produce from a bam file a column with stranded depth for each site with CT reads detected and T base frequencies before and after the model correction (CT reads filtering via thresholding on model probabilities output")
    parser.add_argument("-dfCTpred",
                        required=True,
                        type=str,
                        help="-df_CT_predicted_tsv_filepaht: \t a <str> indicating the fullpath for the df_CT_predicted.tsv table produced by the model.")
    parser.add_argument("-bam",
                        required=True,
                        type=str,
                        help="-bam_file_fullpath: \t a <str> with the full_path for the input bam file")

    parser.add_argument("-Tproba_thr",
                        required=False,
                        type=float,
                        default=0.8,
                        help="-T_proba_threshold_value: \t a <float> indicating a probability threshold to filter out CT reads with lower probability given by model.")

    parser.add_argument("-ref_filepath",
                        required=False,
                        type=str,
                        default=None,
                        help="-reference_filepath: \t a <str> indicating the full-path for the reference to be used when predicting CtoU label at genome-space level.")

    parser.add_argument("-model_type",
                        required=False,
                        type=str,
                        default=None,
                        help="-model_type: \t a <str> indicating the model type used to make correction of CT frequencies, and that will be used to load specific custom threshold for CtoU genome-space predicitons.")

    args = parser.parse_args()
    df_CT_predicted_filepath = args.dfCTpred
    bam_filepath = args.bam
    Tproba_thr = args.Tproba_thr
    ref_filepath = args.ref_filepath
    model_type = args.model_type

    # launch main function
    aggregate_results_genomic_dimension(df_CT_predicted_filepath=df_CT_predicted_filepath, 
                                        bam_filepath=bam_filepath,
                                        Tproba_thr=Tproba_thr,
                                        ref_filepath=ref_filepath,
                                        model_type=model_type)