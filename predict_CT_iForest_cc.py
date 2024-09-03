import os, joblib, shutil
from datetime import datetime
import numpy as np
import pandas as pd
from aggregate_results import aggregate_results_genomic_dimension
import argparse
from tqdm import tqdm
import sys

def reduce_basecalling_features(dfCT):
    '''
    Function to encode and reduce basecalling features in Tqual, MeanQual, n_mismatches, n_del, n_ins 
    (it also reduce dataframe dimension).
    '''
    columns = ["region", "position", "read_name", "strand", "Tqual", "MeanQual", "n_mismatches", "n_dels", "n_ins"]
    encoded = []
    '''
    Reduce basecalling features with encoding them into 
    Tqual, MeanQual, n_mismatches, n_del, n_ins.
    '''
    with tqdm(total=dfCT.shape[0], file=sys.stdout) as pbar:
        for s in dfCT.itertuples():
            site = list(s)[1:]
            # retrieve encoded basecalling features
            vector = np.array(site[4:-1])
            Tqual = vector[3] # quality of T central base
            MeanQual = round(np.mean(np.abs(vector[vector != 0])),1) # mean quality of the interval excluding 0 qualities being these deletions and absolute values since negative values are mismatches
            n_mismatches = vector[vector<0].shape[0]
            n_dels = vector[vector==0].shape[0]
            n_ins = int(site[-1])
            output = site[:4] + [Tqual, MeanQual, n_mismatches, n_dels, n_ins]
            encoded.append(output)
            pbar.update()
    dfCT_enc = pd.DataFrame(encoded, columns=columns)
    
    return dfCT_enc



def predict_CT_reads_and_sites_by_iForestcc1_model(df_CT_native_folderpath, bam_filepath, aligner="minimap2", ref_filepath=None, organism=None):
    '''
    Function to predict from basecalling features of CT reads, T real reads 
    by the use of iForest model tranined on syntethic construct cc1.
    It will produce a new folder with a tsv for per-read predictions (df_CT_predicted.tsv file) and
    a tsv with per-position aggregated results before and after iForest filtering of CT context reads.
    '''
    
    # correct also into argparse
    available_aligners = ("minimap2", "gmap")
    
    
    # define inputs
    df_CT_native_filepath = os.path.join(df_CT_native_folderpath, "CTcontext_reads_features_forward_rev_chr_whole.tsv")
    model_folderpath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"iForest_pretrained_model")
    
    if aligner in available_aligners:
        sc_filepath = os.path.join(model_folderpath, f"cc1_wt_ko.StandardScaler.{aligner}.joblib")
        iForest_filepath = os.path.join(model_folderpath, f"cc1_wt_ko.iForest.{aligner}.joblib")
    else:
        print(f"ERROR! Please use an available aligners among the following list: {available_aligners}", flush=True)

    # retrieve sample name
    sample_name = os.path.splitext(os.path.basename(df_CT_native_folderpath))[0]

    # define output paths
    output_folderpath = os.path.join(os.path.dirname(os.path.dirname(df_CT_native_folderpath)), f"{sample_name}.model_iForest_pretrained_results")
    output_reads_filepath = os.path.join(output_folderpath, "df_CT_predicted.tsv")

    # create output folder
    if os.path.exists(output_folderpath):
        shutil.rmtree(output_folderpath)
    os.mkdir(output_folderpath)

    # retrieving CT native basecalling data (mergin to whole chr tsv if it does not exist yet).
    if not os.path.exists(df_CT_native_filepath):
        # merge all tsv files of CT basecalling features
        os.system(f"cat {df_CT_native_folderpath}/CTcontext_reads_features_forward_rev_chr*.tsv > {df_CT_native_filepath}")
    # load to pandas dataframe
    print(f"[{datetime.now()}] Loading basecalling features extracted for <{sample_name}> from folder: {df_CT_native_folderpath}", flush=True)
    df_CT_native = pd.read_table(df_CT_native_filepath, header=None)
    df_CT_native.columns = ["region", "position", "read_name", "strand"] + [f"pos{i}" for i in range(-3,4)] + ["ins_count"]

    # remove merged tsv because alredy loaded into a pandas dataframe in-memory
    os.remove(df_CT_native_filepath)

    print(f"[{datetime.now()}] Reduce basecalling features.")
    df_CT_native = reduce_basecalling_features(df_CT_native)
    
    # save native reads coordinates for calculation of native CT count and frequency before the indels and mismatches filtering needed for the iForest prediction step.
    df_CT_native_coords = df_CT_native.iloc[:,0:4].copy()
    
    # retain only reads without indels and mismatches
    df_CT_native = df_CT_native.query("n_ins <= 0").query("n_dels <= 0").query("n_mismatches <=0")

    print(f"[{datetime.now()}] Loading pretrained instances of StandardScaler and iForest model for {aligner} aligner from directory: {model_folderpath}.", flush=True)
    # load pretrained StandardScaler and iForest model from disk
    sc = joblib.load(sc_filepath)
    iForest = joblib.load(iForest_filepath)

    # split data into X and standardize them and make prediction by iForest pretrained model (-1: anomaly (real T); +1: normal observation (T false/base ONT error on CT context))
    X = df_CT_native.iloc[:,4:6].values
    print(f"[{datetime.now()}] Predicting real T and error T from native ONT data (basecalling features) by iForest pretrained model [1: T real; 0: false T/error].", flush=True)
    y_hat = iForest.predict(sc.transform(X))

    # encode predicted label (1.0 --> real T; 0.0 --> false T) and put into final read-space dataframe
    y_hat = [1.0 if (i == -1) else 0.0 for i in y_hat]
    
    # eliminating basecalling features from native dataframe and adding predicted class vector
    df_CT_native = df_CT_native.iloc[:,0:4].copy()
    df_CT_native["pred"] = y_hat
    
    # merge with native coords and set to 0 the class for reads with indels and/or mismatches not predicted by iForest model
    df_CT_native = pd.merge(df_CT_native_coords, df_CT_native, how="left", on=["region", "position", "read_name", "strand"])
    df_CT_native.fillna(0, inplace=True) # fill nan values for 'pred' columns on reads without iForest prediction (set to 0...not reliable true T reads)

    print(df_CT_native)

    # save to disk per-read predictions
    print(f"[{datetime.now()}] Saving per-read prediction on disk at path: {output_reads_filepath}", flush=True)
    df_CT_native.to_csv(output_reads_filepath, sep="\t", index=False)

    # aggregate using aggregation function and save to disk inside the same folder of df_CT_native
    aggregate_results_genomic_dimension(df_CT_predicted_filepath=output_reads_filepath, 
                                        bam_filepath=bam_filepath,
                                        ref_filepath=ref_filepath,
                                        model_type="iforest",
                                        organism=organism)


if __name__ == "__main__":

    available_aligners = ("minimap2", "gmap")

    parser = argparse.ArgumentParser(description="Function to predict CT reads basing on basecalling features and pretrained iForest model, and to aggregate them on genomic dimension. It also produce from a bam file a column with stranded depth for each site with CT reads detected and T base frequencies before and after the model correction (T reads with class:1 are those predicted as anomalies by iForest model).")
    parser.add_argument("-dfCT_native_fp",
                        required=True,
                        type=str,
                        help="-df_CT_native_folderpath: \t a <str> indicating the fullpath for the folder containing the extracted basecalling features of CT contexts reads.")
    parser.add_argument("-bam",
                        required=True,
                        type=str,
                        help="-bam_file_fullpath: \t a <str> with the full_path for the input bam file")
    parser.add_argument("-aligner",
                        required=False,
                        default="minimap2",
                        type=str,
                        help=f"-aligner_used: \t a <str> indicating the aligner used to produce BAM file for run in analysis. Available models for aligners: {available_aligners}.")

    parser.add_argument("-ref_filepath",
                        required=False,
                        default=None,
                        type="str",
                        help="-reference_filepath: \t a <str> indicating the full-paht for the reference used for aligments.")
    
    parser.add_argument("-O",
                        "--organism",
                        required=False,
                        default=None,
                        type=str,
                        help=f"-organism: \t a <str> indicating the organims analyzed in order to produce a p-value of the APOBEC1 signature resemblance after aggregation of genome space of the sites.")
    

    args = parser.parse_args()
    df_CTnative_folderpath = args.dfCT_native_fp
    bam_filepath = args.bam
    aligner = args.aligner
    ref_filepath = args.ref_filepath
    organism = args.organism
    
    predict_CT_reads_and_sites_by_iForestcc1_model(df_CT_native_folderpath=df_CTnative_folderpath, 
                                                   bam_filepath=bam_filepath, 
                                                   aligner=aligner,
                                                   ref_filepath=ref_filepath,
                                                   organism=organism)