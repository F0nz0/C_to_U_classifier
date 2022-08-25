from ast import arg
from datetime import datetime
from weakref import ref
from TT_CT_CC_basecalling_features_retriever import T_and_C_reads_features_extractor_forward_rev_main_functions
from predict_CT_iForest_cc import predict_CT_reads_and_sites_by_iForestcc1_model
import os
import argparse
from __init__ import __version__

def pipe_basecalling_main_function(bam_file_path, ref_path, min_depth=50, threads_n=1, aligner="minimap2", organism=None):
    '''
    Whole pipeline for basecalling features extraction and prediction on those via iForestcc1 pretrained model.
    '''

    # 1st step: extract basecalling features from BAM
    T_and_C_reads_features_extractor_forward_rev_main_functions(bam_file_path=bam_file_path, 
                                                                ref_path=ref_path,
                                                                min_depth=min_depth,
                                                                threads_n=threads_n)

    # define output folder path that will be the input folder path for 2nd step
    bam_file_name = os.path.basename(bam_file_path).split(".")[0]
    bam_file_base_path = os.path.join(os.path.dirname(bam_file_path), bam_file_name)
    output_folderpath = os.path.join(bam_file_base_path + ".basecalling_features")
    df_CT_native_folderpath = os.path.join(output_folderpath, bam_file_name + ".CTcontext_reads_features_forward_rev")
    # anyway the path for df_CT_native_folderpath is:
    df_CT_native_folderpath = os.path.join(output_folderpath, bam_file_name + ".CTcontext_reads_features_forward_rev")

    # 2nd step: make prediction on per-read level observation and aggregate on genome space
    predict_CT_reads_and_sites_by_iForestcc1_model(df_CT_native_folderpath=df_CT_native_folderpath,
                                                   bam_filepath=bam_file_path, 
                                                   aligner=aligner, 
                                                   ref_filepath=ref_path,
                                                   organism=organism)
    
    print(f"[{datetime.now()}] Iteration for basecalling pipeline finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Whole pipeline for basecalling features extraction and prediction on those via iForestcc1 pretrained model.")
    
    parser.add_argument("-B",
                        "--bam_filepath",
                        required=True,
                        type=str,
                        help="--bam_filepath: \t a <str> indicating the fullpath for the bam file to be used to extract basecalling features.")
    
    parser.add_argument("-R",
                        "--reference_filepath",
                        required=True,
                        type=str,
                        help="--reference_file_fullpath: \t a <str> with the full_path for the input reference fasta file used to align reads into the BAM file.")
    
    parser.add_argument("-T",
                        "--threshold", 
                        required=False,
                        type=int,
                        default=50,
                        help="--threshold_value_depth: \t a <int> indicating a threshold to filter out position with lower depth.")

    parser.add_argument("-threads", 
                        required=False,
                        type=int,
                        default=1,
                        help="-threads_number: \t a <int> indicating the number of threads to be used (it should be equal to the number of available CPUs).")
    
    parser.add_argument("-aligner",
                        required=False,
                        default="minimap2",
                        type=str,
                        help=f"-aligner_used: \t a <str> indicating the aligner used to produce BAM file for run in analysis.")
    
    parser.add_argument("-O",
                        "--organism",
                        required=False,
                        default=None,
                        type=str,
                        help=f"-organism: \t a <str> indicating the organims analyzed in order to produce a p-value of the APOBEC1 signature resemblance after aggregation of genome space of the sites.")
    

    args = parser.parse_args()
    bam_file_path = args.bam_filepath
    ref_path = args.reference_filepath
    depth_threshold = args.threshold
    threads_n = args.threads
    aligner = args.aligner
    organism = args.organism

    # print some starting info related to the C_to_U_classifier version, to the used program and to the input arguments
    print(f"[{datetime.now()}] C_to_U_classifier version: {__version__}", flush=True)
    print(f"[{datetime.now()}] Pipe Basecalling features extraction. Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launch main function of the pipe
    pipe_basecalling_main_function(bam_file_path=bam_file_path, 
                                   ref_path=ref_path, 
                                   min_depth=depth_threshold,
                                   threads_n=threads_n,
                                   aligner=aligner,
                                   organism=organism)