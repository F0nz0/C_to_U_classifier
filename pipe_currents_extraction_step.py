# step one - extraction basecalling and currents/dwells features and prediction by iForestcc1 model on basecalling features only.
# for the step n° two use train_model_from_scratch.py scripts.
from datetime import datetime
import os, argparse
from pipe_basecalling import pipe_basecalling_main_function
from currents_dwell_retriever_multiprocessing import currents_dwells_retriever
from filters_currents_dwells_features import filters_currents_dwells_features
from __init__ import __version__


def pipe_currents_step1_main_function(bam_file_path, ref_path, eventalign_collapsed_folder_path, min_depth=50, threads_n = 1, aligner="minimap2", organism=None):
    '''
    Perform the first step for the pipe with currents and dwells features to extract features to train
    the CNN WaveNet model on CC and TT currents/dwells features and to predict Cproba and Tproba for each read on
    CT context. It also produce aggregated genome-space related results in addition to per-read predictions.
    Furthermore, during this first step, it will be run internally the pipe_basecalling.py and thus will be produced
    a folder with predictions (both on read and genome-space aggregated level) of CT context by the use of basecalling 
    features only via a pretrained iForestcc1 model (trained on synthetic data of curlcake run).
    '''

    # step 1: basecalling features retrieving
    pipe_basecalling_main_function(bam_file_path=bam_file_path, 
                                   ref_path=ref_path, 
                                   min_depth=min_depth,
                                   threads_n=threads_n,
                                   aligner=aligner,
                                   organism=organism)

    print("\n##################################################################################################################", flush=True)
    print(f"[{datetime.now()}] Starting extraction of currents and dwell times from eventalign collapsed reads from folder: {eventalign_collapsed_folder_path}. Progresses will be written into a log file into the output folder at the same level of the input BAM file.", flush=True)
    # step 2: currents and dwell times retrieving from eventalign.collapsed folder
    currents_dwells_retriever(bam_file_path=bam_file_path,
                              eventalign_collapsed_folder_path=eventalign_collapsed_folder_path,
                              min_depth=min_depth, 
                              threads_n=threads_n)
    
    print("\n##################################################################################################################", flush=True)
    # define expected output of previous steps to be used as input for the last step of this pipe_currents_extraction_step.py
    bam_file_name = os.path.basename(bam_file_path).split(".")[0]
    bam_file_base_path = os.path.join(os.path.dirname(bam_file_path), bam_file_name)
    output_folderpath = os.path.join(bam_file_base_path + ".basecalling_features")
    TTcontext_folder_path = os.path.join(output_folderpath, bam_file_name + ".TTcontext_reads_features_forward")
    CTcontext_folder_path = os.path.join(output_folderpath, bam_file_name + ".CTcontext_reads_features_forward_rev") # it will store also CT reads from reverse strand
    CCcontext_folder_path = os.path.join(output_folderpath, bam_file_name + ".CCcontext_reads_features_forward")
    T_C_reference_currents_dwells_features_reads = bam_file_base_path + ".T_C_reference_currents_dwells_features_reads"

    # step 3: filtering C and T ref currents/dwells using guppy basecalling features extracted in 1st step of this pipeline
    filters_currents_dwells_features(TTcontext_forward_reads_folder_path=TTcontext_folder_path,
                                     CCcontext_forward_reads_folder_path=CCcontext_folder_path,
                                     CTcontext_forward_reads_folder_path=CTcontext_folder_path,
                                     T_C_reference_currents_dwells_features_reads=T_C_reference_currents_dwells_features_reads)

    print(f"[{datetime.now()}] Iteration for currents extraction pipeline finished.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform the first step for the pipe with currents and dwells features to extract features to train the CNN WaveNet model on CC and TT currents/dwells features and to predict Cproba and Tproba for each read on CT context. It also produce aggregated genome-space related results in addition to per-read predictions. Furthermore, during this first step, it will be run internally the pipe_basecalling.py and thus will be produced a folder with predictions (both on read and genome-space aggregated level) of CT context by the use of basecalling features only via a pretrained iForest model.")

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
    
    parser.add_argument("-EC",
                        "--eventalign_collapsed",
                        required=True,
                        type=str,
                        help="--eventalign_collapsed_folder_path: \t a <str> indicating the folder path of the directory in which collapsed eventalign reads are stored together with the 0.collapsed_reads.idx index file.")

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
    eventalign_collapsed_folder_path = args.eventalign_collapsed
    depth_threshold = args.threshold
    thread_n = args.threads
    aligner = args.aligner
    organism = args.organism
    
    # print some starting info related to the C_to_U_classifier version, to the used program and to the input arguments
    print(f"[{datetime.now()}] C_to_U_classifier version: {__version__}", flush=True)
    print(f"[{datetime.now()}] Pipe Currents (and basecalling) features extraction step. Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)
    
    # launch main function of the pipe
    pipe_currents_step1_main_function(bam_file_path=bam_file_path, 
                                      ref_path=ref_path, 
                                      eventalign_collapsed_folder_path=eventalign_collapsed_folder_path, 
                                      min_depth=depth_threshold,
                                      threads_n=thread_n,
                                      aligner=aligner,
                                      organism=organism)