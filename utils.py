import os, sys
import numpy as np
import pandas as pd
import pysam
import matplotlib.pyplot as plt
import seaborn as sn
from tqdm import tqdm
from datetime import datetime
import logomaker
from scipy.stats import chi2
import pickle

def encode_strand(strand_str):
    # encode strand in case of numeric strand (reditools like)
    strand = int(strand_str)
    d = {0:"-", 1:"+", 2:"+"}
    converted_strand = d[strand]
    return converted_strand


def remove_snps(input_table_filepath, idx_type=1, stranded=True, pos_col=1, chr_col=0, strand_col=3, mask_filepath="/lustre/bio_running/refs/snp151.sorted.gtf.gz"):
    """
    Filter known SNP positions from a tabular file.
    """
    # open with pysam mask file
    snps = pysam.TabixFile(mask_filepath)
    
    # output file where putative SNPs were stored
    match_filepath = os.path.join(os.path.dirname(input_table_filepath), os.path.basename(input_table_filepath)+".filtered_out")
    
    # output file where positions without SNPs evidence are saved
    nomatch_filepath = os.path.join(os.path.dirname(input_table_filepath), os.path.basename(input_table_filepath)+".filtered")
        
    # open output files
    match_file = open(match_filepath, "w")
    nomatch_file = open(nomatch_filepath, "w")
    
    with open(input_table_filepath, "r") as input_file:
        for l in input_file:
            line = l.rstrip().split("\t")
            chr_ = line[chr_col]
            pos_ = int(line[pos_col])
            # convert in 1-based if input table is in 0-based format
            if idx_type == 0:
                pos_ = pos_ + 1
            if stranded == True:
                strand_ = line[strand_col]
                if strand_ in ("0","1","2"):
                    strand_ = encode_strand(strand_)
            # check if the site is into the mask file
            # fetch the mask files at current position
            TabixFetch = [str(i).split("\t") for i in snps.fetch(chr_, pos_-1, pos_, parser=pysam.asGTF())]
            if len(TabixFetch) > 0:
                # is a probable snp
                if stranded == True:
                    matched = False
                    for tabix in TabixFetch:
                        if tabix[6] == strand_:
                            matched = True
                    if not matched:
                        nomatch_file.write(l)
                    elif matched:
                        match_file.write(l)
                elif stranded == False:
                    match_file.write(l)
            else:
                # is not a snp
                nomatch_file.write(l)
    match_file.close()
    nomatch_file.close()
    snps.close()


def get_rev_compl(seq):
    seq = seq.upper()
    alphabet = {"A":"T", "T":"A", "C":"G", "G":"C"}
    rev = seq[::-1]
    rev_compl = "".join([alphabet[b] for b in rev])
    return rev_compl


def fastqRNAtoDNA(input_fastq_filepath):
    '''
    Convert U based RNA FastQ to T based DNA FastQ.
    '''
    output_fastq_filepath = os.path.join(os.path.dirname(input_fastq_filepath), 
                                         f"{os.path.basename(input_fastq_filepath).split('.fastq')[0]}.dna.fastq")

    # calculate n of rows into input fastq
    n_rows=0
    with open(input_fastq_filepath, "r") as inp:
        for l in inp:
            n_rows += 1

    print(f"[{datetime.now()}] Converting input fastq {input_fastq_filepath} to dna fastq (UtoT replacement).\n[{datetime.now()}] File will be saved as {output_fastq_filepath}")
    with tqdm(total=n_rows, file=sys.stdout) as pbar:
        with open(output_fastq_filepath, "w") as out:
            with open(input_fastq_filepath, "r") as inp:
                counter = 1
                expected_counter = 2
                for l in inp:
                    line = l
                    if counter == expected_counter:
                        line = l.replace("U", "T")
                        expected_counter += 4
                    out.write(line)
                    counter += 1
                    pbar.update(1)


def plot_frequencies(reference_filepath, 
                     df_aggregated_filepath, 
                     title=None, 
                     freq_threshold=0.02, 
                     native=True, 
                     corrected=True, 
                     save_folderpath=None, 
                     log=True, 
                     offset=0, 
                     native_color="gray",
                     corrected_color="black",
                     native_size=20,
                     corrected_size=8,
                     strand=None,
                     dpi=300):
    # load results of predicted 
    df = pd.read_table(df_aggregated_filepath)
    if strand:
        # filter if strand preferences is expressed, else use both plus and minus strand
        df = df.query(f"strand == '{strand}'")
    
    # go to psedo frequency to avoid 0 freq. sites.
    if log==True:
        df["Tfreq_native"] = df["Tfreq_native"] + offset
        df["Tfreq_corrected"] = df["Tfreq_corrected"] + offset
        freq_threshold += offset
    
    # retrieve min and max coords per contig
    ref = pysam.FastaFile(reference_filepath)
    x_axes = {}
    
    print(f"[{datetime.now()}] Retrieving Contigs and lengths...")
    for contig in ref.references:
        df_positions = pd.DataFrame(np.arange(1,ref.get_reference_length(contig)+1), columns=["position"])
        
        print(f"[{datetime.now()}] Merging dataframes for region {contig}...")
        df_merged = pd.merge(df_positions, df.query(f"region == '{contig}'"), how="left", on=["position"]).fillna(0)
        
        print(f"[{datetime.now()}] Merging finished....start plotting...")
        plt.figure(figsize=(15,5))
        if native == True:
            plt.scatter(df_merged["position"], 
                        df_merged["Tfreq_native"],
                        label="ONT native", c=native_color, s=native_size)
                     
        if corrected == True:
            plt.scatter(df_merged["position"], 
                        df_merged["Tfreq_corrected"],
                        label="ONT corrected", c=corrected_color, s=corrected_size)
            plt.xlabel("Genomic Position (bp)")
        
        plt.ylabel("Tfreq")
        plt.title(f"{title} - {contig}")
        plt.legend(loc="upper right")
        plt.axhline(y=freq_threshold, color='r', linestyle='--', linewidth=1)
        
        if log == True:
            plt.yscale("log")
            yticks = plt.yticks()[0]
            plt.ylabel(f"Log (T freq)")
            plt.yticks(yticks, np.round(yticks-offset, 2))
            
        
        sn.despine(offset=0.00001, trim=True, )
        # to plot better figures
        plt.tight_layout()
        
        # save on disk if required with a given file-path to save in.
        if save_folderpath != None:
            plt.savefig(os.path.join(save_folderpath, f"{title} - {contig}.tiff"), dpi=dpi)
        
        # show figure eventually
        plt.show()
    
    # close reference
    ref.close()


# function to plot a complete confusion matrix with common metrics.
def make_confusion_matrix(cf,
                          group_names=None,
                          categories='auto',
                          count=True,
                          percent=True,
                          cbar=True,
                          xyticks=True,
                          xyplotlabels=True,
                          sum_stats=True,
                          figsize=None,
                          cmap='Blues',
                          title=None,
                          path=None):
    '''
    ###############################################################################################
    CITATION: taken from: https://github.com/DTrimarchi10/confusion_matrix/blob/master/cf_matrix.py
    ###############################################################################################
    
    This function will make a pretty plot of an sklearn Confusion Matrix cm using a Seaborn heatmap visualization.
    Arguments
    ---------
    cf:            confusion matrix to be passed in
    group_names:   List of strings that represent the labels row by row to be shown in each square.
    categories:    List of strings containing the categories to be displayed on the x,y axis. Default is 'auto'
    count:         If True, show the raw number in the confusion matrix. Default is True.
    normalize:     If True, show the proportions for each category. Default is True.
    cbar:          If True, show the color bar. The cbar values are based off the values in the confusion matrix.
                   Default is True.
    xyticks:       If True, show x and y ticks. Default is True.
    xyplotlabels:  If True, show 'True Label' and 'Predicted Label' on the figure. Default is True.
    sum_stats:     If True, display summary statistics below the figure. Default is True.
    figsize:       Tuple representing the figure size. Default will be the matplotlib rcParams value.
    cmap:          Colormap of the values displayed from matplotlib.pyplot.cm. Default is 'Blues'
                   See http://matplotlib.org/examples/color/colormaps_reference.html
                   
    title:         Title for the heatmap. Default is None.
    '''


    # CODE TO GENERATE TEXT INSIDE EACH SQUARE
    blanks = ['' for i in range(cf.size)]

    if group_names and len(group_names)==cf.size:
        group_labels = ["{}\n".format(value) for value in group_names]
    else:
        group_labels = blanks

    if count:
        group_counts = ["{0:0.0f}\n".format(value) for value in cf.flatten()]
    else:
        group_counts = blanks

    if percent:
        group_percentages = ["{0:.2%}".format(value) for value in cf.flatten()/np.sum(cf)]
    else:
        group_percentages = blanks

    box_labels = [f"{v1}{v2}{v3}".strip() for v1, v2, v3 in zip(group_labels,group_counts,group_percentages)]
    box_labels = np.asarray(box_labels).reshape(cf.shape[0],cf.shape[1])


    # CODE TO GENERATE SUMMARY STATISTICS & TEXT FOR SUMMARY STATS
    if sum_stats:
        #Accuracy is sum of diagonal divided by total observations
        accuracy  = np.trace(cf) / float(np.sum(cf))

        #if it is a binary confusion matrix, show some more stats
        if len(cf)==2:
            #Metrics for Binary Confusion Matrices
            precision = cf[1,1] / sum(cf[:,1])
            recall    = cf[1,1] / sum(cf[1,:])
            f1_score  = 2*precision*recall / (precision + recall)
            stats_text = "\n\nAccuracy={:0.3f}\nPrecision={:0.3f}\nRecall={:0.3f}\nF1 Score={:0.3f}".format(
                accuracy,precision,recall,f1_score)
        else:
            stats_text = "\n\nAccuracy={:0.3f}".format(accuracy)
    else:
        stats_text = ""


    # SET FIGURE PARAMETERS ACCORDING TO OTHER ARGUMENTS
    if figsize==None:
        #Get default figure size if not set
        figsize = plt.rcParams.get('figure.figsize')

    if xyticks==False:
        #Do not show categories if xyticks is False
        categories=False


    # MAKE THE HEATMAP VISUALIZATION
    plt.figure(figsize=figsize)
    sn.heatmap(cf,annot=box_labels,fmt="",cmap=cmap,cbar=cbar,xticklabels=categories,yticklabels=categories)

    if xyplotlabels:
        plt.ylabel('True label')
        plt.xlabel('Predicted label' + stats_text)
    else:
        plt.xlabel(stats_text)
    
    if title:
        plt.title(title)
    
    if path:
        plt.tight_layout()
        plt.savefig(path)


# plot different version of consensus graphs
def create_sequence_logo(df, color_scheme=None, ax=None, title=None):
    '''
    Produce Consensus plot.
    '''
    if color_scheme == None:
        color_scheme={"*":"black",
                      "T": "red",
                      "A":"green",
                      "C":"blue",
                      "G":"orange"}
    else:
        color_scheme = color_scheme
    
    if ax == None:
        crp_logo = logomaker.Logo(df,
                                  shade_below=.5,
                                  fade_below=.5,
                                  font_name='Arial Rounded MT Bold', 
                                  color_scheme=color_scheme)
        # style using Logo methods
        crp_logo.style_spines(visible=False)
        crp_logo.style_spines(spines=['left', 'bottom'], visible=True)

        # style using Axes methods
        crp_logo.ax.set_ylabel("Frequency", labelpad=-1)
        crp_logo.ax.xaxis.set_ticks_position('none')
        crp_logo.ax.xaxis.set_tick_params(pad=-1)
        plt.title(title)


def compute_apobec1_signature_pvalue(ref_filepath, region, pos1based, strand, verbosity=0, H0_filepath=None, H1_filepath=None, organism=None):
    '''
    Perform a Likelihood Ratio Test followed by Chi-squared test to compute associated p-values 
    between:
        - H0 model: product of the base frequency observed around non edited positions (with support of Wt and Ko murine illumina samples)
        - H1 model: product of the base frequency observed around edited position (with support of Wt and Ko murine illumina samples)
    '''
    ref = pysam.FastaFile(ref_filepath)
    if H0_filepath == None:
        H0_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"apobec1_models", f"HO_No_Editing_Sites_Consensus.apobec1.{organism}.txt")
        
    if H1_filepath == None:
        H1_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"apobec1_models", f"H1_Editing_Sites_Consensus.apobec1.{organism}.txt")
        
    h0 = pd.read_table(H0_filepath, index_col=0)
    h1 = pd.read_table(H1_filepath, index_col=0)
    
    if verbosity > 2:
        print("H0 model matrix")
        print(h0)
        print()
        print("H1 model matrix")
        print(h1)
    
    if h0.shape[1] == h1.shape[1]:
        interval = int((h0.shape[1]-1)/2)
        
    pos0based=pos1based-1
    start = pos0based-interval
    stop = pos0based+interval+1
    
    # get query position sequence of reference
    # assert if strand is correct (either + or -)
    if not (strand == "+" or strand == "-"):
        print("Ops! Strange strand!")
        
    # load reference sequence interval around central base position
    reference = ref.fetch(region, start, stop)
    
    if strand == "-":
        reference = get_rev_compl(reference)
    try:
        if reference[interval] == "C":
            if len(reference) == (interval*2)+1:
                h0_probs = []
                h1_probs = []
                for rel_pos, ref_base in enumerate(reference):
                    if verbosity > 1:
                        print(rel_pos, ref_base, h0.loc[ref_base, str(rel_pos)], h1.loc[ref_base, str(rel_pos)])
                    h0_probs.append(h0.loc[ref_base, str(rel_pos)])
                    h1_probs.append(h1.loc[ref_base, str(rel_pos)])

                h0_probs_prod = np.prod(h0_probs)
                h1_probs_prod = np.prod(h1_probs)
                llr = -2*np.log(h0_probs_prod/h1_probs_prod)
                pvalue = chi2.sf(llr, 1)
                if verbosity > 0:
                    print("H0 product:", h0_probs_prod)
                    print("H1 product:", h1_probs_prod)
                    print("llr:", llr)
                    print("P-value via Chi-2 test:", pvalue)

                ref.close()

                return pvalue
            else:
                print(f"[{datetime.now()}] OPS! Given position {region}:{pos1based} ({strand}) generated an interval too short! Something got wrong!")
                return None
        else:
            print(f"[{datetime.now()}] OPS! Wrong position {region}:{pos1based} ({strand}) provided: No C central base! Base found: {reference[interval]}")
            return None
    except KeyError:
        print(f"[{datetime.now()}] OPS! An error occured for position {region}:{pos1based} ({strand}) provided. It was not possible to compute APOBEC1 signature pvalue. None was returned for this position.")
        return None


def df_aggr_apobec1_annot(df_aggr, ref_filepath, organism=None):
    '''
    Function to annotate with apobec1 signature p-value the genome
    space aggregated dataset after model correction.
    '''
    # define presumed h1 and h0 model the compute_apobec1_signature_pvalue function is going to use.

    H0_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"apobec1_models", f"HO_No_Editing_Sites_Consensus.apobec1.{organism}.txt")        
    H1_filepath = os.path.join(os.path.dirname(os.path.realpath(__file__)),"apobec1_models", f"H1_Editing_Sites_Consensus.apobec1.{organism}.txt")

    print(f"[{datetime.now()}] Annotating aggregated sites with APOBEC1 signature for organism {organism} using H1/H0 models from {H1_filepath} and {H0_filepath}, respectively.", flush=True)
    df_aggr = df_aggr.copy()
    p_values = []
    with tqdm(total=df_aggr.shape[0], file=sys.stdout) as pbar:
        for s in df_aggr.itertuples():
            p_values.append(compute_apobec1_signature_pvalue(ref_filepath, s[1], s[2], s[3], organism=organism))
            pbar.update(1)
    df_aggr["apobec1_pvalue"] = p_values
    return df_aggr

        
def predict_editing_custom_thrs(dfCTaggr, custom_thrs_filepath, ref_filepath, min_thr):
    '''
    Function to retrieve 5mers of each site of the aggregated CT putative edited positions
    and to use a custom threshold for each 5mer computed from corrected model trained on 
    cc1 sample.
    '''
    # open input needed files
    # open custom thrs file only for cnn model. iForest doesn't need it since it was more stringent.
    if not "iforest" in custom_thrs_filepath:
        custom_thrs = pd.read_table(custom_thrs_filepath, index_col=0)
    refs = []
    y_hat = [] # predictions
    ref = pysam.FastaFile(ref_filepath)
    dfCTaggr = dfCTaggr
    
    # compute 5mer for each site inside dfCTaggregated table
    refs = []
    ref = pysam.FastaFile(ref_filepath)
    with tqdm(total=dfCTaggr.shape[0], file=sys.stdout) as pbar:
        for s in dfCTaggr.itertuples():
            site = f"{s.region}:{s.position}{s.strand}"
            pos1based=s.position
            pos0based=pos1based-1
            start=pos0based-2
            stop=pos0based+3
            strand = s.strand
            reference = ref.fetch(s.region, start, stop)
            if strand == "+":
                refs.append(reference)
            elif strand == "-":
                reference = get_rev_compl(reference)
                refs.append(reference)
            else:
                # stange strand. append unkwnown value as "."
                refs.append(".")
            
            # if it is iforest model no need for 5mer specific threshold
            if "iforest" in custom_thrs_filepath:
                if s.Tfreq_corrected > min_thr:
                    y_hat.append(1)
                else:
                    y_hat.append(0)
                pbar.update(1)
            else:
                # retrieve theshold from custom_thrs dataframe and predict if the current site is predicted or not for CNN-Wavenet model
                try:
                    thr = custom_thrs.loc[reference].values[0]
                    if thr <= min_thr:
                        thr = min_thr
                        if s.Tfreq_corrected > thr:
                            y_hat.append(1)
                        else:
                            y_hat.append(0)
                    else:
                        if s.Tfreq_corrected > thr:
                            y_hat.append(1)
                        else:
                            y_hat.append(0)
                    pbar.update(1)
                except Exception as error:
                    print(f"[{datetime.now()}] PROBLEMS ON SITE:", s, "5-mer:", error, "(it will be discarded)\n")
                    y_hat.append(np.nan)
                    pbar.update(1)

    ref.close()
    dfCTaggr["5mer"] = refs
    dfCTaggr["y_hat"] = y_hat
    
    return dfCTaggr