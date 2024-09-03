'''
VERSION TO EXTRACT TT context, CT context and CC contex reads.
Here is implemented a features extractor from bam files that retrieve basecalling features at 
read-level for sites with a depth higher than threshold on TT, CC and CT contexts 
(reads calling a T on T, a C on C or C of reference, respectively).
Features extracted for each site are encoded into a single vector at per-read level, positional 
information of matches, mismatches, basecalling qualities, deletions and insertions.
Here an explicative example of features extracted from a CT context read:

Fields:            contig pos_1_based  read_name  strand  pos-3  pos-2  pos-1  pos0(site of interest)  pos+1  pos+2  pos+3  ins_count
Ref. seq.:         chr1   654613       c465-a4..  +       A      C      G      C  (always a C)         C      T      C      
Read Seq.:                                                A      G+1    C      T  (always a T)         C+2    *      C         
Basecall Qual:                                            30     25     15     22 (always positive)    12     *      7       
Enc.Vector:        chr1   654613       c465-a4..  +       30    -25    -15     22 (always positive)    12     0      7      2

Here matching bases have positive quality values while mismatches have negative quality values. 
Deletion are encoded as 0 basecalling quality value. The central position (pos0) is insted alway 
positive since this features are used to build an ML model to predict real T from false T in CT contexts, 
thus there is the need to have positive values also when on CT contexts there should be negative values 
since the CT mismatch. The ins_count features counts how many position contains at least one insertion. 
The fields "contig" and "pos_1_based" indicate the coordinates of "pos0" base with respect to the reference.
This version also retrive intervals of sites related to reads mapping on reverse (-) strand and having a 
CT substituition (i.e. GA substitutions).
'''

# importing basic modules
from distutils.dir_util import copy_tree
import subprocess, os, shutil, traceback
import numpy as np
import pandas as pd
import pysam
from datetime import datetime
from multiprocessing import Queue, Process
import argparse
import pandas.io.common

# ---- utils functions ---- #
def retrieve_basecalling_tables(bam_file, contig, site_0_based):
    region_seqs = []
    region_quals = []
    region_qnames = []
    if site_0_based-4 >= 0:
        for pileupcolumn in bam_file.pileup(contig, site_0_based-4, site_0_based+4, truncate=True, max_depth=1000000, min_base_quality=0):
            try:
                region_seqs.append(pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True))
                region_qnames.append(pileupcolumn.get_query_names())
                region_quals.append(pileupcolumn.get_query_qualities())
            except ValueError:
                print(f"[{datetime.now()}] Error! At position {contig}:{site_0_based+1}. Error: {traceback.print_exc()}", flush=True)

    
    # extract called bases
    region_df_seqs = []
    region_df_quals = []
    for n, s in zip(region_qnames, region_seqs):
        region_df_seqs.append(pd.DataFrame(s, index=n))
    try:
        region_df_seqs = pd.concat(region_df_seqs, axis=1)
        region_df_seqs.columns = list(range(site_0_based-3, site_0_based+5)) # columns in format 1-based

        # extract insertions
        region_df_ins = region_df_seqs.copy()
        region_df_ins.columns = region_df_ins.columns + 1
        region_df_ins.drop(columns=region_df_ins.columns[-1], inplace=True)
        region_df_ins = region_df_ins.applymap(lambda x: np.nan if not pd.notnull(x) else 1 if ("+" in x) else 0)

        # drop the first column-position in region_df_seqs used only for insertion calculation
        region_df_seqs.drop(columns=region_df_seqs.columns[0], inplace=True)

        # eliminate +/- terms if present for each position into region_df_seqs
        region_df_seqs = region_df_seqs.applymap(lambda x: x[0] if pd.notnull(x) else np.nan)

        # extract qualities
        for n, s in zip(region_qnames, region_quals):
            region_df_quals.append(pd.DataFrame(s, index=n))
        region_df_quals = pd.concat(region_df_quals, axis=1)
        region_df_quals.columns = list(range(site_0_based-3, site_0_based+5)) # coluns in format 1-based
        region_df_quals.drop(columns=region_df_quals.columns[0], inplace=True)

        # let's eliminates quality scores for reads and positions with deletion
        region_df_quals = region_df_quals - region_df_quals[(region_df_seqs == "*")].fillna(0)
        return region_df_seqs, region_df_ins, region_df_quals
    except:
        return pd.DataFrame([]), pd.DataFrame([]), pd.DataFrame([])
    

def create_final_mask(mask):
    final_mask = []
    for i in mask:
        if i == True:
            final_mask.append(1)
        elif i == False:
            final_mask.append(-1)
    return final_mask


def T_and_C_reads_features_extractor_forward_rev(base_bam_file_path, ref_path, chrom, TTcontext_folder_path, CTcontext_folder_path, CCcontext_folder_path, min_depth=50):
    '''
    Function to retrieve and extracts basecalling features from one of the splitted bam files T reads on 
    TT and CT contexts on references (and CC context reads). It iterate over the depth file limiting 
    the search only to positions with a depth higher than the min_depth argument. 
    To note: for reverse (-) strand, only CT reads will be retrieved! (GA substituition from bam file) and
    will be saved together with forward (+) CT reads into the same output file.
    '''

    # create path to output folders
    output_file_TTcontext = open(os.path.join(TTcontext_folder_path, "TTcontext_reads_features_forward_"+chrom+".tsv"), "w")
    output_file_CTcontext = open(os.path.join(CTcontext_folder_path, "CTcontext_reads_features_forward_rev_"+chrom+".tsv"), "w")
    output_file_CCcontext = open(os.path.join(CCcontext_folder_path, "CCcontext_reads_features_forward_"+chrom+".tsv"), "w")
    
    start_time = datetime.now()
    print(f"[{datetime.now()}] [Global {chrom} thread message] Starting Processing on {chrom}", flush=True)

    bam_file = pysam.AlignmentFile(f"{os.path.splitext(base_bam_file_path)[0]}.{chrom}.sorted.bam")
    ref = pysam.FastaFile(ref_path)
    pos_well_covered = pd.read_table(f"{os.path.splitext(base_bam_file_path)[0]}.{chrom}.sorted.bam.depth", header=None)
    pos_well_covered.columns=["contig", "pos_1_based", "depth"]
    pos_well_covered = pos_well_covered[pos_well_covered["depth"] >= min_depth]
    pos_well_covered.reset_index(inplace=True, drop=True)
    
    # retrieve total number of sites with a T or a C on reference (also reverse)
    total_TTcontext_T_bases_counter = 0
    total_CTcontext_T_bases_counter = 0
    total_CTcontext_rev_T_bases_counter = 0
    for pos in pos_well_covered.itertuples():
        reference_base = ref.fetch(pos[1], pos[2]-1, pos[2])
        if reference_base == "T":
            total_TTcontext_T_bases_counter += 1
        elif reference_base == "C":
            total_CTcontext_T_bases_counter += 1
        elif reference_base == "G":
            total_CTcontext_rev_T_bases_counter += 1
    print(f"[{datetime.now()}] [Forward {chrom} thread message] Total TT context sites on forward strand with a depth higher than {min_depth} found: {total_TTcontext_T_bases_counter}", flush=True)
    print(f"[{datetime.now()}] [Forward {chrom} thread message] Total CC/CT context sites on forward strand with a depth higher than {min_depth} found: {total_CTcontext_T_bases_counter}", flush=True)
    print(f"[{datetime.now()}] [Reverse {chrom} thread message] Total CC/CT context sites on reverse strand with a depth higher than {min_depth} found: {total_CTcontext_rev_T_bases_counter}", flush=True)
    
    # retrieve for the same sites basecalling features for the dataset
    print(f"\n[{datetime.now()}] [Global {chrom} thread message] Starting retrieving basecalling features of sites with T or C on reference from {bam_file.filename}", flush=True)
    TTcontext_T_bases_counter = 0
    CTcontext_T_bases_counter = 0
    CTcontext_rev_T_bases_counter = 0
    CCcontext_T_bases_counter = 0
    evaluated_sites_counter = 0
    for pos in pos_well_covered.itertuples():
        reference_base = ref.fetch(pos[1], pos[2]-1, pos[2])
        
        # retrieving TT context reads features
        if reference_base == "T":
            # some print to monitor the process progress
            evaluated_sites_counter += 1
            if evaluated_sites_counter % 1000 == 0:
                print(f"[{datetime.now()}] [Global {chrom} thread message] Total Sites currently evaluated: {evaluated_sites_counter} on a total of {total_TTcontext_T_bases_counter+total_CTcontext_T_bases_counter+total_CTcontext_rev_T_bases_counter}. Elapsed time {datetime.now() - start_time}", flush=True)
            seqs, ins, quals = retrieve_basecalling_tables(bam_file=bam_file, contig=pos[1], site_0_based=pos[2]-1)
            # selecting only reads with T mapped on the central position
            if not seqs.empty and seqs.shape[1] == 7:
                seqs = seqs[seqs.iloc[:,3] == "T"]
                quals = quals.loc[seqs.index]
                # drop reads that don't covers the whole region
                seqs.dropna(inplace=True)
                quals.dropna(inplace=True)
                ins = ins.loc[seqs.index]
                if not seqs.empty:
                    TTcontext_T_bases_counter += 1
                    if TTcontext_T_bases_counter % 2000 == 0:
                        print(f"[{datetime.now()}] [Forward {chrom} thread message] Positons with T reads on TT context evaluated: {TTcontext_T_bases_counter}. Elapsed Time: {datetime.now() - start_time}", flush=True)
                    # retrieve reference context to evaluate mismatches
                    ref_context = ref.fetch(pos[1], pos[2]-4, pos[2]+3)
                    for read in seqs.itertuples():
                        read_final = []
                        mask = np.array(list(read[1:])) == np.array(list(ref_context))
                        q = quals.loc[read[0]].values
                        i = ins.loc[read[0]].sum()
                        read_final = q * create_final_mask(mask)
                        read_final = np.append(read_final, i)
                        read_final = [pos[1], str(pos[2]), read[0], "+"] + [str(i) for i in read_final.tolist()]
                        line = "\t".join(read_final) + "\n"
                        output_file_TTcontext.write(line)
        
        # retrieving CT anc CC context reads features
        # TO NOTE: here the central T base will be ever positive because we don't know if it is an editing event
        elif reference_base == "C":
            # some print to monitor the process progress
            evaluated_sites_counter += 1
            if evaluated_sites_counter % 1000 == 0:
                print(f"[{datetime.now()}] [Global {chrom} thread message] Total Sites currently evaluated: {evaluated_sites_counter} on a total of {total_TTcontext_T_bases_counter+total_CTcontext_T_bases_counter+total_CTcontext_rev_T_bases_counter}. Elapsed time {datetime.now() - start_time}", flush=True)
            seqs, ins, quals = retrieve_basecalling_tables(bam_file=bam_file, contig=pos[1], site_0_based=pos[2]-1)
            # selecting only reads with T mapped on the central position
            if not seqs.empty and seqs.shape[1] == 7:
                # retrieve CT context reads
                seqs_T = seqs[seqs.iloc[:,3] == "T"].copy()
                quals_T = quals.loc[seqs.index].copy()
                # drop reads that don't covers the whole region
                seqs_T.dropna(inplace=True)
                quals_T.dropna(inplace=True)
                ins_T = ins.loc[seqs.index].copy()
                if not seqs_T.empty:
                    CTcontext_T_bases_counter += 1
                    if CTcontext_T_bases_counter % 2000 == 0:
                        print(f"[{datetime.now()}] [Forward {chrom} thread message] Positons with T reads on CT context evaluated: {CTcontext_T_bases_counter}. Elapsed Time: {datetime.now() - start_time}", flush=True)
                    # retrieve reference context to evaluate mismatches
                    ref_context = ref.fetch(pos[1], pos[2]-4, pos[2]+3)
                    for read in seqs_T.itertuples():
                        read_final = []
                        mask_T = np.array(list(read[1:])) == np.array(list(ref_context))
                        q_T = quals_T.loc[read[0]].values
                        i_T = ins_T.loc[read[0]].sum()
                        read_final = q_T * create_final_mask(mask_T)
                        read_final[3] = read_final[3] * -1 # because it will ever be minus since it is on CT context. Positivize it to get possibility to compair to TT context reads
                        read_final = np.append(read_final, i_T)
                        read_final = [pos[1], str(pos[2]), read[0], "+"] + [str(i) for i in read_final.tolist()]
                        line = "\t".join(read_final) + "\n"
                        output_file_CTcontext.write(line)
                
                # retrieve CC context reads
                seqs_C = seqs[seqs.iloc[:,3] == "C"].copy()
                quals_C = quals.loc[seqs.index].copy()
                # drop reads that don't covers the whole region
                seqs_C.dropna(inplace=True)
                quals_C.dropna(inplace=True)
                ins_C = ins.loc[seqs.index].copy()
                if not seqs_C.empty:
                    CCcontext_T_bases_counter += 1
                    if CCcontext_T_bases_counter % 2000 == 0:
                        print(f"[{datetime.now()}] [Forward {chrom} thread message] Positons with C reads on CC context evaluated: {CCcontext_T_bases_counter}. Elapsed Time: {datetime.now() - start_time}", flush=True)
                    # retrieve reference context to evaluate mismatches
                    ref_context = ref.fetch(pos[1], pos[2]-4, pos[2]+3)
                    for read in seqs_C.itertuples():
                        read_final = []
                        mask_C = np.array(list(read[1:])) == np.array(list(ref_context))
                        q_C = quals_C.loc[read[0]].values
                        i_C = ins_C.loc[read[0]].sum()
                        read_final = q_C * create_final_mask(mask_C)
                        read_final = np.append(read_final, i_C)
                        read_final = [pos[1], str(pos[2]), read[0], "+"] + [str(i) for i in read_final.tolist()]
                        line = "\t".join(read_final) + "\n"
                        output_file_CCcontext.write(line)

        # retrieving CT context reads features from reverse strand (see as G)
        # TO NOTE: here the central T base will be ever positive because we don't know if it is an editing event
        elif reference_base == "G":
            evaluated_sites_counter += 1
            if evaluated_sites_counter % 1000 == 0:
                print(f"[{datetime.now()}] [Global {chrom} thread message] Total Sites currently evaluated: {evaluated_sites_counter} on a total of {total_TTcontext_T_bases_counter+total_CTcontext_T_bases_counter+total_CTcontext_rev_T_bases_counter}. Elapsed time {datetime.now() - start_time}", flush=True)
            seqs, ins, quals = retrieve_basecalling_tables(bam_file=bam_file, contig=pos[1], site_0_based=pos[2]-1)
            # selecting only reads with T mapped on the central position on reverse strand (that are seen as minus a char by pysam)
            if not seqs.empty and seqs.shape[1] == 7:
                # retrieve CT context reads
                seqs_T = seqs[seqs.iloc[:,3] == "a"].copy()
                quals_T = quals.loc[seqs.index].copy()
                # drop reads that don't covers the whole region
                seqs_T.dropna(inplace=True)
                quals_T.dropna(inplace=True)
                ins_T = ins.loc[seqs.index].copy()
                if not seqs_T.empty:
                    CTcontext_rev_T_bases_counter += 1
                    if CTcontext_rev_T_bases_counter % 2000 == 0:
                        print(f"[{datetime.now()}] [Reverse {chrom} thread message] Positons with T (a) reads on CT (ga) context evaluated: {CTcontext_rev_T_bases_counter}. Elapsed Time: {datetime.now() - start_time}", flush=True)
                    # retrieve reference context to evaluate mismatches
                    # (convert to lower-case because pysam give back lower case string for reads mapping on minus strand of genome)
                    ref_context = (ref.fetch(pos[1], pos[2]-4, pos[2]+3)).lower()
                    for read in seqs_T.itertuples():
                        read_final = []
                        mask_T = np.array(list(read[1:])) == np.array(list(ref_context))
                        q_T = quals_T.loc[read[0]].values
                        i_T = ins_T.loc[read[0]].sum()
                        read_final = q_T * create_final_mask(mask_T)
                        read_final[3] = read_final[3] * -1 # because it will ever be minus since it is on CT context. Positivize it to get possibility to compair to TT context reads
                        # invert since it is on reverse strand
                        read_final = read_final[::-1]
                        # append insertions information
                        read_final = np.append(read_final, i_T)
                        # append coordinates 
                        read_final = [pos[1], str(pos[2]), read[0], "-"] + [str(i) for i in read_final.tolist()]
                        line = "\t".join(read_final) + "\n"
                        output_file_CTcontext.write(line)



    output_file_TTcontext.close()
    output_file_CTcontext.close()
    output_file_CCcontext.close()
    bam_file.close()
    
    stop_time = datetime.now()
    print(f"\n[{datetime.now()}] [Forward {chrom} message] Total number of T reads successfully extracted from TTcontext sites: {TTcontext_T_bases_counter}", flush=True)
    print(f"[{datetime.now()}] [Forward {chrom} message] Total number of T reads successfully extracted from CTcontext sites: {CTcontext_T_bases_counter}", flush=True)
    print(f"[{datetime.now()}] [Forward {chrom} message] Total number of C reads successfully extracted from CCcontext sites: {CCcontext_T_bases_counter}", flush=True)
    print(f"[{datetime.now()}] [Reverse {chrom} message] Total number of T reads successfully extracted from CTcontext sites: {CTcontext_rev_T_bases_counter}", flush=True)
    print(f"[{datetime.now()}] [Global {chrom} message] Elapsed time: {stop_time - start_time}", flush=True)

    
# main function retriever for TT and CT context T reads features for sites with depth > min_depth
def T_and_C_reads_features_extractor_forward_rev_main_functions(bam_file_path, ref_path, output_folderpath=None, threads_n=1, min_depth=50):
    # defining paths to files and opening them
    bam_file_path = bam_file_path
    ref_path = ref_path

    # retrieve bam file name and bam file path without extension
    bam_file_name = os.path.splitext(os.path.basename(bam_file_path))[0]
    bam_file_base_path = os.path.join(os.path.dirname(bam_file_path), bam_file_name)
    
    # open reference file
    ref = pysam.Fastafile(ref_path)

    start_time = datetime.now()

    # retrieve reference contigs
    contigs = []
    for c in ref.references:
        if "chr" in c:
            contigs.append(c)
    # if contings are not denominated as chr1 ....chrn thus it will takes all contigs
    if len(contigs) != 0:
        print(f"\n[{datetime.now()}] [Main thread message] Retrieved contigs: {contigs}", flush=True)
    else:
        print(f"\n[{datetime.now()}] [Main thread message] WARNING: No retrieved contigs with 'chr' string. Retrying with no filters for fasta file contigs.", flush=True)
        contigs = []
        for c in ref.references:
            contigs.append(c)
        print(f"\n[{datetime.now()}] [Main thread message] Retrieved contigs (without 'chr' string filter): {contigs}", flush=True)

    # split bam files file by contigs and filter not primary and supplementary alignments
    for chrom in contigs:
        print(f"[{datetime.now()}] [Main thread message] Splitting, sorting and indexing bam file {bam_file_path} on {chrom}.", flush=True)
        subprocess.run(f"samtools view -b -F 2304 {bam_file_path} {chrom}".split(" "), stdout=open(f"{bam_file_base_path}.{chrom}.bam", "w"))
        subprocess.run(f"samtools sort {bam_file_base_path}.{chrom}.bam -o {bam_file_base_path}.{chrom}.sorted.bam".split(" "))
        subprocess.run(f"samtools index {bam_file_base_path}.{chrom}.sorted.bam".split(" "))
        print(f"[{datetime.now()}] [Main thread message] Removing unsorted bam file for {chrom}", flush=True)
        os.remove(f"{bam_file_base_path}.{chrom}.bam")
        print(f"[{datetime.now()}] [Main thread message] Producing depth file for {chrom}", flush=True)
        subprocess.run(f"samtools depth {bam_file_base_path}.{chrom}.sorted.bam".split(" "), stdout=open(f"{bam_file_base_path}.{chrom}.sorted.bam.depth", "w"))
    
    if output_folderpath == None:
        output_folderpath = os.path.join(bam_file_base_path + ".basecalling_features")
    TTcontext_folder_path = os.path.join(output_folderpath, bam_file_name + ".TTcontext_reads_features_forward")
    CTcontext_folder_path = os.path.join(output_folderpath, bam_file_name + ".CTcontext_reads_features_forward_rev") # it will store also CT reads from reverse strand
    CCcontext_folder_path = os.path.join(output_folderpath, bam_file_name + ".CCcontext_reads_features_forward")

    # create output folders
    if os.path.exists(output_folderpath):
        shutil.rmtree(output_folderpath)
    os.mkdir(output_folderpath)

    if os.path.exists(TTcontext_folder_path):
        shutil.rmtree(TTcontext_folder_path)
    os.mkdir(TTcontext_folder_path)

    if os.path.exists(CTcontext_folder_path):
        shutil.rmtree(CTcontext_folder_path)
    os.mkdir(CTcontext_folder_path)

    if os.path.exists(CCcontext_folder_path):
        shutil.rmtree(CCcontext_folder_path)
    os.mkdir(CCcontext_folder_path)

    # start multiprocessing pool on contigs
    def consumer_worker(q):
        while True:
            contig = q.get()
            if contig != None:
                try:
                    T_and_C_reads_features_extractor_forward_rev(bam_file_path, ref_path, contig, TTcontext_folder_path, CTcontext_folder_path, CCcontext_folder_path, min_depth=min_depth)
                except pandas.errors.EmptyDataError:
                    print(f"[{datetime.now()}] [Global {contig} thread message] No well covered positions for this region.", flush=True)
                    continue
            else:
                # stopping the loop of the consumer
                print(f"[{datetime.now()}] [consumer worker message] Found end of Queue.\n", flush=True)
                break
                
    print(f"[{datetime.now()}] [Main thread message] Starting Iterations on Chromosomes (Forward Strand).\n", flush=True)

    # create a queue
    q = Queue()

    # create consumers processes
    consumers = []
    for t in range(threads_n):
        consumers.append( Process(target = consumer_worker, args=(q,)))
    print(f"[{datetime.now()}] [Main thread message] N° of Consumers: {len(consumers)}\n", flush=True)

    # start consumers
    for c in consumers:
        c.start()

    # fill the queue with the contigs names (chromosomes)
    for contig in contigs:
        q.put(contig)
    for t in range(threads_n):
        q.put(None) # put end signals to the queue for threads
    
    q.close()
    q.join_thread()
    
    # join consumers
    for c in consumers:
        c.join()
    
    # delete splitted bam files and related bai and depth files
    print(f"[{datetime.now()}] [Main thread message] Removing splitted BAM files and its related depth and bai files.", flush=True)
    for chrom in contigs:
        os.remove(f"{bam_file_base_path}.{chrom}.sorted.bam")
        os.remove(f"{bam_file_base_path}.{chrom}.sorted.bam.bai")
        os.remove(f"{bam_file_base_path}.{chrom}.sorted.bam.depth")
        
    
    stop_time = datetime.now()
    print(f"\n[{datetime.now()}] [Main thread message] Iteration on chromosomes finished. Elapsed time: {stop_time - start_time}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve TT and CT context T reads from a BAM file and reference fasta file and extract basecalling features for uForest algorithm.")
    parser.add_argument("-b", 
                        required=True, 
                        type=str, 
                        help="-bam_file_fullpath: \t a <str> with the full_path for the input bam file")
    parser.add_argument("-r", 
                        required=True, 
                        type=str, 
                        help="-reference_file_fullpath: \t a <str> with the full_path for the input reference fasta file used to align reads into the BAM file.")
    parser.add_argument("-o",
                        required=False,
                        default=None,
                        help="-output_folder_path: \t a <str> with the full_path for the output directory where results will be saved. By default output folder will be created into the same directory of the given bam file.")
    
    parser.add_argument("-threads", 
                        required=False,
                        type=int,
                        default=1,
                        help="-threads_number: \t a <int> indicating the number of threads to be used (it should be equal to the number of available CPUs).")

    parser.add_argument("-thr", 
                        required=False,
                        type=int,
                        default=50,
                        help="-threshold_value: \t a <int> indicating a threshold to filter out position with lower depth.")

    args = parser.parse_args()
    bam_file_path = args.b
    ref_path = args.r
    output_folder_path = args.o
    depth_threshold = args.thr
    threads_n = args.threads
    
    # launch main command function
    T_and_C_reads_features_extractor_forward_rev_main_functions(bam_file_path=bam_file_path, 
                                                                ref_path=ref_path, 
                                                                output_folderpath=output_folder_path, 
                                                                threads_n=threads_n, 
                                                                min_depth=depth_threshold)



    
