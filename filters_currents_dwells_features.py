'''
A script to filter currents and dwells extracted by currents_dwell_retriever.py that are not
on TT, CC or CT context. It uses list of sites-reads from the basecalling features retrieved 
by TT_CT_CC_basecalling_features_retriever.py
'''

import os, sys
from datetime import datetime
from tqdm import tqdm
import argparse
import gzip


### --- utils functions --- ###

def filter_TTcontext_currents_dwells_from_Treferece_reads_forward(TTcontext_forward_reads_folder_path, T_C_reference_currents_dwells_features_reads):
    '''
    Function to filter among Treference currents and dwells reads 
    those that have a TT context basing on TT_context basecalling reads mapping on forward strand 
    table. It saves the filtered version into the same folder of T_C_reference_currents_dwells_features_reads.
    '''
    
    TTcontext_forward_reads_merged_filepath = os.path.join(TTcontext_forward_reads_folder_path, "TTcontext_reads_features_forward_chr_whole.tsv")
    T_reference_reads_currents_dwells_filepath = os.path.join(T_C_reference_currents_dwells_features_reads, "T_reference_reads_currents_dwells_reads.tsv.gz")
    TTcontext_reads_currents_dwells_output_filepath = os.path.join(T_C_reference_currents_dwells_features_reads, "TTcontext_reads_currents_dwells_reads_forward.tsv") # output file with filtered reads

    if not os.path.exists(TTcontext_forward_reads_merged_filepath):
        print(f"[{datetime.now()}] Merging tsv files into a unique file.", flush=True)
        os.system(f"cat {TTcontext_forward_reads_folder_path}/TTcontext_reads_features_forward_chr* > {TTcontext_forward_reads_folder_path}/TTcontext_reads_features_forward_chr_whole.tsv")
        print(f"[{datetime.now()}] Merging Finished.", flush=True)
    
    TTcontext_reads_merged_file = open(TTcontext_forward_reads_merged_filepath, "r")
    T_reference_reads_currents_dwells_file = gzip.open(T_reference_reads_currents_dwells_filepath, "r") 
    TTcontext_reads_currents_dwells_output_file = open(TTcontext_reads_currents_dwells_output_filepath, "w")

    # set where to load target sites to be filtered
    TTcontext_reads_sites = set()

    print(f"[{datetime.now()}] Producing a set containing a list of reads which should be used to filter currents_dwells of T reference.", flush=True)
    
    # calculate number of rows into TTcontext table
    TTcontext_rows = 0
    for r in TTcontext_reads_merged_file:
        TTcontext_rows += 1
    TTcontext_reads_merged_file.close()
    chunks = 10000
    n_reads_per_chunk = int(TTcontext_rows / chunks)
    
    # reopen TTcontext merged file
    TTcontext_reads_merged_file = open(TTcontext_forward_reads_merged_filepath, "r")
    with tqdm(total=TTcontext_rows, file=sys.stdout) as pbar:
        counter = 0
        for r in TTcontext_reads_merged_file:
            read = r.rstrip().split("\t")
            # retrieve strand of the file
            if counter == 0:
                strand = read[3]
            read_id = read[0]+":"+read[1]+":"+read[2]
            TTcontext_reads_sites.add(read_id)
            counter += 1
            if counter % n_reads_per_chunk == 0:
                pbar.update(n_reads_per_chunk)
                pbar.refresh()

    print(f"[{datetime.now()}] A total of {len(TTcontext_reads_sites)} unique sites have been retrieved from {counter} reads from TTcontext_basecalling_features tsv table.", flush=True)

    # after the set was retrieved the tsv can be closed
    TTcontext_reads_merged_file.close()
    # remove merged basecalling_features file (whole file)
    os.remove(TTcontext_forward_reads_merged_filepath)

    # filter currents and dwells reads from T reference using the retrieve list of sites
    print(f"[{datetime.now()}] Filtering Treference currents and dwells...", flush=True)
    Treference_rows = 0
    for r in T_reference_reads_currents_dwells_file:
        Treference_rows += 1
    T_reference_reads_currents_dwells_file.close()
    chunks = 10000
    n_reads_per_chunk = int(Treference_rows / chunks)

    # reopen T reference currents and dwells files
    print(f"[{datetime.now()}] Filtering currents and dwells of reads on T reference by the use of sites list of basecalling features.", flush=True)
    T_reference_reads_currents_dwells_file = gzip.open(T_reference_reads_currents_dwells_filepath, "r") 
    with tqdm(total=Treference_rows, file=sys.stdout) as pbar:
        counter = 0
        counter_retained = 0
        for r in T_reference_reads_currents_dwells_file:
            read = r.decode("utf-8").rstrip().split("\t")
            read_id = read[0]+":"+read[1]+":"+read[2]
            counter += 1
            if counter % n_reads_per_chunk == 0:
                pbar.update(n_reads_per_chunk)

            if read_id in TTcontext_reads_sites:
                output_read = read[:3]+[strand]+read[3:]
                output_read = "\t".join(output_read) + "\n"
                TTcontext_reads_currents_dwells_output_file.write(output_read)
                counter_retained += 1

    print(f"[{datetime.now()}] On a total amount of {counter} Treference_reads_currents_dwells, a total of {counter_retained} have been retained and saved on {TTcontext_reads_currents_dwells_output_filepath}", flush=True)

    # close all tsv files
    T_reference_reads_currents_dwells_file.close()
    TTcontext_reads_currents_dwells_output_file.close()


def filter_CC_and_CT_context_currents_dwells_from_Treferece_reads_forward(CCcontext_forward_reads_folder_path, CTcontext_forward_reads_folder_path, T_C_reference_currents_dwells_features_reads):
    '''
    Function to filter among Creference currents and dwells reads 
    those that have a CC or CT context basing on CC and CT context basecalling reads mapping on forward strand (and reverse only for CT contexts) 
    tables. It saves the filtered version into the same folder of T_C_reference_currents_dwells_features_reads.
    '''
    CCcontext_forward_reads_merged_filepath = os.path.join(CCcontext_forward_reads_folder_path, "CCcontext_reads_features_forward_chr_whole.tsv")
    CTcontext_forward_rev_reads_merged_filepath = os.path.join(CTcontext_forward_reads_folder_path, "CTcontext_reads_features_forward_rev_chr_whole.tsv")

    C_reference_reads_currents_dwells_filepath = os.path.join(T_C_reference_currents_dwells_features_reads, "C_reference_reads_currents_dwells_reads.tsv.gz")
    C_reference_reads_putative_reverse_currents_dwells_filepath  =os.path.join(T_C_reference_currents_dwells_features_reads, "C_reference_reads_putative_reverse_currents_dwells_reads.tsv.gz")
    
    CCcontext_reads_currents_dwells_output_filepath = os.path.join(T_C_reference_currents_dwells_features_reads, "CCcontext_reads_currents_dwells_reads_forward.tsv") # output file with filtered reads for CC contexts
    CTcontext_reads_currents_dwells_output_filepath = os.path.join(T_C_reference_currents_dwells_features_reads, "CTcontext_reads_currents_dwells_reads_forward_rev.tsv") # output file with filtered reads for CT contexts

    # merging CC context basecalling features files
    if not os.path.exists(CCcontext_forward_reads_merged_filepath):
        print(f"[{datetime.now()}] Merging CC_context tsv files into a unique file.", flush=True)
        os.system(f"cat {CCcontext_forward_reads_folder_path}/CCcontext_reads_features_forward_chr* > {CCcontext_forward_reads_folder_path}/CCcontext_reads_features_forward_chr_whole.tsv")
        print(f"[{datetime.now()}] Merging Finished.", flush=True)

    # merging CT context basecalling features files
    if not os.path.exists(CTcontext_forward_rev_reads_merged_filepath):
        print(f"[{datetime.now()}] Merging CT_context tsv files into a unique file.", flush=True)
        os.system(f"cat {CTcontext_forward_reads_folder_path}/CTcontext_reads_features_forward_rev_chr* > {CTcontext_forward_reads_folder_path}/CTcontext_reads_features_forward_rev_chr_whole.tsv")
        print(f"[{datetime.now()}] Merging Finished.", flush=True)
    

    CCcontext_reads_merged_file = open(CCcontext_forward_reads_merged_filepath, "r")
    CTcontext_reads_merged_file = open(CTcontext_forward_rev_reads_merged_filepath, "r")
    
    C_reference_reads_currents_dwells_file = gzip.open(C_reference_reads_currents_dwells_filepath, "r")
    C_reference_reads_putative_reverse_currents_dwells_file = gzip.open(C_reference_reads_putative_reverse_currents_dwells_filepath, "r")
    
    CCcontext_reads_currents_dwells_output_file = open(CCcontext_reads_currents_dwells_output_filepath, "w")
    CTcontext_reads_currents_dwells_output_file = open(CTcontext_reads_currents_dwells_output_filepath, "w")

    # set where to load target sites to be filtered
    CCcontext_reads_sites = set()
    CTcontext_reads_sites_forw = set()
    CTcontext_reads_sites_rev = set()

    print(f"[{datetime.now()}] Producing three sets containing a list of reads (CC and CT contexts, the latter with forw and rev info) which should be used to filter currents_dwells of C reference and putative reverse C.", flush=True)
    
    # calculate number of rows into CCcontext table
    CCcontext_rows = 0
    for r in CCcontext_reads_merged_file:
        CCcontext_rows += 1
    CCcontext_reads_merged_file.close()
    chunks = 10000
    n_reads_per_chunk = int(CCcontext_rows / chunks)
    
    # reopen CCcontext merged file
    CCcontext_reads_merged_file = open(CCcontext_forward_reads_merged_filepath, "r")
    with tqdm(total=CCcontext_rows, file=sys.stdout) as pbar:
        counter = 0
        for r in CCcontext_reads_merged_file:
            read = r.rstrip().split("\t")
            # retrieve strand of the file
            if counter == 0:
                strand = read[3]
            read_id = read[0]+":"+read[1]+":"+read[2]
            CCcontext_reads_sites.add(read_id)
            counter += 1
            if counter % n_reads_per_chunk == 0:
                pbar.update(n_reads_per_chunk)

    print(f"[{datetime.now()}] A total of {len(CCcontext_reads_sites)} unique sites have been retrieved from {counter} reads from CCcontext_basecalling_features tsv table.", flush=True)

    # after the set was retrieved the tsv can be closed
    CCcontext_reads_merged_file.close()
    # remove merged basecalling_features file (whole file)
    os.remove(CCcontext_forward_reads_merged_filepath)


    # calculate number of rows into CTcontext table
    CTcontext_rows = 0
    for r in CTcontext_reads_merged_file:
        CTcontext_rows += 1
    CTcontext_reads_merged_file.close()
    if CTcontext_rows > 10000:
        chunks = 10000
    else:
        chunks = CTcontext_rows
    n_reads_per_chunk = int(CTcontext_rows / chunks)
    
    # reopen CTcontext merged file
    CTcontext_reads_merged_file = open(CTcontext_forward_rev_reads_merged_filepath, "r")
    with tqdm(total=CTcontext_rows, file=sys.stdout) as pbar:
        counter = 0
        counter_forw = 0
        counter_rev = 0
        for r in CTcontext_reads_merged_file:
            read = r.rstrip().split("\t")
            # retrieve strand of the file
            strand = read[3]
            read_id = read[0]+":"+read[1]+":"+read[2]
            if strand == "+":
                CTcontext_reads_sites_forw.add(read_id)
                counter_forw += 1
            elif strand == "-":
                CTcontext_reads_sites_rev.add(read_id)
                counter_rev += 1
            counter += 1
            if counter % n_reads_per_chunk == 0:
                pbar.update(n_reads_per_chunk)

    print(f"[{datetime.now()}] A total of {len(CTcontext_reads_sites_forw)} unique sites mapping on forward (+) strand have been retrieved from {counter} reads from CTcontext_basecalling_features tsv table.", flush=True)
    print(f"[{datetime.now()}] A total of {len(CTcontext_reads_sites_rev)} unique sites mapping on reverse (-) strand have been retrieved from {counter} reads from CTcontext_basecalling_features tsv table.", flush=True)

    # after the set was retrieved the tsv can be closed
    CTcontext_reads_merged_file.close()
    # remove merged basecalling_features file (whole file)
    os.remove(CTcontext_forward_rev_reads_merged_filepath)


    # filter currents and dwells reads from C reference using the retrieve list of sites from forward (+) strand
    print(f"[{datetime.now()}] Filtering Creference currents and dwells mapping on forward (+) strand...", flush=True)
    
    # calculating the total amount of rows into Creference file of forward (+) reads
    Creference_rows = 0
    for r in C_reference_reads_currents_dwells_file:
        Creference_rows += 1
    C_reference_reads_currents_dwells_file.close()
    chunks = 10000
    n_reads_per_chunk = int(Creference_rows / chunks)

    # reopen C reference currents and dwells files
    print(f"[{datetime.now()}] Filtering currents and dwells of forward reads on C reference by the use of the two sites lists of basecalling features retrieved from CC and CT context tables.", flush=True)
    C_reference_reads_currents_dwells_file = gzip.open(C_reference_reads_currents_dwells_filepath, "r") 
    with tqdm(total=Creference_rows, file=sys.stdout) as pbar:
        counter = 0
        counter_retained_CC = 0
        counter_retained_CT = 0
        strand  = "+"
        for r in C_reference_reads_currents_dwells_file:
            read = r.decode("utf-8").rstrip().split("\t")
            read_id = read[0]+":"+read[1]+":"+read[2]
            counter += 1
            if counter % n_reads_per_chunk == 0:
                pbar.update(n_reads_per_chunk)
            if read_id in CCcontext_reads_sites:
                output_read = read[:3]+[strand]+read[3:]
                output_read = "\t".join(output_read) + "\n"
                CCcontext_reads_currents_dwells_output_file.write(output_read)
                counter_retained_CC += 1
            elif read_id in CTcontext_reads_sites_forw:
                output_read = read[:3]+[strand]+read[3:]
                output_read = "\t".join(output_read) + "\n"
                CTcontext_reads_currents_dwells_output_file.write(output_read)
                counter_retained_CT += 1

    print(f"[{datetime.now()}] On a total amount of {counter} Creference_reads_currents_dwells, a total of {counter_retained_CC} CCcontext reads of forward (+) strand have been retained and saved on {CCcontext_reads_currents_dwells_output_filepath}", flush=True)
    print(f"[{datetime.now()}] On a total amount of {counter} Creference_reads_currents_dwells, a total of {counter_retained_CT} CTcontext reads of forward (+) strand have been retained and saved on {CTcontext_reads_currents_dwells_output_filepath}", flush=True)

    # close C_reference_reads_currents_dwells_file tsv
    C_reference_reads_currents_dwells_file.close()

    # filter currents and dwells reads from putative C reference using the retrieve list of sites from reverse (-) strand
    print(f"[{datetime.now()}] Filtering putative Creference currents and dwells mapping on reverse (-) strand...", flush=True)
    
    # calculating the total amount of rows into Creference file of forward (+) reads
    Creference_rows = 0
    for r in C_reference_reads_putative_reverse_currents_dwells_file:
        Creference_rows += 1
    C_reference_reads_putative_reverse_currents_dwells_file.close()
    chunks = 10000
    n_reads_per_chunk = int(Creference_rows / chunks)

    # reopen putative reverse (-) C reference currents and dwells files
    print(f"[{datetime.now()}] Filtering currents and dwells of putative reverse (-) reads on C reference by the use of the sites list of basecalling features retrieved from CT reverse (-) context tables.", flush=True)
    C_reference_reads_putative_reverse_currents_dwells_file = gzip.open(C_reference_reads_putative_reverse_currents_dwells_filepath, "r")
    with tqdm(total=Creference_rows, file=sys.stdout) as pbar:
        counter = 0
        counter_retained_CT = 0
        strand  = "-"
        for r in C_reference_reads_putative_reverse_currents_dwells_file:
            read = r.decode("utf-8").rstrip().split("\t")
            read_id = read[0]+":"+read[1]+":"+read[2]
            counter += 1
            if counter % n_reads_per_chunk == 0:
                pbar.update(n_reads_per_chunk)
            if read_id in CTcontext_reads_sites_rev:
                output_read = read[:3]+[strand]+read[3:]
                output_read = "\t".join(output_read) + "\n"
                CTcontext_reads_currents_dwells_output_file.write(output_read)
                counter_retained_CT += 1
    
    print(f"[{datetime.now()}] On a total amount of {counter} putative reverse (-) strand Creference_reads_currents_dwells, a total of {counter_retained_CT} CTcontext reads of reverse (-) strand have been retained and saved on {CTcontext_reads_currents_dwells_output_filepath}", flush=True)

    CCcontext_reads_currents_dwells_output_file.close()
    CTcontext_reads_currents_dwells_output_file.close()


### --- main function --- ###
def filters_currents_dwells_features(TTcontext_forward_reads_folder_path, 
                                     CCcontext_forward_reads_folder_path, 
                                     CTcontext_forward_reads_folder_path, 
                                     T_C_reference_currents_dwells_features_reads):
    
    start_time =datetime.now()

    print(f"[{datetime.now()}] Starting filtering T reference currents and dwells reads which have TT context.", flush=True)
    filter_TTcontext_currents_dwells_from_Treferece_reads_forward(TTcontext_forward_reads_folder_path, T_C_reference_currents_dwells_features_reads)

    print(f"[{datetime.now()}] Starting filtering C reference currents and dwells reads which have CC or CT context.", flush=True)
    filter_CC_and_CT_context_currents_dwells_from_Treferece_reads_forward(CCcontext_forward_reads_folder_path, CTcontext_forward_reads_folder_path, T_C_reference_currents_dwells_features_reads)
    
    print(f"[{datetime.now()}] Iteration finished.", flush=True)
    print(f"[{datetime.now()}] Elapsed time: {datetime.now() - start_time}", flush=True)




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="A script to filter currents and dwells extracted by currents_dwell_retriever.py that are not on TT, CC or CT context. It uses list of sites-reads from the basecalling features retrieved by TT_CT_CC_basecalling_features_retriever.py")
    
    parser.add_argument("-TT",
                        required=True,
                        type=str,
                        help="-TTcontext_forward_reads_folder_path: \t a <str> indicating the folder path where are stored the basecalling features of reads extracted from TT contexts.")

    parser.add_argument("-CC",
                        required=True,
                        type=str,
                        help="-CCcontext_forward_reads_folder_path: \t a <str> indicating the folder path where are stored the basecalling features of reads extracted from CC contexts.")

    parser.add_argument("-CT",
                        required=True,
                        type=str,
                        help="-CTcontext_forward_reads_folder_path: \t a <str> indicating the folder path where are stored the basecalling features of reads extracted from CT contexts.")

    parser.add_argument("-TCref",
                        required=True,
                        type=str,
                        help="T_C_reference_currents_dwells_features_reads_folder_path: \t a <str> indicating the folder path where are stores the currents and dwells features extracted from eventalign collapsed reads mapping on T or C reference bases.")

    args = parser.parse_args()
    TTcontext_forward_reads_folder_path = args.TT
    CCcontext_forward_reads_folder_path = args.CC
    CTcontext_forward_reads_folder_path = args.CT
    T_C_reference_currents_dwells_features_reads = args.TCref

    # launch main function
    filters_currents_dwells_features(TTcontext_forward_reads_folder_path=TTcontext_forward_reads_folder_path,
                                     CCcontext_forward_reads_folder_path=CCcontext_forward_reads_folder_path,
                                     CTcontext_forward_reads_folder_path=CTcontext_forward_reads_folder_path,
                                     T_C_reference_currents_dwells_features_reads=T_C_reference_currents_dwells_features_reads)
