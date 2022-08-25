import os, subprocess, shutil
from typing_extensions import Required
import pandas as pd
from datetime import date, datetime
import argparse
from multiprocessing import Queue, Process, Value

'''
Script to retrieve currents and dwells time from eventalign collapsed table produced by Eventalign_splitter.py script.
The focus will be done on sites with a depth > min_depth and only reads mapping on T or C positions on reference.
The output will be splitted into a folder within the same directory of the input bam file (or if specified into a specific path) with three files:
    
    - output_basefolder_path/T_C_reference_currents_dwells_features_reads/T_reference_reads_currents_dwells_reads.tsv (all -3/+3 regions of reads mapping around a well covered T on reference)
    - output_basefolder_path/T_C_reference_currents_dwells_features_reads/C_reference_reads_currents_dwells_reads.tsv (all -3/+3 regions of reads mapping around a well covered C on reference)
    - output_basefolder_path/T_C_reference_currents_dwells_features_reads/C_reference_reads_putative_reverse_currents_dwells_reads.tsv (all -3/+3 regions of reads mapping around a well 
            covered position with a G on the most-right base of related event 5-mer, as putative reverse strand read with a mapping C or a T: it should be filtered with the supporting informations
            of the basecalling features script).

File are tsv tables with the following structure, where the pos_1_based filed is referred to the position 0 of the context with respect 
the reference (cur0=mean current at position 0; dw0=dwell at position 0):

# COLUMNS IN OUTPUT FILES
contig  pos_1_based read_name   cur-3  cur-2    cur-1   cur0    cur+1   cur+2   cur+3   dw-3    dw-2    dw-1    dw0 dw+1    dw+2    dw+3

To note: Output files must be filtered using files produced by TT_CT_CC_basecalling_features_retriever.py in order to retain only
rows mapping on forward strand (and reverse/minus strand for CT context reads) and which have a given context (e.g. TT, CC, and CT). 
'''


### --- main function --- ###
def currents_dwells_retriever(bam_file_path, eventalign_collapsed_folder_path, min_depth=50, output_basefolder_path=None, continue_iterations=False, threads_n=1):
    '''
    Retrieve currents and dwells from -3/+3 regions of reads mapping the central position to a T or C on the reference.
    Continue iteration allow to continue an interrupted iteration from the last reads into collapsed_reads index file.
    '''
    sample_name = os.path.basename(bam_file_path).split(".")[0]

    # Define outputs path and files
    if output_basefolder_path == None:
        output_folder_path = os.path.join(os.path.dirname(bam_file_path), sample_name + ".T_C_reference_currents_dwells_features_reads")
    else:
        output_folder_path = os.path.join(output_basefolder_path, sample_name + ".T_C_reference_currents_dwells_features_reads")
    
    T_reference_reads_filepath = os.path.join(output_folder_path, "T_reference_reads_currents_dwells_reads.tsv")
    C_reference_reads_filepath = os.path.join(output_folder_path, "C_reference_reads_currents_dwells_reads.tsv")
    C_reference_reads_putative_reverse_filepath = os.path.join(output_folder_path, "C_reference_reads_putative_reverse_currents_dwells_reads.tsv")

    start_time = datetime.now()

    if continue_iterations == False:
        if os.path.exists(output_folder_path):
            shutil.rmtree(output_folder_path)
        os.mkdir(output_folder_path)
        log_file = open(os.path.join(output_folder_path, "log_file.log"), "w")
    elif continue_iterations == True:
        log_file = open(os.path.join(output_folder_path, "log_file.log"), "a")

    # produce depth file
    if not os.path.exists(f"{bam_file_path + '.depth'}"):
        print(f"[{datetime.now()}] Retrieving depths for current bam file: {bam_file_path}.", flush=True, file=log_file)
        subprocess.run(f"samtools depth {bam_file_path}".split(" "), stdout=open(f"{bam_file_path + '.depth'}", "w"))

    # retrieve list of well covered positions
    pos_well_covered = pd.read_table(f"{bam_file_path + '.depth'}", header=None)
    pos_well_covered.columns=["contig", "pos_1_based", "depth"]
    pos_well_covered = pos_well_covered[pos_well_covered["depth"] >= min_depth]
    pos_well_covered.reset_index(inplace=True, drop=True)
    pos_well_covered_coords_0_based = [i[1]+":"+str(i[2]-1) for i in pos_well_covered.itertuples()]
    del(pos_well_covered)
    pos_well_covered_coords_0_based = set(pos_well_covered_coords_0_based)
    print(f"[{datetime.now()}] [Main Thread Message] A total of {len(pos_well_covered_coords_0_based)} sites with a depth higher than {min_depth} have been found for {bam_file_path}", flush=True, file=log_file)

    # remove depth file
    os.remove(f"{bam_file_path + '.depth'}")

    # open reads collapsed index file
    print(f"[{datetime.now()}] [Main Thread Message] Opening the eventalign collapsed reads index file from {eventalign_collapsed_folder_path} directory.", flush=True, file=log_file)
    eventalign_reads_index = pd.read_table(os.path.join(eventalign_collapsed_folder_path, "0.collapsed_reads.idx"), header=None)
    eventalign_reads_index.columns = ["contig", "read_name", "eventalign_collapsed_path"]
    # if we want to continue previous interrupted iteration
    if continue_iterations == True:
        last_read_filepath = os.path.join(os.path.dirname(T_reference_reads_filepath), "last_read")
        os.system(f"cat {T_reference_reads_filepath} | tail | cut -f 3 | sort | uniq > {last_read_filepath}")
        with open(last_read_filepath, "r") as f:
            last_read = f.read().rstrip()
        os.remove(last_read_filepath)
        # discard reads until the last reads into index file
        eventalign_reads_index = eventalign_reads_index.iloc[eventalign_reads_index[eventalign_reads_index["read_name"]==last_read].index[0]:,:]

    if continue_iterations == False:
        print(f"[{datetime.now()}] [Main Thread Message] A total of {eventalign_reads_index.shape[0]} reads with collapsed eventalign events have been found into index file.", flush=True, file=log_file)
    elif continue_iterations == True:
        print(f"[{datetime.now()}] [Main Thread Message] A total of {eventalign_reads_index.shape[0]} reads with collapsed eventalign events have been found into index file after the last read {last_read}.", flush=True, file=log_file)

    # go to the right position
    os.chdir(eventalign_collapsed_folder_path)
    os.chdir("..")


    def currents_dwells_retriever_function(collapsed_read, counter, failed_reads_counter, dim_chunk, consumer_n, T_reference_reads_file, C_reference_reads_file, C_reference_reads_putative_reverse_file, pos_well_covered_coords_0_based):
        contig = collapsed_read[1]
        read_path = collapsed_read[3]
        # open eventalign collapsed related to current read
        try:
            read_eventalign_collapsed = pd.read_table(read_path, header=None)
            read_eventalign_collapsed.columns = ["contig", "pos_0_based", "ref_kmer", "read_name", "mean_current", "dwell_time"]
            read_eventalign_collapsed = read_eventalign_collapsed.query(f"contig == '{contig}'")
            if not read_eventalign_collapsed.empty:
                coord_min = read_eventalign_collapsed.pos_0_based.min()
                coord_max = read_eventalign_collapsed.pos_0_based.max()
                region_df = pd.DataFrame([i for i in range(coord_min, coord_max+1)], columns=["pos_0_based"])
                region_df["contig"] = [read_eventalign_collapsed.contig.value_counts().index.values[0] for i in range(region_df.shape[0])]
                read_eventalign_collapsed = pd.merge(read_eventalign_collapsed, region_df, how="right", on=["pos_0_based", "contig"])
                for pos in read_eventalign_collapsed.iloc[3:-3,:].itertuples():
                    coord = pos[1]+":"+str(pos[2])
                    if coord in pos_well_covered_coords_0_based:
                        if not pd.isna(pos.ref_kmer):
                            if pos.ref_kmer[0] == "T":
                                filtered = read_eventalign_collapsed[(read_eventalign_collapsed["pos_0_based"] >= pos.pos_0_based-3 )&(read_eventalign_collapsed["pos_0_based"] <= pos.pos_0_based+3)]
                                # create an output and write to file of T contexts (pos-1-based)
                                output = [pos.contig, str(pos.pos_0_based+1), pos.read_name] + [str(round(cur,3)) for cur in filtered.mean_current.to_list()] + [str(round(dwell,0)) for dwell in filtered.dwell_time]
                                output = "\t".join(output) + "\n"
                                T_reference_reads_file.write(output)
                            elif pos.ref_kmer[0] == "C":
                                filtered = read_eventalign_collapsed[(read_eventalign_collapsed["pos_0_based"] >= pos.pos_0_based-3 )&(read_eventalign_collapsed["pos_0_based"] <= pos.pos_0_based+3)]
                                # create an output and write to file of C contexts (pos-1-based)
                                output = [pos.contig, str(pos.pos_0_based+1), pos.read_name] + [str(round(cur,3)) for cur in filtered.mean_current.to_list()] + [str(round(dwell,0)) for dwell in filtered.dwell_time]
                                output = "\t".join(output) + "\n"
                                C_reference_reads_file.write(output)
                            if pos.ref_kmer[-1] == "G":
                                filtered = read_eventalign_collapsed[(read_eventalign_collapsed["pos_0_based"] >= pos.pos_0_based-3 )&(read_eventalign_collapsed["pos_0_based"] <= pos.pos_0_based+3)]
                                # create an output and write to file of putative reverse G contexts contexts (pos-1-based) but reversed since it could be on reverse strand
                                # next filtering step will assess that question
                                # to note: position will be implemented by +4 and converted to 1-based since the reference kmer is the reverse complement of the real one if 
                                # the read maps on reverse (-) strand.
                                output = [pos.contig, str(pos.pos_0_based+1+4), pos.read_name] + [str(round(cur,3)) for cur in filtered.mean_current.to_list()][::-1] + [str(round(dwell,0)) for dwell in filtered.dwell_time][::-1]
                                output = "\t".join(output) + "\n"
                                C_reference_reads_putative_reverse_file.write(output)

        except FileNotFoundError:
            with failed_reads_counter.get_lock():
                failed_reads_counter.value += 1
        with counter.get_lock():
            counter.value += 1
        if counter.value % dim_chunk == 0:
            # update progress of whole process
            print(f"[{datetime.now()}] [consumer worker {consumer_n} message] {q.qsize()} remaining items to process..", flush=True, file=log_file)
    
    def consumer_worker(q, counter, failed_reads_counter, dim_chunk, consumer_n, T_reference_reads_filepath, C_reference_reads_filepath, C_reference_reads_putative_reverse_filepath, pos_well_covered_coords_0_based, continue_iterations):
        if continue_iterations == False:
            mode = "w"
        elif continue_iterations == True:
            mode = "a"            
        T_reference_reads_file_consumer = open(T_reference_reads_filepath+f".tmp{consumer_n}", mode)
        C_reference_reads_file_consumer = open(C_reference_reads_filepath+f".tmp{consumer_n}", mode)
        C_reference_reads_putative_reverse_file_consumer = open(C_reference_reads_putative_reverse_filepath+f".tmp{consumer_n}", mode)
        while True:
            try:
                collapsed_read = q.get()
                if collapsed_read != None:
                    currents_dwells_retriever_function(collapsed_read, counter, failed_reads_counter, dim_chunk, consumer_n, T_reference_reads_file_consumer, C_reference_reads_file_consumer, C_reference_reads_putative_reverse_file_consumer, pos_well_covered_coords_0_based)
                else:
                    # stopping the loop of the consumer
                    print(f"[{datetime.now()}] [consumer worker {consumer_n} message] Found end of Queue.\n", flush=True, file=log_file)
                    T_reference_reads_file_consumer.close()
                    C_reference_reads_file_consumer.close()
                    C_reference_reads_putative_reverse_file_consumer.close()
                    break
            except Exception as error:
                print(f"[{datetime.now()}] [consumer worker {consumer_n} message] ERROR {error}", flush=True, file=log_file)
                break

    counter = Value("i", 0)
    failed_reads_counter = Value("i", 0)
    n_chunks = 100
    dim_chunk = int(eventalign_reads_index.shape[0] / n_chunks)

    # define threads number
    print(f"[{datetime.now()}] [Main Thread Message] Working with {threads_n} threads in parallel.", flush=True, file=log_file)

    # create a queue
    q = Queue()

    # create consumers processes
    consumers = []
    for t in range(threads_n):
        consumer_n = t+1
        consumers.append( Process(target = consumer_worker, args=(q, counter, failed_reads_counter, dim_chunk, consumer_n, T_reference_reads_filepath, C_reference_reads_filepath, C_reference_reads_putative_reverse_filepath, pos_well_covered_coords_0_based, continue_iterations,)))

    # start consumers
    for c in consumers:
        c.start()

    # fill the queue with the names
    for collapsed_read in eventalign_reads_index.itertuples(name=None):
        q.put(collapsed_read)
    
    for t in range(threads_n):
        q.put(None) # put end signals to the queue for threads
    
    q.close()
    q.join_thread()

    # join consumers
    for c in consumers:
        c.join()

    # merging outputs produced by each thread and compressing with gzip
    print(f"[{datetime.now()}] [Main Thread Message] Merging output files.", flush=True, file=log_file)
    os.system(f"cat {T_reference_reads_filepath}.* > {T_reference_reads_filepath}")
    os.system(f"rm {T_reference_reads_filepath}.*")
    os.system(f"gzip {T_reference_reads_filepath}")
    os.system(f"cat {C_reference_reads_filepath}.* > {C_reference_reads_filepath}")
    os.system(f"rm {C_reference_reads_filepath}.*")
    os.system(f"gzip {C_reference_reads_filepath}")
    os.system(f"cat {C_reference_reads_putative_reverse_filepath}.* > {C_reference_reads_putative_reverse_filepath}")
    os.system(f"rm {C_reference_reads_putative_reverse_filepath}.*")
    os.system(f"gzip {C_reference_reads_putative_reverse_filepath}")

    # TROVARE UN MODO PER FARE IL SORTING SULLA BASE DELLE PRIME 3 COLONNE (opzionale)

    # print some final logs and statistics
    print(f"[{datetime.now()}] [Main Thread Message] Computations finished.", flush=True, file=log_file)
    print(f"[{datetime.now()}] [Main Thread Message] Total reads evaluated from index file: {counter.value}", flush=True, file=log_file)
    print(f"[{datetime.now()}] [Main Thread Message] Total reads failed to be retrieved: {failed_reads_counter.value}", flush=True, file=log_file)
    print(f"[{datetime.now()}] [Main Thread Message] All output files were saved inside the folder: {output_folder_path}", flush=True, file=log_file)
    print(f"[{datetime.now()}] [Main Thread Message] Elapsed time: {datetime.now() - start_time}", flush=True, file=log_file)
    log_file.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Retrieve currents and dwells time from eventalign collapsed reads produced by Eventalign_splitter.py script.")
    
    parser.add_argument("-b",
                        required=True,
                        type=str,
                        help="-bam_file_path: \t a <str> indicating the full path of the bam file to retrieve well covered sites wiht a depth higher than 'min_depth'")
    
    parser.add_argument("-efp",
                        required=True,
                        type=str,
                        help="-eventalign_collapsed_folder_path: \t a <str> indicating the folder path of the directory in which collapsed eventalign reads are stored together with the 0.collapsed_reads.idx index file.")
    
    parser.add_argument("-d",
                        required=False,
                        type=int,
                        default=50,
                        help="-min_depth: \t a <int> indicating the theshold value to be used to retrieve well covered position from bam file. This will be used to filter out position in eventalign collapsed files which are not well covered in order to spare time.")

    parser.add_argument("-o",
                        required=False,
                        type=str,
                        default=None,
                        help="-output_basefolder_path: \t a <str> indicating the base folder path where to save the otput folder named as 'output_basefolder_path/T_C_reference_currents_dwells_features_reads'.")

    parser.add_argument("-threads", 
                        required=False,
                        type=int,
                        default=1,
                        help="-threads_number: \t a <int> indicating the number of threads to be used (it should be equal to the number of available CPUs).")

    parser.add_argument("-ci",
                        required=False,
                        type=bool,
                        default=False,
                        help="-continue_iterations: \t a <bool> indicating if one would to continue iterations after an early stop, starting from the last read processed.")

    args = parser.parse_args()
    bam_file_path = args.b
    eventalign_collapsed_folder_path = args.efp
    min_depth = args.d
    output_basefolder_path = args.o
    continue_iterations = args.ci
    threads_n = args.threads

    currents_dwells_retriever(bam_file_path=bam_file_path, 
                              eventalign_collapsed_folder_path=eventalign_collapsed_folder_path, 
                              min_depth=min_depth, 
                              output_basefolder_path=output_basefolder_path, 
                              continue_iterations=continue_iterations, 
                              threads_n=threads_n)