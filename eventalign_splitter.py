# ---- Importing basic modules --------#
import os, sys, shutil
import pandas as pd
import numpy as np
from pandas import DataFrame
from datetime import datetime
from multiprocessing import Process, Queue
import argparse
from datetime import datetime
from __init__ import __version__

# ---- defining utils functions -------#

# def. aggregation function for samples
def agg_samples(x):
    '''
    Function to aggregate samples feature in the groupby operation
    merging all samples from the same event mapping onto a position
    into a unique list of float values 
    '''
    final = list( map(float, ",".join(x).split(",")))
    
    return final


# defining custom function for preprocessing of filtered eventalign dataframes.
def collapse_per_read_splitted_eventalign_table(df):
    '''
    (ONLY RNA! FOR DNA THE START POSITION SHOULD BE 5!)
    Function to preprocess and collapse events mapping on the same position and kmer 
    like to nanoRMS (Begik et al., 2021).
    It takes in input:
        - df --> dataframe pandas of filtered eventalign output.
    Returns:
        - df_grouped --> dataframe pandas of grouped and collapsed reads using samples values to calculate statistics.
    '''
    # filtering events with NNNNN model_kmer
    df = df[df["model_kmer"] != "NNNNN"]
    #Collapsing multiple observations from the same read on position and ref. kmer
    df_grouped = df.groupby(["contig", "position", "reference_kmer", "read_name"]).agg({"samples":[agg_samples]})
    # resetting index and columns' names
    df_grouped = df_grouped.reset_index()
    df_grouped.columns =  df_grouped.columns.droplevel(-1)
    # create dwell, weighted average of event_level_mean and event_level_std colums from elements in the merged samples list
    df_grouped["event_level_mean"] = df_grouped["samples"].apply(np.mean)
    #df_grouped["event_level_std"] = df_grouped["samples"].apply(np.std)
    df_grouped["dwell"] = df_grouped["samples"].apply(len).values
    df_grouped.drop("samples", axis=1, inplace=True)
    df_grouped.sort_values(["contig", "read_name"], inplace=True)
    df_grouped.reset_index(drop=True, inplace=True)
    
    return df_grouped


# defining custom function for 2nd rounds preprocessing of filtered eventalign dataframes with more than one resulting collapsed dataframe.
def collapse_per_read_splitted_eventalign_table_2ndRound(df):
    '''
    (ONLY RNA! FOR DNA THE START POSITION SHOULD BE 5!)
    Function to preprocess and collapse events mapping on the same position and kmer 
    like to nanoRMS but used id the df of the current reads already exists.
    It takes in input:
        - df --> dataframe pandas of filtered collapsed eventalign output with more than one round.
    Returns:
        - df_grouped --> dataframe pandas of grouped and collapsed reads using samples values to calculate statistics.
    '''
    # filtering events with NNNNN model_kmer
    df = df[df[2] != "NNNNN"]
    # create new features to implement the average mean of means and stds. Then another for the sum of the dwell times.
    df["mean_*_dwell"] = df[4] * df[5]
    #Collapsing multiple observations from the same read on position and ref. kmer
    df_grouped = df.groupby([0, 1, 2, 3]).agg({"mean_*_dwell":["sum"], 5:["sum"]})
    # resetting index and columns' names
    df_grouped = df_grouped.reset_index()
    df_grouped.columns =  df_grouped.columns.droplevel(-1)
    # create dwell, weighted average of event_level_mean and event_level_std colums from elements in the merged samples list
    df_grouped[4] = df_grouped["mean_*_dwell"] / df_grouped[5]
    df_grouped.drop(["mean_*_dwell"], axis=1, inplace=True)
    df_grouped = df_grouped.filter([0,1,2,3,4,5])
    df_grouped.sort_values(by=[0, 1, 2])
    
    return df_grouped


# ---- defining producers and consumers workers -------------------------------#
def producer(eventalign_path, out_dir, q, threads_n):
    print(f"[{datetime.now()}] [Producer Message] Starting Iteration on the input eventalign file...\n", flush=True)
    start = datetime.now()

    reads_idx_file = open( os.path.join(out_dir, "0.collapsed_reads.idx"), "wt" )

    with open(eventalign_path, "rt") as f:
        read_counter = 0
        last_read_id = None
        for line_count, line in enumerate(f):
            if not line_count == 0:
                line_list = line.rstrip().split("\t")
                current_read_id = line_list[0] + " - " + line_list[3]
                if last_read_id != current_read_id:
                    if read_counter % 10000 == 0:
                        print(f"[{datetime.now()}] [Producer Message] Elapsed Time {datetime.now() - start}", flush=True)
                        print(f"[{datetime.now()}] [Producer Message] Total Reads currently processed by producer:", read_counter, flush=True)
                        print(f"[[{datetime.now()}] [Producer Message] Total Eventalign file's lines currently processed by producer:", line_count ,"\n" , flush=True)
                    last_read_id = current_read_id
                    # only the first read will not have an opened file of the previous read. Close only on reads with count > 1
                    try:
                        q.put({"read_list_id":out_file_path, "read_list":read_list})
                    except:
                        pass
                    read_counter += 1
                    contig = line_list[0]
                    # assign a unique output file name with the counter of read number detected, contig name assigned at the first event of the read and its read-name 
                    out_file_path = os.path.join(out_dir, f"{contig}.{line_list[3]}.eventalign.collapsed")
                    if not os.path.exists(out_file_path):
                        # write a line on the reads index file relative to the new current read only if the reads doesn't exists yet.
                        reads_idx_file.write(f"{contig}\t{line_list[3]}\t{out_file_path}\n")
                    # create a new list where to store every event of that read.
                    read_list = []
                    # write the first event of the current new line into the new output file
                    read_list.append(line_list)
                else:
                    # write the current event/line into the output file of the current read.
                    read_list.append(line_list)
        # for last read list, adding it to queue
        q.put({"read_list_id":out_file_path, "read_list":read_list})
    
    # put ending elements in queue for stopping signal to consumers
    for t in range(threads_n):
        q.put({"read_list_id":None, "read_list":None})

    # close reads_idx_file.
    reads_idx_file.close()
    
    # drop duplicates if exist.
    read_index_path = os.path.join(out_dir, "0.collapsed_reads.idx")
    read_idx = pd.read_table(read_index_path, header=None)
    read_idx.drop_duplicates(inplace=True)
    read_idx.to_csv(read_index_path, sep="\t", header=None, index=None)

    # print producer's statistics.
    stop = datetime.now()
    print(f"[{datetime.now()}] [Producer Message] Producer finished to split the eventalign file.", flush=True)
    print(f"[{datetime.now()}] [Producer Message] Producer Elapsed time:", stop - start, flush=True)
    print(f"[{datetime.now()}] [Producer Message] Producer Total processed lines: {line_count + 1}", flush=True)
    print(f"[{datetime.now()}] [Producer Message] Producer Total processed reads: {read_counter}", flush=True)
    
    q.close()
    q.join_thread()

    print(f"[{datetime.now()}] [Producer Message] Finished", flush=True)


def consumer_worker(q, id_consumer):
    columns = "contig position reference_kmer read_name strand event_index event_level_mean event_stdv event_length model_kmer model_mean model_stdv standardized_level start_idx end_idx samples".split(" ")
    while True:
        eventalign_read_dict = q.get()
        eventalign_read_path = eventalign_read_dict["read_list_id"]
        eventalign_read_list = eventalign_read_dict["read_list"]
        if eventalign_read_path != None: # END SIGNAL AT THE END OF THE QUEUE
            df_input = pd.DataFrame(eventalign_read_list)
            df_input.columns = columns
            output_path = eventalign_read_path
            df_processed = collapse_per_read_splitted_eventalign_table(df_input)
            df_processed.to_csv(output_path, sep="\t", header=None, index=None, mode='a') # if exist append to it and then preprocess 2nd round again to recollapse events.
            if os.path.exists(output_path):
                df_input_2 = pd.read_table(output_path, header=None)
                df_processed_2 = collapse_per_read_splitted_eventalign_table_2ndRound(df_input_2)
                df_processed_2.to_csv(output_path, sep="\t", header=None, index=None)
        else:
            # Stopping the loop of the consumer if found a end-signal (None) in the Queue.
            print(f"[{datetime.now()}] [Consumer {id_consumer} Message] Found end of Queue", flush=True)
            break


def eventalign_splitter(eventalign_path, threads_n=1):
    '''
    Function that split a unique eventalign file produced by f5c or nanopolish eventalign into a folder containing
    a file per read detected into the same contig. Each file contains a table of the read with collapsed currents statistics 
    per site-kmer (event_level_mean and dwell time). Each file has a semantic filename containing read_counter, contig and 
    read-name for a successive easier retrieve of the data. 
    In parallel, an index file will be produced in the same output folder, containing for each read a row with the corresponding 
    contig and its path for the collapsed file into the folder.
    The program leverages on Queue and Process classes of the Python multiprocessing package. 
    One producer will parse, line by line the eventalign file and will put into a Queue the events of the same reads mapping on a contig.
    Then a pool of consumer processes, will read from that Queue and will process each bunch of read-related events writing it as a tsv table
    on disk with its corresponding semantic name. 
    Please launch the script from the same folder where you want to to be produced the eventalign.collapsed folder containing all the processed read-tables.
    
    Inputs:
        -eventalign_path --> a string with the full_path for the eventalign file to be splitted.
        -threads_n --------> an <int> with the number of consumers to be used to dequeue. (default 1)
    
    To Note: f5c or nanopolish eventalign MUST be lauched with flags --samples --scale-events --print-read-names --signal-index.

    The output folder with collapsed tables is produced in the same folder where the eventalign table is present.
    '''
    start_global = datetime.now()
    threads_n = threads_n
    print(f"[{datetime.now()}] [Main Process Message] Thread used: {threads_n}", flush=True)
    eventalign_path = eventalign_path
    out_dir = f"{eventalign_path}.collapsed"
    print(f"[{datetime.now()}] [Main Process Message] Input eventalign file: {eventalign_path}", flush=True)
    print(f"[{datetime.now()}] [Main Process Message] Output folder for eventalign collapsed reads: {out_dir}", flush=True)

    # creating a new folder if it doesn't exist and if it exist eliminate that and create a new directory
    if os.path.exists(out_dir):
        shutil.rmtree(out_dir)
    os.makedirs(out_dir)
    
    # create a queue where producer will put the path of each splitted read file and the consumer will read them afterwards
    q = Queue(maxsize=2000)

    # create consumers processes
    consumers = []
    for t in range(threads_n):
        consumers.append( Process(target = consumer_worker, args=(q,t+1,)))
    print(f"[{datetime.now()}] [Main Process Message] N° of Consumers: {len(consumers)}", flush=True)
    # start consumers
    for c in consumers:
        c.start()

    # create a producer process and start it
    p = Process(target=producer, args=(eventalign_path, out_dir, q, threads_n))
    p.start()

    # join consumers
    for c in consumers:
        c.join()
    # join producer
    p.join()
    
    stop_global = datetime.now()

    print(f"[{datetime.now()}] Global Elapsed time: {stop_global - start_global}", flush=True)
    print(f"[{datetime.now()}] EXITING...Queue size is:", q.qsize(), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    Program name: eventalign_splitter (multiprocessing version)\n\n
    
    Function that split a unique eventalign file produced by f5c or nanopolish eventalign into a folder containing
    a file per read detected into the same contig. Each file contains a table of the read with collapsed currents statistics 
    per site-kmer (event_level_mean, std and dwell time). Each file has a semantic filename containing read_counter, contig and 
    read-name for a successive easier retrieve of the data.\n
    In parallel, an index file will be produced in the same output folder, containing for each read a row with the corresponding 
    contig and its path for the collapsed file into the folder. \n
    The program leverage on Queue and Process classes of the Python multiprocessing package.\n
    One producer will parse, line by line the eventalign file and will put into a Queue the events of the same reads mapping on a contig.
    Then a pool of consumer processes, will read from that Queue and will process each bunch of read-related events writing it as a tsv table
    on disk with its corresponding semantic name.\n
    Please launch the script from the same folder where you want to to be produced the eventalign.collapsed folder containing all the processed read-tables.\n\n

    To Note: f5c or nanopolish eventalign MUST be lauched with flags --samples --scale-events --print-read-names --signal-index.

    The output folder with collapsed tables is produced in the same folder where the eventalign table is present.
    ''')

    parser.add_argument('-e', required=True, type=str, help='-eventalign_path: \t a <str> with the full_path for the eventalign file to be splitted.')
    parser.add_argument('-t', default=1, type=int, help='-threads_n: \t an <int> with the number of consumers to be used to dequeue.')

    args = parser.parse_args()

    eventalign_path = args.e
    threads_n = args.t

    # print some starting info related to the C_to_U_classifier version, to the used program and to the input arguments
    print(f"[{datetime.now()}] C_to_U_classifier version: {__version__}", flush=True)
    print(f"[{datetime.now()}] Eventalign Splitter. Input arguments:", flush=True)
    for argument in args.__dict__.keys():
        print(f"\t- {argument} --> {args.__dict__[argument]}", flush=True)

    # launching main function of program if executed as console script.
    eventalign_splitter(eventalign_path, threads_n)