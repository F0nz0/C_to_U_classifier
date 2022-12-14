# CtoUclassifier

A Machine and Deep Learning based Python package for the Nanopore noise attenuation, useful for the amelioration of C-to-U editing signal in direct-RNA sequencing experiments.

**Required Softwares**:\
C_to_U_classifiers uses internally (and not) some software that should be installed preferably into a new conda enviroment. \
After the activation of the conda enviroment install the following softwares:
1) Python >= 3.7
2) f5c >= 0.7
3) Samtools >= 1.3.1
4) Minimap2 == 2.24
5) Guppy == 5.0.11 (please see documentation on https://community.nanoporetech.com/downloads)

**Installation**:
1) Download the source code from GitHub repository at the url:
        
    https://github.com/F0nz0/C_to_U_classifier

2) Enter the C_to_U_classifier folder and decompress the iForest_pretrained_model folder:
		
		tar -xf iForest_pretrained_model.tar.gz

3) Create a new virtual environment (it's suggested to create, use and activate a base conda environment with all the required software):

		# create a new conda environment
        conda create --name CtoU python=3.8

		# activate the conda env
		conda activate CtoU

		# install samtools
		conda install -c bioconda samtools

		# install minimap2
		conda install -c bioconda minimap2

		# install f5c (optimised re-implementation of the call-methylation and eventalign modules in Nanopolish)
		conda install -c bioconda f5c

		# create virtual environment inside Conda CtoU env
		python3 -m venv venv

4) Activate the venv:
	
	    source venv/bin/activate

5) Upgrade pip version:
	
	    python3 -m pip install --upgrade pip

6) Install wheel package via pip:
	
	    pip install wheel

7) Install required Python packages using the requirements.txt file:

        python -m pip install -r requirements.txt

**Basic Usage**:\
The C_to_U_classifier is able to ameliorate the C-to-U editing signal in direct-RNA Nanopore runs. 
It does this by mean of two main pipelines:

A)  *Basecalling Pipeline*:\
	It extracts from the BAM file C-to-U signals at per-read level for every position with a depth higher than a given threshold (default 50 reads).\
	Then it uses a pretrained iForest anomaly detection algorithm to correct the C-to-U signal ameliorating the signal-to-noise ratio using the extracted basecalling features (central base quality, mean interval quality and discarding reads intervals with indels and/or mismatches).

B) *Currents Pipeline*:\
	This pipeline includes the basecalling one since it needs basecalling features to work. 
	Firts of all, it extracts the currents intensities and dweel times retreived by the use of the f5c eventalign software and then, groups each position-read pair as CC, CT and TT contexts basing on Guppy data extracted from the BAM file.
	Then, it performs a CNN-WaveNet model training on the CC and TT context reads and tries to correct the C-to-U signal classifying CT context reads as CC or TT (error-T or real-T, respectively) basing on currents features.

The pipelines give back as output the results at both per-read level and genome-space aggregated level.\
If you are interested in testing the packages on a test dataset please download at the url (and unzip the archive):

https://drive.google.com/file/d/1oayc3gttcoB7bHwuoTAshb93lWrXVCG9/view?usp=sharing

The first step is to perform the basecalling of the raw data contained into Fast5 files using the Guppy basecalling software version 5.0.11 with the following command (use -x flag if cuda is already installed on your system, counterwise you can avoid it):
	
    guppy_basecaller -c rna_r9.4.1_70bps_hac.cfg -i [FAST5_DIR] -s [OUTPUT_DIR] -r --fast5_out -x "cuda:0"

The reference genome that should be used can be downloaded from GENECODE website and should be filtered from non-chromosome regions. The software search for chr* format to detect a region-contig. Thus, if needed, modify in advance accordingly reference fasta file headers of each chromosome in order to be: >chr1, >chr2, >chr3 and so on.

Then, fastQ files generated by Guppy software that have passed the quality control step (inside the pass forlder, see Guppy documentation for more detail), have to be mapped against the reference genome (tested on Human hg38 and Murine mm39 genomes) using minimap2 software using the following command:
	
	minimap2 -t [THREADS] -ax splice -uf -k14 --secondary=no [REFERENCE_FILEPATH] [FASTQ_FILEPATH] > [SAM_FILEPATH]

Now, the next step is to perform some filtering and the conversion of SAM to BAM file using samtools, then index it and remove the SAM file:
	
	# filtering of unmapped reads sorting and conversion to binary alignment file
	samtools view -b -F 4 [SAM_FILEPATH] | samtools sort -O BAM > [BAM_FILEPATH]
	
	# index the BAM file
	samtools index [BAM_FILEPATH]
	
	# remove unused SAM file
	rm [SAM_FILEPATH]

If your are interested in currents pipeline, there is the need for the remapping of raw currents against the reference genome using f5c software. To do that, first of all you need to index fastQ and fast5 files as following:
	
	f5c index -d [FAST5_DIR] [FASTQ_FILEPATH]

Now it is possible to launch f5c eventalign program with the following parameters (it's of pivotal importance to set all these flags):
	
	f5c eventalign --iop [IOP] -t [THREADS] -r [FASTQ_FILEPATH] \
		-b [BAM_FILEPATH] \
		-g [REFERENCE_FILEPATH] \
		--rna \
		--scale-events \
		--print-read-names \
		--samples \
		--signal-index > [EVENTALIGN_FILEPATH]

Since the output of the f5c eventalign command is huge (TB scale also for minion devices), repetitive and redundant, the C_to_classifier package has a Python script useful for collapsing events mapping on the same genome position decreasing thus the amount of data to be analyzed. Furthemore, this step allows for a more flexible indexing strategy of currents data at per-read level. Launch thus, the following command (to note, C_to_U_classifier scripts need for full-path of the input files in order to work finely):

	python eventalign_splitter.py -e [EVENTALIGN_FILEPATH]-t [THREADS]

Now the script should have produced a folder named as *[EVENTALIGN_FILEPATH].collapsed* containing an index file *0.collapsed_reads.idx* and a file for each contig-readname pair with currents intensity values and dwell times for every retrieved genomic coordinate.
At this point, it is possible to choose either the basecalling or currents pipeline.

**Usage of basecalling pipeline**:

	python pipe_basecalling.py \
		-B [BAM_FILEPATH] \
		-R [REFERENCE_FILEPATH] \
		-threads [THREADS] \
		-O [ORGANISM: human OR murine] # <-- used for the computation of the APOBEC1 signature with a Likelihood-ratio test

The script will produce a folder with the extracted basecalling features and a folder containing the resulting iForest model correction data at both per-read and genomic-space levels.

**Usage of currents pipeline**:

	python pipe_currents_extraction_step.py \
		-B [BAM_FILEPATH] \
		-R [REFERENCE_FILEPATH] \
		-EC [EVENTALIGN_COLLAPSED_FOLDER] \
		-threads [THREADS] \
		-O [ORGANISM: human OR murine] # <-- used for the computation of the APOBEC1 signature with a Likelihood-ratio test

After this, in addition to the basecalling pipeline related outputs, a directory named *[BAM_FILE_ROOT_NAME].T_C_reference_currents_dwells_features_reads* containg currents and dwells times features will be produced.
The final step is to train the CNN-WaveNet model and to make prediction on CT context reads for the C-to-U signal correction:

	python train_model_from_scratch.py \
		-mn [MODEL_NAME_SUFFIX] \
		-ffp [CURRENTS_DWELLS_FEATURES_DIRECTORY] \
		-bam [BAM_FILEPATH] \
		-ref_filepath [REFERENCE_FILEPATH] \
		-O [ORGANISM: human OR murine] # <-- used for the computation of the APOBEC1 signature with a Likelihood-ratio test

Also in this case, after the training and prediction steps, a folder named *[BAM_FILE_ROOT_NAME].model_CNN_[BAM_FILE_ROOT_NAME]_[MODEL_NAME_SUFFIX]* will be produced as output containing both per-read and genome-space aggregated data for the CNN-WaveNet model.