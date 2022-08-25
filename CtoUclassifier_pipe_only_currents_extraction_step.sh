#!/bin/bash
source /lustrehome/afonzino/anaconda3/bin/activate nanocompore
source /lustrehome/afonzino/C_to_U_classifier_plus_min/venv/bin/activate

python3 /lustrehome/afonzino/C_to_U_classifier_plus_min/pipe_only_currents_extraction_step.py \
    --bam_filepath $1 \
    --eventalign_collapsed $2 \
    -threads $3