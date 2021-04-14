#!/usr/bin/env bash

# Obtain ID of 5UTRs with note: Potential AUG annotation error
grep 'annotation error' data/original_data/Nagalakshmi_2008_5UTRs_V64.gff3 | cut -f9 | awk '{split($0,a,";"); print a[1]}' | awk '{split($0,b,"="); print b[2]}' > data/derived_data/utr5_aug_annotation_error.txt

# Retain only 5UTRs with no conflicts with annotated CDS
grep -v -F -f data/derived_data/utr5_aug_annotation_error.txt data/original_data/Nagalakshmi_2008_5UTRs_V64.bed > data/derived_data/Nagalakshmi_2008_5UTRs_V64_wo_conflict_UTRs.bed

# Eliminate YLL054C from 3'UTRs annotation since it is located inside CDS
grep -v 'YLL054C' data/original_data/Nagalakshmi_2008_3UTRs_V64.bed > data/derived_data/Nagalakshmi_2008_3UTRs_V64_wo_conflict_UTRs.bed

# Process UTRs annotations
python3 scr/process_UTR_annotations.py

# Union of the annotations ensuring no duplicated entries
cat data/derived_data/Nagalakshmi_2008_3UTRs_V64_formatted.bed \
data/derived_data/Nagalakshmi_2008_5UTRs_V64_formatted.bed \
| sort | uniq > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_w_UTRs.bed

# Sort the annotation
sortBed -i data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_w_UTRs.bed > data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_w_UTRs_sorted.bed

# Obtain FASTA sequences from the annotation and write it to a TAB file (reverse complement if feature in antisense)
bedtools getfasta -name -tab -s -fi data/original_data/Saccharomyces_cerevisiae.R64-1-1_release-89.dna.chromosome.I.fa \
-bed data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_w_UTRs_sorted.bed \
-fo data/derived_data/Saccharomyces_cerevisiae.R64-1-1.89_w_UTRs_sequences.tab
