#! /usr/bin/env python

import gzip
import sys
import re

# flags for development
debug = True

# Data location - will eventually replace with command line argument
fastq_dir = "/data/projects/WEHI_CoViD_RNAseq/DB CloudStor 20200812"
fastq_filename = fastq_dir + "/" + "Multiplex_Virus_30_R1.fastq.gz"

with gzip.open(fastq_filename, 'r') as fastq_file:
    seq_id_line = fastq_file.readline().decode('ascii').rstrip()
    if seq_id_line[0] != '@':
        sys.exit("Expected first line of a read, beginning with @. '" + seq_id_line[0] + "' seen. Exiting.")
    sequence = fastq_file.readline().decode('ascii').rstrip()
    plus_line = fastq_file.readline().decode('ascii').rstrip()
    if plus_line[0] != '+':
        sys.exit("Expected third line of a read, beginning with '+'. '" + seq_id_line[0] + "' seen. Exiting.")
    quality = fastq_file.readline().decode('ascii').rstrip()
    if debug:
        print(seq_id_line, sequence, quality, sep = "\n")
    amplicon_match = re.search(":([^:]+)$", seq_id_line)
    if amplicon_match == None:
        sys.exit("No amplicon ID found at end of sequence ID line:\n" + seq_id_line + "\nExiting.")
    amplicon_id = amplicon_match.group(1)
    umi = sequence[0:15]
    well_id = sequence[16:23]
    if debug:
        print('amplicon ID:\t' + amplicon_id, 'UMI:\t\t' + umi, 'Well ID:\t' + well_id, sep = '\n')