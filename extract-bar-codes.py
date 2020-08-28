#! /usr/bin/env python

import gzip
import re
import sys

if len(sys.argv) != 2:
    sys.exit('A single command line argument specifying the fastq.gz file to process is required. Exiting.')
# first and only command line argument is the fastq.gz file to process
fastq_filename = sys.argv[1]

# create empty nested dictionary for counts
counts = {}

# Go through FASTQ, one four-line block at a time
with gzip.open(fastq_filename, 'r') as fastq_file:
    while True: # Loop until we don't find another seq_id_line
        seq_id_line: str = fastq_file.readline().decode('ascii').rstrip()
        if not seq_id_line:
            break
        if seq_id_line[0] != '@':
            sys.exit("Expected first line of a read, beginning with @. '" + seq_id_line[0] + "' seen. Exiting.")
        sequence = fastq_file.readline().decode('ascii').rstrip()
        plus_line: str = fastq_file.readline().decode('ascii').rstrip()
        if plus_line[0] != '+':
            sys.exit("Expected third line of a read, beginning with '+'. '" + seq_id_line[0] + "' seen. Exiting.")
        quality: str = fastq_file.readline().decode('ascii').rstrip()
        amplicon_match = re.search(':([^:]+)$', seq_id_line)
        if amplicon_match == None:
            sys.exit('No amplicon ID found at end of sequence ID line:\n' + seq_id_line + '\nExiting.')
        amplicon_id: str = amplicon_match.group(1)
        umi = sequence[0:15]
        well_id = sequence[16:23]
        counts_key = umi + ':' + well_id + ':' + amplicon_id
        if counts_key in counts:
            counts[counts_key] += 1
        else:
            counts[counts_key] = 1

# Print result as a csv file
# TODO use a proper csv writer for this
print('UMI,Well_ID,Amplicon_ID,Count')
sorted_keys = sorted(counts, key = counts.get, reverse = True)
for key in sorted_keys:
    ids_match = re.search('^([^:]+):([^:]+):([^:]+)$', key)
    if ids_match == None:
        sys.exit('counts key does not match expected format. Exiting.')
    print(ids_match.group(1) + ',' + ids_match.group(2) + ',' + ids_match.group(3) + ',' + str(counts[key]))
