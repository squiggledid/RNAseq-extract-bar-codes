#! /usr/bin/env python

import gzip
import re
import sys
import pickle

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from UMIHashTrie import UMIHashTrie
from UMIData import UMIData

if len(sys.argv) != 2:
    sys.exit('A single command line argument specifying the fastq.gz file to process is required. Exiting.')
# first and only command line argument is the fastq.gz file to process
fastq_filename = sys.argv[1]

# constants
umi_start = 0
umi_length = 16
well_id_start = 16
well_id_length = 8

report_every = 100000

# create empty UMIHashTrie to store and search all id count data
umi_trie = UMIHashTrie()

# Go through FASTQ, one four-line block at a time
n_read = 0
with gzip.open(fastq_filename, 'r') as fastq_file:
    while True: # Loop until we don't find another seq_id_line
        seq_id_line: str = fastq_file.readline().decode('ascii').rstrip()
        if not seq_id_line:
            break
        # Check the four lines
        if seq_id_line[0] != '@':
            sys.exit("Expected first line of a read, beginning with @. '" + seq_id_line[0] + "' seen. Exiting.")
        sequence: str = fastq_file.readline().decode('ascii').rstrip()
        plus_line: str = fastq_file.readline().decode('ascii').rstrip()
        if plus_line[0] != '+':
            sys.exit("Expected third line of a read, beginning with '+'. '" + seq_id_line[0] + "' seen. Exiting.")
        quality: str = fastq_file.readline().decode('ascii').rstrip()
        amplicon_match = re.search(':([^:]+)$', seq_id_line)
        if amplicon_match == None:
            sys.exit('No amplicon ID found at end of sequence ID line:\n' + seq_id_line + '\nExiting.')
        # Extract the IDs
        amplicon_id: str = amplicon_match.group(1)
        umi = sequence[umi_start:(umi_start + umi_length - 1)]
        well_id = sequence[well_id_start:(well_id_start + well_id_length - 1)]
        # add to UMIHashTrie
        umi_trie.record_read(umi, well_id, amplicon_id)
        # write a progress indicator
        n_read += 1
        if (n_read % report_every) == 0:
            print(f'%d reads added to the UMIHashTrie' % n_read)

# save UMITrie data structure with pickle
umi_trie_filename = fastq_filename + '_UMIHashTrie.pkl'
print(f'Saving UMIHashTrie to %s...' % umi_trie_filename)
umi_trie_file = open(umi_trie_filename, 'wb')
pickle.dump(umi_trie, umi_trie_file)
umi_trie_file.close()

