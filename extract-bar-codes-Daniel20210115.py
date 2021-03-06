#! /usr/bin/env python

import gzip
import re
import sys

# local classes
# project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
# sys.path.insert(0, project_dir)
# from InterleavedUMIReadDataDaniel20210115 import InterleavedUMIReadDataDaniel20210115

if len(sys.argv) != 2:
    sys.exit('A single command line argument specifying the fastq.gz file to process is required. Exiting.')
# first and only command line argument is the fastq.gz file to process
fastq_filename = sys.argv[1]

# cut-and-paste definitions from InterleavedUMIReadDataDaniel20210115 so this can be a simple single file script
def _extract_by_pos_and_length(seq, pos_length):
  return(''.join([seq[pos:pos+length] for pos, length in pos_length]))
  
class InterleavedUMIReadDataDaniel20210115:
  '''
  Information associated with a Read from an R1 file from with interleaved UMI and Well ID as used by 15 base version in experiment Daniel sent 20210115
  '''
  # constants
  umi_pos_length = [(0, 4), (9, 4)]
  well_id_pos_length = [(4, 5), (13, 2)]
  seq_length = sum([length for pos, length in umi_pos_length] + [length for pos, length in well_id_pos_length])
  
  def __init__(self, R1_seq):
    # Extract the IDs
    self.umi_well_seq = R1_seq[0:self.seq_length] # discard padding 'N's

  def __str__(self):
    string_rep = f'InterleavedUMIReadData: umi_well_seq: {self.umi_well_seq}'
    return string_rep
  
  @classmethod
  def extract_umi_and_well_id(cls, umi_well_seq):
    umi = _extract_by_pos_and_length(umi_well_seq, cls.umi_pos_length)
    well_id = _extract_by_pos_and_length(umi_well_seq, cls.well_id_pos_length)
    return(umi, well_id)
    
# create empty dictionary for counts
counts = {}

# Go through FASTQ, one four-line block at a time
ignore_Ns = True
n_skipped = 0
with gzip.open(fastq_filename, 'r') as fastq_file:
    while True: # Loop until we don't find another seq_id_line
        seq_id_line: str = fastq_file.readline().decode('ascii').rstrip()
        if not seq_id_line:
            break
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
        amplicon_id: str = amplicon_match.group(1)
        umi_well_seq: str = sequence[0:InterleavedUMIReadDataDaniel20210115.seq_length]
        if ignore_Ns and ('N' in umi_well_seq):
            n_skipped += 1
            continue
        (umi, well_id) = InterleavedUMIReadDataDaniel20210115.extract_umi_and_well_id(umi_well_seq)
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
