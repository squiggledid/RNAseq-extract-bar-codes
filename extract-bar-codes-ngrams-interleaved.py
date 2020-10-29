#! /usr/bin/env python

import gzip
import re
import sys
import os
import pickle
import colorama
colorama.init()
import numpy as np
import pandas as pd
import collections
import Levenshtein
from scipy.spatial.distance import pdist

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

import squtils as sq
from InterleavedUMIReadData import InterleavedUMIReadData
from ReadNgramHash import ReadNgramHash, seq_target_query

###############################################################################

mismatch_colours = (colorama.Fore.GREEN, colorama.Fore.YELLOW, colorama.Fore.MAGENTA, colorama.Fore.RED)
def highlight_well_id(umi_well_seq, position, mismatch):
  mismatch = min(mismatch, 3) # >2 should never happen as presently constructed
  return(f'{umi_well_seq[0:position]}{mismatch_colours[mismatch]}{umi_well_seq[position:position + InterleavedUMIReadData.well_id_length]}{colorama.Style.RESET_ALL}{umi_well_seq[position + InterleavedUMIReadData.well_id_length:]}')

###############################################################################
if (len(sys.argv)) < 2 or (len(sys.argv) > 3):
  sys.exit('A single command line argument specifying the fastq.gz file to process is required, followed by an optional argument specifying experimental metadata. Exiting.')
# first and only command line argument is the fastq.gz file to process
read_filename = sys.argv[1]
if len(sys.argv) == 3:
  metadata_filename = sys.argv[2]
else:
  metadata_filename = None

### Read metadata

if (metadata_filename != None) and os.path.exists(metadata_filename):
  well_id_colname = 'Well_Barcode_sequence'
  viral_copies_colname = 'Viral_copies'
  metadata = pd.read_csv(metadata_filename)
  well_ids = metadata[well_id_colname].to_list()
  well_id_counts = metadata[viral_copies_colname].to_list() # number of each well_id should be proportion to the number of viral copies in the well
  # well_ids = [(well_id) for count, well_id in sorted(zip(well_id_counts, well_ids), reverse = True)] # this should make finding well_ids faster, as we will look for the common ones first
  sorted_well_ids_counts = [(well_id, count) for count, well_id in sorted(zip(well_id_counts, well_ids), reverse = True)]
  known_well_id_counts_hash = dict(sorted_well_ids_counts)
  well_ids = [well_id for well_id, _ in sorted_well_ids_counts]
else:
  well_ids = []
  known_well_id_counts_hash = {}

### Read Interleaved reads
interleaved_read_ngrams = ReadNgramHash(InterleavedUMIReadData.seq_length) # create an object so we can access the ngram_length
interleaved_read_ngram_hash_filename = f'{read_filename}_ReadNgramHash_{interleaved_read_ngrams.ngram_length}_{InterleavedUMIReadData.umi_well_padding}.pkl'
if os.path.exists(interleaved_read_ngram_hash_filename):
  sq.log(f'Reading data from {interleaved_read_ngram_hash_filename}...')
  interleaved_read_ngram_hash_file = open(interleaved_read_ngram_hash_filename, 'rb')
  interleaved_read_ngrams = pickle.load(interleaved_read_ngram_hash_file)
  interleaved_read_ngram_hash_file.close()
else:
  sq.log(f'Building ReadNgramHash from {read_filename}...')
  # Go through FASTQ file, one four-line block at a time
  ignore_Ns = True
  n_read = 0
  n_skipped = 0
  report_every = 100000
  max_to_read = None # None for no limit :)
  # max_to_read = 1000 # For testing
  if max_to_read and (max_to_read < report_every):
    report_every = max_to_read
  print(f'InterleavedUMIReadData.seq_length: {InterleavedUMIReadData.seq_length}')
  interleaved_read_ngrams = ReadNgramHash(InterleavedUMIReadData.seq_length);
  print(f'interleaved_read_ngrams.seq_length: {interleaved_read_ngrams.seq_length}')
  with gzip.open(read_filename, 'r') as read_file:
    while (max_to_read == None) or (n_read < max_to_read): # Loop until we don't find another read_line, or have reached max_to_read
      read_line: str = read_file.readline().decode('ascii').rstrip()
      if not read_line:
        break
      # Extract R1 from ihe line
      R1_seq = read_line.split()[1]
      # create a InterleavedUMIReadData object
      interleaved_read = InterleavedUMIReadData(R1_seq)
      # # add to ngram hash
      interleaved_read_ngrams.insert(interleaved_read)
      # write a progress indicator
      n_read += 1
      if (n_read % report_every) == 0:
        print(f'%d items read from read_filename (%d skipped)' % (n_read, n_skipped))
  
  # save ReadNgramHash data structure with pickle
  sq.log(f'Saving ReadNgramHash to %s...' % interleaved_read_ngram_hash_filename)
  interleaved_read_ngram_hash_file = open(interleaved_read_ngram_hash_filename, 'wb')
  pickle.dump(interleaved_read_ngrams, interleaved_read_ngram_hash_file)
  interleaved_read_ngram_hash_file.close()

### Load or build well_id hash
max_well_id_offset = 4
max_dist = 2
interleaved_well_id_hash_filename = f'{read_filename}_InterleavedWellIDHash_{max_well_id_offset}_{max_dist}.pkl'
if os.path.exists(interleaved_well_id_hash_filename):
  sq.log(f'Reading data from {interleaved_well_id_hash_filename}...')
  interleaved_well_id_hash_file = open(interleaved_well_id_hash_filename, 'rb')
  interleaved_well_id_hash = pickle.load(interleaved_well_id_hash_file)
  interleaved_well_id_hash_file.close()
else:
  sq.log(f'Building interleaved_well_id_hash')
  # exhaustive looping is fine, as we need to visit every umi_well_seq, and the number of well_ids is low (i.e. 4)
  interleaved_well_id_hash = {}
  for umi_well_seq in interleaved_read_ngrams.umi_well_seq_hash:
    (umi, well_id) = InterleavedUMIReadData.extract_umi_and_well_id(umi_well_seq)
    interleaved_well_id_hash[umi_well_seq] = well_id
  # save well_id hash data structure
  sq.log(f'Saving interleaved_well_id_hash to %s...' % interleaved_well_id_hash_filename)
  interleaved_well_id_hash_file = open(interleaved_well_id_hash_filename, 'wb')
  pickle.dump(interleaved_well_id_hash, interleaved_well_id_hash_file)
  interleaved_well_id_hash_file.close()

### Summary
print(f'interleaved_read_ngrams.umi_well_seq_hash: {len(interleaved_read_ngrams.umi_well_seq_hash)} items')        
print(f'interleaved_well_id_hash: {len(interleaved_well_id_hash)} items')

well_ids = [well_id for well_id in interleaved_well_id_hash.values()]
well_id_counts_hash = collections.Counter(well_ids)
print(f'{len(well_id_counts_hash)} unique well_ids seen')
max_well_ids_to_print = 1000
well_ids_to_analyse = sorted(well_id_counts_hash, key = well_id_counts_hash.get, reverse = True)[0:max_well_ids_to_print]
for well_id in well_ids_to_analyse:
  if known_well_id_counts_hash:
    print(f'{well_id}: {well_id_counts_hash[well_id]} items (viral copies: {known_well_id_counts_hash[well_id]})')
  else:
    print(f'{well_id}: {well_id_counts_hash[well_id]} items')
well_id_array = np.array(well_ids_to_analyse).reshape(-1, 1)[1:len(well_ids_to_analyse), ] # we need a 2D array for pdist
well_id_distances = pdist(well_id_array, lambda well_id1, well_id2 : Levenshtein.hamming(well_id1[0], well_id2[0]))
well_id_distance_counts = collections.Counter(well_id_distances)
print(f'Average well_id Hamming distance:  {np.average(well_id_distances)}')
print(f'Maxzimum well_id Hamming distance: {np.max(well_id_distances)}')
print(f'Minimum well_id Hamming distance:  {np.min(well_id_distances)}')
print('Hamming distance distribution:')
print('\n'.join(f'\t{dist}: {well_id_distance_counts[dist]}' for dist in sorted(well_id_distance_counts)))
well_id_distances = pdist(well_id_array, lambda well_id1, well_id2 : Levenshtein.distance(well_id1[0], well_id2[0]))
well_id_distance_counts = collections.Counter(well_id_distances)
print(f'Average well_id Levenshtein distance:  {np.average(well_id_distances)}')
print(f'Maxzimum well_id Levenshtein distance: {np.max(well_id_distances)}')
print(f'Minimum well_id Levenshtein distance:  {np.min(well_id_distances)}')
print('Levenshtein distance distribution:')
print('\n'.join(f'\t{dist}: {well_id_distance_counts[dist]}' for dist in sorted(well_id_distance_counts)))

# quit()

### Tests

sq.log(f'Sorting umi_well_seqs')
sorted_umi_well_seqs = sorted(interleaved_read_ngrams.umi_well_seq_hash, key = lambda umi_well_seq : interleaved_read_ngrams.num_reads(umi_well_seq), reverse = True)#[0:10]
# find matchs for all umi_well_seqs
query_num = 1
for umi_well_seq in sorted_umi_well_seqs:
  sq.log(f'Querying with {umi_well_seq}...')
  ngram_matches = interleaved_read_ngrams.umi_well_seq_query(umi_well_seq)
  # # sort matches by total number of times each inexact match seen
  ngram_matches = sorted(ngram_matches, key = lambda ngram_match : interleaved_read_ngrams.num_reads(ngram_match[0]), reverse = True)
  num_inexact_matches = sum([interleaved_read_ngrams.num_reads(match_umi_well_seq) for match_umi_well_seq, num_ngram_matches in ngram_matches if match_umi_well_seq != umi_well_seq])
  print(f'{query_num}. Summary: {umi_well_seq}  Exact: {interleaved_read_ngrams.num_reads(umi_well_seq)} Inexact: {num_inexact_matches}')
  padding_len = len(f'{query_num}. Summary: ')
  for match in ngram_matches:
    (umi, well_id) = InterleavedUMIReadData.extract_umi_and_well_id(match[0])
    print(f'{match[0]} (UMI: {umi} Well: {well_id}) {interleaved_read_ngrams.num_reads(match[0]):5} (sim: {match[1]})')
  query_num += 1


