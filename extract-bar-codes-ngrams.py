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
from FastqReadData import FastqReadData
from FastqReadNgramHash import FastqReadNgramHash, seq_target_query

###############################################################################

mismatch_colours = (colorama.Fore.GREEN, colorama.Fore.YELLOW, colorama.Fore.MAGENTA, colorama.Fore.RED)
def highlight_well_id(umi_well_seq, position, mismatch):
  mismatch = min(mismatch, 3) # >2 should never happen as presently constructed
  return(f'{umi_well_seq[0:position]}{mismatch_colours[mismatch]}{umi_well_seq[position:position + FastqReadData.well_id_length]}{colorama.Style.RESET_ALL}{umi_well_seq[position + FastqReadData.well_id_length:]}')

###############################################################################
if (len(sys.argv)) < 2 or (len(sys.argv) > 3):
  sys.exit('A single command line argument specifying the fastq.gz file to process is required, followed by an optional argument specifying experimental metadata. Exiting.')
# first and only command line argument is the fastq.gz file to process
fastq_filename = sys.argv[1]
if len(sys.argv) == 3:
  metadata_filename = sys.argv[2]
else:
  metadata_filename = None

### Read metadata

if os.path.exists(metadata_filename):
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

### Read FASTQ reads

fastq_read_ngram_hash_filename = fastq_filename + '_FastqReadNgramHash_' + str(FastqReadNgramHash.ngram_length) + '.pkl'
if os.path.exists(fastq_read_ngram_hash_filename):
  sq.log(f'Reading data from {fastq_read_ngram_hash_filename}...')
  fastq_read_ngram_hash_file = open(fastq_read_ngram_hash_filename, 'rb')
  fastq_read_ngrams = pickle.load(fastq_read_ngram_hash_file)
  fastq_read_ngram_hash_file.close()
else:
  sq.log(f'Building FastqReadNgramHash from {fastq_read_ngram_hash_filename}...')
  # Go through FASTQ file, one four-line block at a time
  ignore_Ns = True
  n_read = 0
  n_skipped = 0
  report_every = 100000
  max_to_read = None # None for no limit :)
  # max_to_read = 1000 # For testing
  if max_to_read and (max_to_read < report_every):
    report_every = max_to_read
  umi_well_seq_end = FastqReadData.umi_start + FastqReadData.umi_length + FastqReadData.well_id_length + FastqReadData.umi_well_padding - 1
  fastq_read_ngrams = FastqReadNgramHash();
  with gzip.open(fastq_filename, 'r') as fastq_file:
    while (max_to_read == None) or (n_read < max_to_read): # Loop until we don't find another read_id_line, or have reached max_to_read
      read_id_line: str = fastq_file.readline().decode('ascii').rstrip()
      if not read_id_line:
        break
      # Check the four lines
      if read_id_line[0] != '@':
        sys.exit("Expected first line of a read, beginning with @. '" + read_id_line[0] + "' seen. Exiting.")
      sequence: str = fastq_file.readline().decode('ascii').rstrip()
      plus_line: str = fastq_file.readline().decode('ascii').rstrip()
      if plus_line[0] != '+':
        sys.exit("Expected third line of a read, beginning with '+'. '" + read_id_line[0] + "' seen. Exiting.")
      quality: str = fastq_file.readline().decode('ascii').rstrip()
      # bail out now if we're ignoring Ns in the IDs, and there is one
      if ignore_Ns and ('N' in sequence[FastqReadData.umi_start:umi_well_seq_end]):
        n_skipped += 1
        continue
      # create a FastqReadData object
      fastq_read = FastqReadData(read_id_line, sequence, quality)
      # # add to ngram hash
      fastq_read_ngrams.insert(fastq_read)
      # write a progress indicator
      n_read += 1
      if (n_read % report_every) == 0:
        print(f'%d items read from fastq_filename (%d skipped)' % (n_read, n_skipped))
  
  # save FastqReadNgramHash data structure with pickle
  sq.log(f'Saving FastqReadNgramHash to %s...' % fastq_read_ngram_hash_filename)
  fastq_read_ngram_hash_file = open(fastq_read_ngram_hash_filename, 'wb')
  pickle.dump(fastq_read_ngrams, fastq_read_ngram_hash_file)
  fastq_read_ngram_hash_file.close()

### Load or build well_id hash
max_well_id_offset = 4
max_dist = 4
fastq_well_id_hash_filename = f'{fastq_filename}_FastqWellIDHash_{max_well_id_offset}_{max_dist}.pkl'
if os.path.exists(fastq_well_id_hash_filename):
  sq.log(f'Reading data from {fastq_well_id_hash_filename}...')
  fastq_well_id_hash_file = open(fastq_well_id_hash_filename, 'rb')
  fastq_well_id_hash = pickle.load(fastq_well_id_hash_file)
  fastq_well_id_hash_file.close()
else:
  sq.log(f'Building fastq_well_id_hash')
  # exhaustive looping is fine, as we need to visit every umi_well_seq, and the number of well_ids is low (i.e. 4)
  fastq_well_id_hash = {}
  for umi_well_seq in fastq_read_ngrams.umi_well_seq_hash:
    best_well_id = None
    best_pos = None
    best_dist = FastqReadData.well_id_length
    for well_id in well_ids:
      (this_pos, this_dist) = seq_target_query(well_id, umi_well_seq, FastqReadData.well_id_start, max_well_id_offset)
      if this_dist < best_dist:
        best_pos = this_pos
        best_dist = this_dist
        best_well_id = well_id
        if best_dist <= max_dist:
          fastq_well_id_hash[umi_well_seq] = (well_id, best_pos, best_dist)
          if best_dist == 0 and best_pos == FastqReadData.well_id_start:
            break # no need to try other well_ids if we've found a perfect match
  # save FastqReadNgramHash data structure with pickle
  sq.log(f'Saving fastq_well_id_hash to %s...' % fastq_well_id_hash_filename)
  fastq_well_id_hash_file = open(fastq_well_id_hash_filename, 'wb')
  pickle.dump(fastq_well_id_hash, fastq_well_id_hash_file)
  fastq_well_id_hash_file.close()

### Summary 
print(f'fastq_read_ngrams.umi_well_seq_hash: {len(fastq_read_ngrams.umi_well_seq_hash)} items')        
print(f'fastq_well_id_hash: {len(fastq_well_id_hash)} items')

well_id_counts_hash = collections.Counter(well_id for well_id, best_pos, best_dist in fastq_well_id_hash.values())
for well_id in sorted(well_id_counts_hash, key = well_id_counts_hash.get, reverse = True):
  print(f'{well_id}: {well_id_counts_hash[well_id]} items (viral copies: {known_well_id_counts_hash[well_id]})')   

well_id_array = np.array(well_ids).reshape(-1, 1)[1:len(well_ids), ] # we need a 2D array for pdist
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
sorted_umi_well_seqs = sorted(fastq_read_ngrams.umi_well_seq_hash, key = lambda umi_well_seq : fastq_read_ngrams.num_reads(umi_well_seq), reverse = True)#[0:10]
# find matchs for all umi_well_seqs
for umi_well_seq in sorted_umi_well_seqs:
  sq.log(f'Querying with {umi_well_seq}...')
  ngram_matches = fastq_read_ngrams.umi_well_seq_query(umi_well_seq, max_mismatches = FastqReadNgramHash.ngram_length + 2)
  # # sort matches by total number of times each inexact match seen
  ngram_matches = sorted(ngram_matches, key = lambda ngram_match : fastq_read_ngrams.num_reads(ngram_match[0]), reverse = True)
  num_inexact_matches = sum([fastq_read_ngrams.num_reads(match_umi_well_seq) for match_umi_well_seq, num_ngram_matches in ngram_matches if match_umi_well_seq != umi_well_seq])
  if umi_well_seq in fastq_well_id_hash:
    print(f'Summary: {highlight_well_id(umi_well_seq, fastq_well_id_hash[umi_well_seq][1], fastq_well_id_hash[umi_well_seq][2])}  Exact: {fastq_read_ngrams.num_reads(umi_well_seq)} Inexact: {num_inexact_matches}')
  else:
    print(f'Summary: {umi_well_seq}  Exact: {fastq_read_ngrams.num_reads(umi_well_seq)} Inexact: {num_inexact_matches}')
  for match in ngram_matches:
    if match[0] in fastq_well_id_hash:
      print(f'{" "*len("Summary: ")}{highlight_well_id(match[0], fastq_well_id_hash[match[0]][1], fastq_well_id_hash[match[0]][2])}: {fastq_read_ngrams.num_reads(match[0])} {fastq_well_id_hash[match[0]][0]} (hist. int.: {match[1]})')
    else:
      print(f'{" "*len("Summary: ")}{match[0]}: {fastq_read_ngrams.num_reads(match[0])}{" "*(FastqReadData.well_id_length + 1)} (hist. int.: {match[1]})')



