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
from ReadNgramHash import ReadNgramHash, seq_target_query

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

ngram_length = ReadNgramHash(0).ngram_length # create a dummy object so we can access default ngram_length
# ngram_length = 4
### Read FASTQ reads
fastq_read_ngrams = ReadNgramHash(FastqReadData.seq_length) # create an object so we can access the ngram_length
fastq_read_ngram_hash_filename = f'{fastq_filename}_ReadNgramHash_{ngram_length}_{FastqReadData.umi_well_padding}.pkl'
if os.path.exists(fastq_read_ngram_hash_filename):
  sq.log(f'Reading data from {fastq_read_ngram_hash_filename}...')
  fastq_read_ngram_hash_file = open(fastq_read_ngram_hash_filename, 'rb')
  fastq_read_ngrams = pickle.load(fastq_read_ngram_hash_file)
  fastq_read_ngram_hash_file.close()
else:
  sq.log(f'Building ReadNgramHash from {fastq_filename}...')
  # Go through FASTQ file, one four-line block at a time
  ignore_Ns = True
  n_read = 0
  n_skipped = 0
  report_every = 100000
  max_to_read = None # None for no limit :)
  # max_to_read = 100000 # For testing
  if max_to_read and (max_to_read < report_every):
    report_every = max_to_read
  umi_well_seq_end = FastqReadData.umi_start + FastqReadData.umi_length + FastqReadData.well_id_length + FastqReadData.umi_well_padding - 1
  fastq_read_ngrams = ReadNgramHash(seq_length = FastqReadData.seq_length, ngram_length = ngram_length);
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
  
  # save ReadNgramHash data structure with pickle
  sq.log(f'Saving ReadNgramHash to %s...' % fastq_read_ngram_hash_filename)
  fastq_read_ngram_hash_file = open(fastq_read_ngram_hash_filename, 'wb')
  pickle.dump(fastq_read_ngrams, fastq_read_ngram_hash_file)
  fastq_read_ngram_hash_file.close()

### Load or build well_id hash
max_well_id_offset = 6
max_dist = 2
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
      (this_pos, this_dist) = seq_target_query(well_id, umi_well_seq, FastqReadData.well_id_start, max_well_id_offset) # , dist_measure = Levenshtein.distance)
      if this_dist < best_dist:
        best_pos = this_pos
        best_dist = this_dist
        best_well_id = well_id
        if best_dist <= max_dist:
          fastq_well_id_hash[umi_well_seq] = (well_id, best_pos, best_dist)
          if best_dist == 0 and best_pos == FastqReadData.well_id_start:
            break # no need to try other well_ids if we've found a perfect match
  # save well_id hash data structure
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

### See if there is any point to doing histogram intersection with this data
# list_of_bin_counts_lists = [list(ngram_hash.values()) for ngram_hash in fastq_read_ngrams.ngram_histogram_cache.values()]
# bin_count_counts = collections.Counter(sq.flatten(list_of_bin_counts_lists))
# print(bin_count_counts)
# for 6-grams, and well_padding of 22, we get:
# Counter({1: 26296303, 2: 130726, 3: 5821, 4: 1916, 5: 1853, 6: 747, 8: 399, 7: 368, 9: 169, 35: 93, 10: 84, 29: 73, 30: 64, 31: 48, 36: 45, 33: 44, 34: 43, 11: 31, 32: 29, 12: 23, 25: 21, 27: 19, 28: 19, 13: 18, 23: 17, 24: 16, 26: 14, 37: 13, 20: 12, 19: 9, 15: 9, 18: 8, 21: 8, 17: 7, 16: 7, 14: 5, 22: 4, 38: 4})
# for 4-grams, and well_padding of 22, we get:
# Counter({1: 23341910, 2: 1935512, 3: 143798, 4: 35539, 5: 10486, 6: 4115, 7: 2297, 8: 1948, 9: 1210, 10: 896, 11: 620, 12: 395, 13: 249, 14: 131, 33: 112, 15: 103, 37: 93, 34: 76, 35: 48, 38: 45, 16: 45, 36: 44, 29: 34, 31: 25, 32: 25, 30: 25, 17: 23, 18: 18, 19: 15, 39: 13, 25: 13, 24: 12, 22: 11, 23: 9, 27: 9, 20: 8, 28: 8, 26: 7, 21: 5, 40: 4})
# ... so perhaps it is worth it. Will have to check to see how often, if even, it causes matches to be excluded
# quit()

### Tests
use_histogram_intersection = True
require_good_well_ids = True
min_exact_match_well_id_frac = 0.5
exclude_seq = 'GACCATTTCACAGATC' # characteristic of the primer dimer problem?
min_exclude_dist = 2
sq.log(f'Sorting umi_well_seqs')
sorted_umi_well_seqs = sorted(fastq_read_ngrams.umi_well_seq_hash, key = lambda umi_well_seq : fastq_read_ngrams.num_reads(umi_well_seq), reverse = True)#[0:10]
# find matchs for all umi_well_seqs
query_num = 1
for umi_well_seq in sorted_umi_well_seqs:
  sq.log(f'Checking {umi_well_seq} against {exclude_seq}...')
  best_exclude_seq_match_pos, best_exclude_seq_match_dist = seq_target_query(exclude_seq, umi_well_seq, \
    expected_pos = 0, max_pos_miss = len(umi_well_seq), dist_measure = Levenshtein.distance)
  print(f'best_exclude_seq_match_pos, best_exclude_seq_match_dist: {best_exclude_seq_match_pos, best_exclude_seq_match_dist}')
  if best_exclude_seq_match_dist < min_exclude_dist:
    print('Skipping due to match with exclude_seq')
    continue # skip umi_well_seqs that contain a sequence sufficiently close to the exclude_seq
  sq.log(f'Querying with {umi_well_seq}...')
  if use_histogram_intersection:
    ngram_matches = fastq_read_ngrams.umi_well_seq_query_with_histogram_intersection(umi_well_seq)
  else:
    ngram_matches = fastq_read_ngrams.umi_well_seq_query(umi_well_seq)
  num_total_matches = sum([fastq_read_ngrams.num_reads(match_umi_well_seq) for match_umi_well_seq, num_ngram_matches in ngram_matches])
  print(f'Num. times histogram intersection made a difference: {fastq_read_ngrams.hist_match_diff}/{fastq_read_ngrams.num_hist_ints}')
  ### check if most of the matches have good well_ids
  if require_good_well_ids:
    match_well_ids_and_counts = [(fastq_well_id_hash[match_umi_well_seq], fastq_read_ngrams.num_reads(match_umi_well_seq)) for match_umi_well_seq, num_ngram_matches in ngram_matches]
    match_well_ids = list(set(match_well_id_and_count[0][0] for match_well_id_and_count in match_well_ids_and_counts))
    well_id_dist_counts = {match_well_id:{} for match_well_id in match_well_ids}
    for well_id_dist_count in well_id_dist_counts:
      well_id_dist_counts[well_id_dist_count] = {dist:0 for dist in range(FastqReadData.well_id_length)}
    for match_well_id_and_count in match_well_ids_and_counts:
      well_id_dist_counts[match_well_id_and_count[0][0]][match_well_id_and_count[0][2]] += match_well_id_and_count[1]
    # find the distance with the maximum count for each matched well_id
    max_dist_counts = [(well_id, dist, well_id_dist_counts[well_id][dist]) for well_id, dist in \
      [(well_id, max(well_id_dist_counts[well_id], key = lambda dist: well_id_dist_counts[well_id][dist])) for well_id in well_id_dist_counts] ]
    mode_well_id, mode_well_id_dist, mode_well_id_count = max(max_dist_counts, key = lambda max_dist_count: max_dist_count[2])
    if mode_well_id_dist == 0:
      exact_match_well_id_frac = mode_well_id_count/num_total_matches
    else:
      exact_match_well_id_frac = 0
    print(f'Fraction of exact well_id matches: {exact_match_well_id_frac}')
    if exact_match_well_id_frac < min_exact_match_well_id_frac:
      print(f'Skipping: exact_match_well_id_frac = {exact_match_well_id_frac}')
      continue
  # sort matches by total number of times each inexact match seen
  ngram_matches = sorted(ngram_matches, key = lambda ngram_match : fastq_read_ngrams.num_reads(ngram_match[0]), reverse = True)
  num_exact_matches = fastq_read_ngrams.num_reads(umi_well_seq)
  num_inexact_matches = num_total_matches - num_exact_matches
  prefix = f'{query_num}. Query: '
  suffix = f'  Exact: {num_exact_matches} Inexact: {num_inexact_matches}'
  if umi_well_seq in fastq_well_id_hash:
    umi_well_seq_string = highlight_well_id(umi_well_seq, fastq_well_id_hash[umi_well_seq][1], fastq_well_id_hash[umi_well_seq][2])
  else:
    umi_well_seq_string = umi_well_seq
  print(f'{prefix}{umi_well_seq_string}{suffix}')
  padding_len = len(prefix)
  for match in ngram_matches:
    if match[0] in fastq_well_id_hash:
      match_umi_well_seq_string = highlight_well_id(match[0], fastq_well_id_hash[match[0]][1], fastq_well_id_hash[match[0]][2])
      well_id_string = fastq_well_id_hash[match[0]][0]
    else:
      match_umi_well_seq_string = match[0]
      well_id_string = f'{" "*(FastqReadData.well_id_length + 1)}'
    print(f'{" "*padding_len}{match_umi_well_seq_string}: {fastq_read_ngrams.num_reads(match[0]):5} {well_id_string} (sim.: {match[1]})')
  query_num += 1


