#! /usr/bin/env python

import gzip
import re
import sys
import os
import pickle
import colorama
colorama.init()
import pandas as pd

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

import squtils as sq
from FastqReadData import FastqReadData
from FastqReadNgramHash import FastqReadNgramHash

###############################################################################

mismatch_colours = (colorama.Fore.GREEN, colorama.Fore.YELLOW, colorama.Fore.MAGENTA, colorama.Fore.RED)
def highlight_well_id(umi_well_seq, well_id, position, mismatch):
  mismatch = min(mismatch, 3) # >2 should never happen as presently constructed
  return(f'{umi_well_seq[0:position]}{mismatch_colours[mismatch]}{umi_well_seq[position:position + len(well_id)]}{colorama.Style.RESET_ALL}{umi_well_seq[position + len(well_id):]}')

###############################################################################
if (len(sys.argv)) < 2 or (len(sys.argv) > 3):
  sys.exit('A single command line argument specifying the fastq.gz file to process is required, followed by an optional argument specifying experimental metadata. Exiting.')
# first and only command line argument is the fastq.gz file to process
fastq_filename = sys.argv[1]
if len(sys.argv) == 3:
  metadata_filename = sys.argv[2]
else:
  metadata_filename = None

if os.path.exists(metadata_filename):
  well_id_colname = 'Well_Barcode_sequence'
  metadata = pd.read_csv(metadata_filename)
  well_ids = metadata[well_id_colname].to_list()
else:
  well_ids = []
  
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
  max_to_read = 10 # For testing
  if max_to_read and (max_to_read < report_every):
    report_every = max_to_read
  umi_well_id_end = FastqReadData.umi_start + FastqReadData.umi_length + FastqReadData.well_id_length + FastqReadData.umi_well_padding - 1
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
      if ignore_Ns and ('N' in sequence[FastqReadData.umi_start:umi_well_id_end]):
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

### Tests

# Test querying with umi_well_seq data, ordered by well_id

max_well_id_offset = 6
fastq_well_id_hash = {}
sq.log(f'Building fastq_well_id_hash')
for well_id in well_ids:
  well_id_matchs = fastq_read_ngrams.seq_query(well_id)
  for fastq_read, position, mismatch in well_id_matchs:
    if abs(position - FastqReadData.well_id_start) > max_well_id_offset:
      continue # assume this is an erroneous match in the UMI
    if fastq_read in fastq_well_id_hash:
      if fastq_well_id_hash[fastq_read][2] > mismatch \
        or ( \
          fastq_well_id_hash[fastq_read][2] == mismatch \
          and abs(fastq_well_id_hash[fastq_read][1] - FastqReadData.well_id_start) > abs(position - FastqReadData.well_id_start \
        ) \
      ):
        # print('Better match found!')
        fastq_well_id_hash[fastq_read] = (well_id, position, mismatch) # we've found a better well_id for this read
    else:
        fastq_well_id_hash[fastq_read] = (well_id, position, mismatch) # first time we've seen this read

print(f'fastq_well_id_hash: {len(fastq_well_id_hash)} items')    
print(f'fastq_read_ngrams.umi_well_id_hash: {len(fastq_read_ngrams.umi_well_id_hash)} items')        

# print('\n'.join([f'{highlight_well_id(fastq_read.umi_well_seq, well_id, position)}: {position}' for fastq_read, position in fastq_well_id_hash.items()]))     
  # if len(well_id_matchs) > 1:
  #   closest_match_index = min(enumerate(well_id_matchs), key = lambda x : abs(x[1][1] - FastqReadData.well_id_start))[0]
  #   well_id_matchs = well_id_matchs[closest_match_index]

sq.log(f'Sorting umi_well_seqs')
sorted_umi_well_seqs = sorted(fastq_read_ngrams.umi_well_id_hash, key = lambda k : len(fastq_read_ngrams.umi_well_id_hash[k]), reverse = True)#[0:10]
sorted_fastq_reads = [fastq_read_ngrams.umi_well_id_hash[umi_well_seq][0] for umi_well_seq in sorted_umi_well_seqs] # first will do, as all have the same umi_well_seq
# only sort those with exact match well_ids for now
# sorted_fastq_reads = sorted(fastq_well_id_hash.keys(), key = lambda k : len(fastq_read_ngrams.umi_well_id_hash[k.umi_well_seq]), reverse = True)#[0:100]
# find matchs for all umi_well_seqs
for fastq_read in sorted_fastq_reads:
  umi_well_seq = fastq_read.umi_well_seq
  sq.log(f'Querying FastqReadNgramHash with read with {umi_well_seq}')
  ngram_matches = fastq_read_ngrams.fastq_read_query(fastq_read, max_mismatches = FastqReadNgramHash.ngram_length + 2)
  # # sort matches by total number of times each inexact match seen
  ngram_matches = sorted(ngram_matches, key = lambda k : len(fastq_read_ngrams.umi_well_id_hash[k[0].umi_well_seq]), reverse = True)
  num_inexact_matches = sum([len(fastq_read_ngrams.umi_well_id_hash[fastq_read.umi_well_seq]) for fastq_read, num_ngram_matches in ngram_matches if fastq_read.umi_well_seq != umi_well_seq])
  if fastq_read in fastq_well_id_hash:
    print(f'UMI-WellID: {highlight_well_id(umi_well_seq, well_id, fastq_well_id_hash[fastq_read][1], fastq_well_id_hash[fastq_read][2])}\tExact: {len(fastq_read_ngrams.umi_well_id_hash[umi_well_seq])}\tInexact: {num_inexact_matches}')
  else:
    print(f'UMI-WellID: {umi_well_seq}\tExact: {len(fastq_read_ngrams.umi_well_id_hash[umi_well_seq])}\tInexact: {num_inexact_matches}')
  for match in ngram_matches:
    if match[0].umi_well_seq != umi_well_seq:
      # print(match)
      # print(match[0].umi_well_seq)
      # print(fastq_well_id_hash[match[0]][1])
      # print(fastq_read_ngrams.umi_well_id_hash[match[0].umi_well_seq])
      if (fastq_read in fastq_well_id_hash) and (match[0] in fastq_well_id_hash):
        print(f'            {highlight_well_id(match[0].umi_well_seq, fastq_well_id_hash[fastq_read][0], fastq_well_id_hash[match[0]][1], fastq_well_id_hash[match[0]][2])}: {len(fastq_read_ngrams.umi_well_id_hash[match[0].umi_well_seq])} (histogram intersection: {match[1]})')
      else:
        print(f'            {match[0].umi_well_seq}: {len(fastq_read_ngrams.umi_well_id_hash[match[0].umi_well_seq])} (histogram intersection: {match[1]})')

# # Test querying with well_id data
# well_id = 'TCCATAGG'
# well_id_matchs = fastq_read_ngrams.exact_match_query(well_id)
# print(f'{len(well_id_matchs)} exact matches for well_id {well_id}:')
# for fastq_read, position in well_id_matchs:
#   umi_well_seq = fastq_read.umi_well_seq
#   print(f'{highlight_well_id(umi_well_seq, well_id, position)}: position')


