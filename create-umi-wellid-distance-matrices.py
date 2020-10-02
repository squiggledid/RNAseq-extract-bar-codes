#! /usr/bin/env python

import sys
import os # for os.path.join(...), for examnple
import time
import pickle
import collections
import numpy as np
from scipy.spatial.distance import pdist, squareform
import Levenshtein # By far the fastest library tried for this - including for Hamming distance. Trued textdistance, scipy.spatial.distance.hamming

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from UMIHashTrie import UMIHashTrie
from UMIData import UMIData

if len(sys.argv) != 2:
    sys.exit('A single command line argument specifying the umi-trie .pkl file to process is required. Exiting.')
# first and only command line argument is the fastq.gz file to process
umi_trie_filename = sys.argv[1]

def log(s):
  print(time.asctime() + ":", s)
  
# load UMIHashTrie from file
log('Loading UMIHashTrie from ' + umi_trie_filename)
umi_trie_file = open(umi_trie_filename, 'rb')
umi_trie = pickle.load(umi_trie_file)
umi_trie_file.close()

# get UMIs in descending order of count
log('Sorting UMIs')
sorted_umis = sorted(umi_trie.entries.items(), key = lambda x : x[1].count, reverse = True)

# get UMI count distribution
umi_counts = [(x[0], x[1].count) for x in sorted_umis]

# find the last row with at least min_count UMI
min_count = 100 # from eyeballing a plot, at present
if umi_counts[-1][1] < min_count:
  max_row = [n for n, x in enumerate( [x[1] for x in umi_counts] ) if x < min_count][0]
else:
  max_row = len(sorted_umis)
# max_row = 4 # test
umis = [umi for (umi, count) in umi_counts if not 'N' in umi]

umi_trie = None # make memory available to GC
# compute distance matrices for UMIs, and save them to .npy files
umis_array = np.array(umis).reshape(-1, 1)[1:max_row, ] # we need a 2D array for pdist
log(f'Computing UMI Hamming distance matrix for %d UMIs' % max_row)
umi_hamming_distance_matrix = pdist(umis_array, lambda umi1, umi2 : Levenshtein.hamming(umi1[0], umi2[0])) # need to wrap the distance function
np.save(umi_trie_filename + '_hamming.npy', squareform(umi_hamming_distance_matrix))
umi_hamming_distance_matrix = None # make memory available to GC
log(f'Computing UMI Levenshtein distance matrix for %d UMIs' % max_row)
umi_levenshtein_distance_matrix = pdist(umis_array, lambda umi1, umi2 : Levenshtein.distance(umi1[0], umi2[0])) # need to wrap the distance function
np.save(umi_trie_filename + '_levenshtein.npy', squareform(umi_levenshtein_distance_matrix))
umi_levenshtein_distance_matrix = None # make memory available to GC
umis_array = None

