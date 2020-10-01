#! /usr/bin/env python

import sys
import pickle
import collections

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from UMIHashTrie import UMIHashTrie
from UMIData import UMIData

if len(sys.argv) != 2:
    sys.exit('A single command line argument specifying the umi-trie .pkl file to process is required. Exiting.')
# first and only command line argument is the fastq.gz file to process
umi_trie_filename = sys.argv[1]

# load UMIHashTrie from file
umi_trie_file = open(umi_trie_filename, 'rb')
umi_trie = pickle.load(umi_trie_file)
umi_trie_file.close()

# get UMIs in descending order of count
sorted_umis = sorted(umi_trie.entries.items(), key = lambda x : x[1].count, reverse = True)

# Example query
matches = umi_trie.full_query_with_errors(sorted_umis[0][0], 1)
sorted_matches = sorted(matches.values(), key = lambda x : x['trie_object'].count, reverse = True)
for match in sorted_matches:
  print(f'Num match errors: ', match['num_errors'])
  print(match['trie_object'])
  
# get UMI count distribution
counts = [(x[0], x[1].count) for x in sorted_umis]
# write this as a CSV
csv_filename = umi_trie_filename + '_UMI-counts.csv'
csv_file = open(csv_filename, 'w')
for key, count in counts:
  csv_file.write(f'%s,%d\n' % (key, count))
csv_file.close()
  
# # write a CSV for a number of neighbours distribution
# max_neighbours = 4
# min_count = 100 # from eyeballing a plot, at present
# if counts[-1][1] < min_count:
#   max_row = [n for n, x in enumerate( [x[1] for x in counts] ) if x < min_count][0]
# else:
#   max_row = len(sorted_umis)
# # max_row = 2 # test
# dist_headers = ['dist' + str(x) for x in range(1, max_neighbours + 1)]
# num_neighbours_filename = umi_trie_filename + '_UMI-num-neighbours.csv'
# num_neighbours_file = open(num_neighbours_filename, 'w')
# # print header
# num_neighbours_file.write('UMI,count,' + ','.join(dist_headers) + '\n')
# 
# # print data rows
# for row in range(max_row):
#   matches = umi_trie.full_query_with_errors(sorted_umis[row][0], max_neighbours)
#   distances = [x['num_errors'] for x in matches.values()]
#   dist_counts_hash = collections.Counter(distances)
#   dist_counts = [str(dist_counts_hash[dist]) for dist in range(1, max_neighbours + 1)] # ignore first element, which is distance 0 from self
#   row_string = ','.join([sorted_umis[row][0]] + [str(counts[row][1])] + dist_counts)
#   num_neighbours_file.write(row_string + '\n')
#   if row % 100 == 0:
#     print(f'Row %d of %d: %s' % (row, max_row, row_string))
#     num_neighbours_file.flush() # for progress monitoring with tail -f
# num_neighbours_file.close()

# write a CSV for a destructive number of neighbours distribution
max_neighbours = 8
min_count = 100 # from eyeballing a plot, at present
if counts[-1][1] < min_count:
  max_row = [n for n, x in enumerate( [x[1] for x in counts] ) if x < min_count][0]
else:
  max_row = len(sorted_umis)
# max_row = 2 # test
dist_headers = ['dist' + str(x) for x in range(1, max_neighbours + 1)]
num_neighbours_filename = umi_trie_filename + '_UMI-num-neighbours-' + str(max_neighbours) + '_destructive.csv'
num_neighbours_file = open(num_neighbours_filename, 'w')
# print header
num_neighbours_file.write('UMI,count,' + ','.join(dist_headers) + '\n')

# compute distances one distance slice at a time, since is destructive (i.e. an n+1_neighbour of a centre could be the n_neighbour of another)
distance_cols = []
for dist in range(1, max_neighbours + 1):
  print(f'Finding %d-neighbours, and destroying...' % dist)
  matches_col = [umi_trie.full_query_with_errors(sorted_umis[row][0], dist, destructive = True) for row in range(max_row)]
  distances_col = [[x['num_errors'] for x in matches.values() if x['num_errors'] != 0] for matches in matches_col]
  dist_counts_col = [str(len(matches)) for matches in distances_col] # no need to count separately as destruction means there can't be any closer than dist
  distance_cols.append(dist_counts_col)
print('Transposing distance data...')
distance_rows = [list(distance_col) for distance_col in zip(*distance_cols)] # transpose
# print data rows
for row in range(len(distance_rows)):
  row_string = ','.join([sorted_umis[row][0]] + \
    [str(counts[row][1])] + \
    distance_rows[row]
  )
  num_neighbours_file.write(row_string + '\n')
  if row % 100 == 0:
    print(f'Row %d of %d: %s' % (row, max_row, row_string))
    num_neighbours_file.flush() # for progress monitoring with tail -f
num_neighbours_file.close()
