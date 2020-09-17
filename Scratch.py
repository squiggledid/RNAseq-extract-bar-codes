#! /usr/bin/env python

import sys
import pickle

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from UMIHashTrie import UMIHashTrie
from UMIData import UMIData

umi_trie_filename = '/Users/davids/OneDrive - The University of Melbourne/WEHI Rory CoVid19 RNA/DB CloudStor 20200812/Standard_R1.fastq.gz_UMIHashTrie.pkl'

# load UMIHashTrie from file
umi_trie_file = open(umi_trie_filename, 'rb')
umi_trie = pickle.load(umi_trie_file)
umi_trie_file.close()

# get UMIs in descending order of count
sorted_umis = sorted(umi_trie.entries.items(), key = lambda x : x[1].count, reverse = True)

# get UMI count distribution
counts = [(x[0], x[1].count) for x in sorted_umis]
# write this as a CSV
csv_filename = umi_trie_filename + '_UMI-counts.csv'
csv_file = open(csv_filename, 'w')
for key, count in counts:
  csv_file.write(f'%s,%d\n' % (key, count))
csv_file.close()
  
# write a CSV for a number of neighbours distribution
max_neighbours = 4
min_count = 100 # from eyeballing a plot, at present
if counts[-1][1] < min_count:
  max_row = [n for n, x in enumerate( [x[1] for x in counts] ) if x < min_count][0]
else:
  max_row = len(sorted_umis)
# max_row = 2 # test
dist_headers = ['dist' + str(x) for x in range(1, max_neighbours + 1)]
num_neighbours_filename = umi_trie_filename + '_UMI-num-neighbours.csv'
num_neighbours_file = open(num_neighbours_filename, 'w')
# print header
num_neighbours_file.write('UMI,count,' + ','.join(dist_headers) + '\n')

# print data rows
for row in range(max_row):
  num_neighbours = [len(umi_trie.full_query_with_errors(sorted_umis[row][0], col)) - 1 for col in range(1, 5)]
  row_string = ','.join([sorted_umis[row][0]] + \
    [str(counts[row][1])] + \
    [str(len(umi_trie.full_query_with_errors(sorted_umis[row][0], col)) - 1) for col in range(1, 5)]
  )
  num_neighbours_file.write(row_string + '\n')
  if row % 100 == 0:
    print(f'Row %d of %d: %s' % (row, max_row, row_string))
    num_neighbours_file.flush() # for progress monitoring
num_neighbours_file.close()
quit()

# print the largest one
print(sorted_umis[0][0])
# get the UMIs within a distance of one of this
print(umi_trie.full_query_with_errors(sorted_umis[0][0], 1))


print('Test with umi_trie2')
result = umi_trie2.full_query_with_errors(UMIs[0], 8)
for umi_id in result.keys():
  print(result[umi_id]['trie_object'])
