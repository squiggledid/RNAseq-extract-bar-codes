#! /usr/bin/env python

import sys
import pickle

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from UMIHashTrie import UMIHashTrie
from UMIData import UMIData
  
UMIs = ('ACACACAC', 'ACACAGAC', 'ACTCACAC', 'ACACACGG', 'ABCDEFGH', 'DEFGHIJK', 'BCDEFGHI')
Well_IDs = ('fred', 'bob', 'jane')
Amplicon_IDs = ('cat', 'dog')

umi_trie = UMIHashTrie()

umi_trie.record_read(UMIs[0], Well_IDs[0], Amplicon_IDs[0])
umi_trie.record_read(UMIs[0], Well_IDs[1], Amplicon_IDs[0])
umi_trie.record_read(UMIs[1], Well_IDs[1], Amplicon_IDs[1])
umi_trie.record_read(UMIs[6], Well_IDs[1], Amplicon_IDs[1])
umi_trie.record_read(UMIs[6], Well_IDs[1], Amplicon_IDs[1])

result = umi_trie.full_query_with_errors(UMIs[0], 8)
for umi_id in result.keys():
  print(result[umi_id]['trie_object'])
  
# save UMITrie data structure with pickle
umi_trie_filename = 'test_umi_trie_save.pkl'
umi_trie_file = open(umi_trie_filename, 'wb')
pickle.dump(umi_trie, umi_trie_file)
umi_trie_file.close()
# reload UMITrie from file
umi_trie_file = open(umi_trie_filename, 'rb')
umi_trie2 = pickle.load(umi_trie_file)
umi_trie_file.close()
# test
print('Test with umi_trie2')
result = umi_trie2.full_query_with_errors(UMIs[0], 8)
for umi_id in result.keys():
  print(result[umi_id]['trie_object'])
