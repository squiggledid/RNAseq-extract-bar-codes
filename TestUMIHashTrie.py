#! /usr/bin/env python

import sys

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
umi_trie.record_read(UMIs[2], Well_IDs[1], Amplicon_IDs[1])
umi_trie.record_read(UMIs[3], Well_IDs[1], Amplicon_IDs[1])
umi_trie.record_read(UMIs[4], Well_IDs[1], Amplicon_IDs[1])
umi_trie.record_read(UMIs[5], Well_IDs[1], Amplicon_IDs[1])


result = umi_trie.full_query_with_errors(UMIs[1], 2, destructive = True)
print(f'Search: %s, max_errors: 2' % UMIs[1])
for umi_id in result.keys():
  print(result[umi_id]['trie_object'])
result = umi_trie.full_query_with_errors(UMIs[0], 8)
print(f'Search: %s, max_errors: 8' % UMIs[0])
for umi_id in result.keys():
  print(result[umi_id]['trie_object'])
