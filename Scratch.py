#! /usr/bin/env python

import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from HashTrie import HashTrie
from UMIData import UMIData 

def record_read(umi_id, well_id, amplicon_id):
  
  umi_data = umi_trie.full_query(umi_id)
  if umi_data:
    umi_data = umi_data[umi_id]['trie_object']
    umi_data.add_read(well_id, amplicon_id)
  else:
    new_umi = UMIData(umi_id, well_id, amplicon_id)
    umi_trie.insert(new_umi.umi, new_umi)
  
  
UMIs = ('ACACACAC', 'ACACAGAC', 'ACTCACAC', 'ACACACGG', 'ABCDEFGH', 'DEFGHIJK', 'BCDEFGHI')
Well_IDs = ('fred', 'bob', 'jane')
Amplicon_IDs = ('cat', 'dog')

umi_trie = HashTrie()

record_read(UMIs[0], Well_IDs[0], Amplicon_IDs[0])
record_read(UMIs[0], Well_IDs[1], Amplicon_IDs[0])
record_read(UMIs[1], Well_IDs[1], Amplicon_IDs[1])
record_read(UMIs[6], Well_IDs[1], Amplicon_IDs[1])
record_read(UMIs[6], Well_IDs[1], Amplicon_IDs[1])

result = umi_trie.full_query_with_errors(UMIs[0], 8)
for umi_id in result.keys():
  print(result[umi_id]['trie_object'])
