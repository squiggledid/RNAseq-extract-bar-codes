#! /usr/bin/env python

import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

# from Trie import Trie
from HashTrie import HashTrie

ID1 = "ACACACAC"
ID2 = "ACACAGAC"
ID3 = "ACTCACAC"
ID4 = "ACACACGG"
ID5 = "ABCDEFGH"
ID6 = "DEFGHIJK"
ID7 = "BCDEFGHI"

IDs = HashTrie()

IDs.insert(ID1)
IDs.insert(ID1)
IDs.insert(ID2)
IDs.insert(ID3)
IDs.insert(ID4)
IDs.insert(ID5)
IDs.insert(ID6)
IDs.insert(ID6)
IDs.insert(ID6)
IDs.insert(ID7)

# result = IDs.full_query("ACACAC")
# print(result)
result = IDs.full_query(ID6)
print(result)
result = IDs.full_query_with_errors(ID1, 2)
print(result)

