#! /usr/bin/env python

import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.append(project_dir)

from ID_Trie import ID_Trie # TODO learn the proper Python way to do this

ID1 = "ACACACAC"
ID2 = "ACACAGAC"
ID3 = "ACTCACAC"
ID4 = "ACACACGG"
ID5 = "ABCDEFGH"

IDs = ID_Trie()

IDs.insert(ID1)
IDs.insert(ID1)
IDs.insert(ID2)
IDs.insert(ID3)
IDs.insert(ID4)
IDs.insert(ID5)

result = IDs.full_query("ACACAC")
print(result)
result = IDs.full_query(ID5)
print(result)
