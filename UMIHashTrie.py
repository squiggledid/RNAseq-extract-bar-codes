import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from HashTrie import HashTrie
from UMIData import UMIData

class UMIHashTrie(HashTrie):
  
    def __init__(self):
      super().__init__()
    
    def record_read(self, umi_id, well_id, amplicon_id):
  
      umi_data = self.full_query(umi_id)
      if umi_data:
        umi_data = umi_data[umi_id]['trie_object']
        umi_data.add_read(well_id, amplicon_id)
      else:
        new_umi = UMIData(umi_id, well_id, amplicon_id)
        self.insert(new_umi.umi, new_umi)
