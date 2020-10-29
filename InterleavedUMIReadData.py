import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from ReadData import ReadData

def _extract_by_pos_and_length(seq, pos_length):
  return(''.join([seq[pos:pos+length] for pos, length in pos_length]))
    
class InterleavedUMIReadData(ReadData):
  '''
  Information associated with a Read from a 'tab' file from with interleaved UMI and Well ID as used by apharseq
  '''
  # constants
  umi_pos_length = [(0, 4), (9, 4), (18, 2)]
  well_id_pos_length = [(4, 5), (13, 5)]
  seq_length = sum([length for pos, length in umi_pos_length] + [length for pos, length in well_id_pos_length])
  
  def __init__(self, R1_seq):
      # Extract the IDs
      self.umi_well_seq = R1_seq[0:self.seq_length] # discard padding 'N's

  def __str__(self):
    string_rep = f'InterleavedUMIReadData: umi_well_seq: {self.umi_well_seq}'
    return string_rep
  
  @classmethod
  def extract_umi_and_well_id(cls, umi_well_seq):
    umi = _extract_by_pos_and_length(umi_well_seq, cls.umi_pos_length)
    well_id = _extract_by_pos_and_length(umi_well_seq, cls.well_id_pos_length)
    return(umi, well_id)
