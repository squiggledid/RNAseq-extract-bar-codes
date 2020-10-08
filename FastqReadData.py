import re

class FastqReadData:
  '''
  Information associated with a Read from a FASTQ file
  '''
  # constants
  umi_start = 0
  umi_length = 16
  well_id_start = 16
  well_id_length = 8
  umi_well_padding = 4

  def __init__(self, read_id, sequence, quality, store_sequence = False, store_quality = False):
      # Extract the IDs
      self.read_id = read_id
      self.umi_well_seq = sequence[self.umi_start:(self.umi_start + self.umi_length + self.well_id_length + self.umi_well_padding - 1)]
      amplicon_match = re.search(':([^:]+)$', read_id)
      self.amplicon_id: str = amplicon_match.group(1)
      if store_sequence:
        self.sequence = sequence
      if store_quality:
        self.quality = quality

  def __str__(self):
    string_rep = f'FastqReadData: {self.read_id}, umi_well_seq: {self.umi_well_seq}'
    return string_rep

