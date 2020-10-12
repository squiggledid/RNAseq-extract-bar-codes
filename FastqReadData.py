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
  
  # dynamic class variables
  ngram_histogram_cache = {}
  
  def __init__(self, read_id, sequence, quality, store_sequence = False, store_quality = False):
      # Extract the IDs
      self.read_id = read_id
      self.umi_well_seq = sequence[self.umi_start:(self.umi_start + self.umi_length + self.well_id_length + self.umi_well_padding)]
      amplicon_match = re.search(':([^:]+)$', read_id)
      self.amplicon_id: str = amplicon_match.group(1)
      if store_sequence:
        self.sequence = sequence
      if store_quality:
        self.quality = quality
      # initialise ngram cache

  def __str__(self):
    string_rep = f'FastqReadData: {self.read_id}, umi_well_seq: {self.umi_well_seq}'
    return string_rep

  def get_ngram_histogram(self, ngram_length): # takes too much space to save them
    if self.umi_well_seq not in FastqReadData.ngram_histogram_cache:
      seq_length = FastqReadData.umi_length + FastqReadData.well_id_length + FastqReadData.umi_well_padding
      ngram_histogram = {}
      for offset in range(seq_length - ngram_length):
        ngram = self.umi_well_seq[offset:(offset + ngram_length)]
        if ngram in ngram_histogram:
          ngram_histogram[ngram] += 1
        else:
          ngram_histogram[ngram] = 1
      FastqReadData.ngram_histogram_cache[self.umi_well_seq] = ngram_histogram
    return(FastqReadData.ngram_histogram_cache[self.umi_well_seq])
