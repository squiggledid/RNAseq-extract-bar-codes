import sys
import collections

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from FastqReadData import FastqReadData

class FastqReadNgramHash:
  '''
  Hash storing FastqReadData objects, indexed by ngrams, with entries being lists of
  tuples of the corresponding FastqReadData object and the offset the ngram for that
  object started from
  '''

  ngram_length = 8
  
  def __init__(self):
    self.umi_well_id_hash = {}
    self.ngram_hash = {}

  def __str__(self):
    # return '\n'.join(f'{key}: {len(self.ngram_hash[key])}{"\n\t".join(self.ngram_hash[key])}' for key, val in self.ngram_hash.items())
    string_rep = ''
    for key in sorted(self.ngram_hash, key = lambda x : len(self.ngram_hash[x]), reverse = True):
      string_rep += f'{key}: {len(self.ngram_hash[key])}\n'
      counts_hash = collections.Counter(offset for fastq_read, offset in self.ngram_hash[key])
      string_rep += '\t\n'.join(f'{key}: {counts_hash[key]}' for key in sorted(counts_hash, key = counts_hash.get, reverse = True))
      # string_rep += '\t\n'.join(f'{fastq_read}, {offset}' for fastq_read, offset in self.ngram_hash[key])
      string_rep += '\n'
    return string_rep
    
  def insert(self, fastq_read):
    # track the reads where we've seen this UMI/well_id combo
    if fastq_read.umi_well_seq in self.umi_well_id_hash:
      self.umi_well_id_hash[fastq_read.umi_well_seq].append(fastq_read)
    else:
      self.umi_well_id_hash[fastq_read.umi_well_seq] = [fastq_read]
      # insert ngrams - only do this for unique umi_well_seq
      for offset in range(FastqReadData.umi_length + FastqReadData.well_id_length - FastqReadNgramHash.ngram_length):
        ngram = fastq_read.umi_well_seq[offset:(offset + self.ngram_length)]
        if ngram in self.ngram_hash:
          self.ngram_hash[ngram].append((fastq_read, offset))
        else:
          self.ngram_hash[ngram] = [(fastq_read, offset)]
  
  def query(self, umi_well_seq, max_mismatches = 4, sort_result = False):
    min_matches = FastqReadData.umi_length + FastqReadData.well_id_length - FastqReadNgramHash.ngram_length - max_mismatches
    match_fastq_reads = {}
    for offset in range(FastqReadData.umi_length + FastqReadData.well_id_length - FastqReadNgramHash.ngram_length):
      ngram = umi_well_seq[offset:(offset + FastqReadNgramHash.ngram_length)]
      new_matches = self.ngram_hash[ngram]
      for (fastq_read, offset) in new_matches:
        if fastq_read in match_fastq_reads:
          match_fastq_reads[fastq_read].append(offset)
        else:
          match_fastq_reads[fastq_read] = [offset]
    final_matches = [(key, val) for (key, val) in match_fastq_reads.items() if len(val) > min_matches] # TODO check for off by one
    if sort_result:
      final_matches = sorted(final_matches, key = lambda x : len(x[1]), reverse = True)
    # TODO to we want/need to check order? Only matters in max_mismatches is high, I guess
    return(final_matches)

