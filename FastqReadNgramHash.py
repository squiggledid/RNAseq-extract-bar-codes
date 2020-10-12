import sys
import collections

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from FastqReadData import FastqReadData

def histogram_intersection(hist1, hist2):
  similarity = 0
  for entry in hist1.keys():
    if entry in hist2:
      similarity += min(hist1[entry], hist2[entry])
  return similarity
    
class FastqReadNgramHash:
  '''
  Hash storing FastqReadData objects, indexed by ngrams, with entries being lists of
  tuples of the corresponding FastqReadData object and the offset the ngram for that
  object started from
  '''

  ngram_length = 6
  
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
    seq_length = FastqReadData.umi_length + FastqReadData.well_id_length + FastqReadData.umi_well_padding
    # track the reads where we've seen this UMI/well_id combo
    if fastq_read.umi_well_seq in self.umi_well_id_hash:
      self.umi_well_id_hash[fastq_read.umi_well_seq].append(fastq_read)
    else:
      self.umi_well_id_hash[fastq_read.umi_well_seq] = [fastq_read]
      # insert ngrams
      ngram_histogram = {}
      for offset in range(seq_length - FastqReadNgramHash.ngram_length):
        ngram = fastq_read.umi_well_seq[offset:(offset + self.ngram_length)]
        if ngram in self.ngram_hash:
          self.ngram_hash[ngram].append((fastq_read, offset))
        else:
          self.ngram_hash[ngram] = [(fastq_read, offset)]
    
  def fastq_read_query(self, query_fastq_read, max_mismatches = 4, sort_result = False):
    seq_length = FastqReadData.umi_length + FastqReadData.well_id_length + FastqReadData.umi_well_padding
    min_similarity = seq_length - FastqReadNgramHash.ngram_length - max_mismatches
    # build a hash of matches
    match_fastq_ngram_hash_entries = {}
    for offset in range(seq_length - FastqReadNgramHash.ngram_length):
      ngram = query_fastq_read.umi_well_seq[offset:(offset + FastqReadNgramHash.ngram_length)]
      match_fastq_ngram_hash_entries.update(dict(self.ngram_hash[ngram]))
    # compute similarities to query
    match_similarities = [(match_fastq_read, histogram_intersection(query_fastq_read.get_ngram_histogram(FastqReadNgramHash.ngram_length), match_fastq_read.get_ngram_histogram(FastqReadNgramHash.ngram_length))) for match_fastq_read in match_fastq_ngram_hash_entries.keys()]
    final_matches = [(match_fastq_read, similarity) for (match_fastq_read, similarity) in match_similarities if similarity > min_similarity]
    if sort_result:
      final_matches = sorted(final_matches, key = lambda x : x[1], reverse = True)
    return(final_matches)

  def exact_match_query(self, query_seq):
    '''
    Return a list of tuples of mqtches, each with a reference to a read that contains the
    query sequence query_seq, and the postion at which it was found in that read's umi_well_seq
    '''
    if len(query_seq) < FastqReadNgramHash.ngram_length:
      sys.stderr.write(f'Cannot query a FastqReadNgramHash with a query string shorter than the ngram length ({FastqReadNgramHash.ngram_length})\n')
    match_fastq_reads = {}
    num_ngrams = len(query_seq) - FastqReadNgramHash.ngram_length + 1
    for offset in range(num_ngrams):
      ngram = query_seq[offset:(offset + FastqReadNgramHash.ngram_length)]
      if ngram in self.ngram_hash:
        new_matches = self.ngram_hash[ngram]
      else:
        new_matches = []
      for (fastq_read, offset) in new_matches:
        if fastq_read in match_fastq_reads:
          match_fastq_reads[fastq_read].append(offset)
        else:
          match_fastq_reads[fastq_read] = [offset]
    final_matches = [(key, min(val)) for (key, val) in match_fastq_reads.items() if len(val) == num_ngrams]  
    return(final_matches)
    
  
