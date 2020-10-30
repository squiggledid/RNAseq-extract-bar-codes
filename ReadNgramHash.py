import sys
import Levenshtein

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from ReadData import ReadData
        
class ReadNgramHash:
  '''
  Hash storing ReadData objects, indexed by ngrams, with entries being lists of
  tuples of the corresponding ReadData object and the offset the ngram for that
  object started from
  '''

  # methods
  def __init__(self, seq_length, ngram_length = 6, store_reads = False, build_ngram_histogram_cache = True):
    self.seq_length = seq_length
    self.ngram_length = ngram_length
    self.store_reads = store_reads
    self.build_ngram_histogram_cache = build_ngram_histogram_cache
    self.umi_well_seq_hash = {}
    self.ngram_hash = {}
    self.ngram_histogram_cache = {}

  def __str__(self):
    import collections
    string_rep = ''
    for key in sorted(self.ngram_hash, key = lambda umi_well_seq : self.num_reads(umi_well_seq), reverse = True):
      string_rep += f'{key}: {self.num_reads(key)}\n'
      counts_hash = collections.Counter(offset for umi_well_seq, offset in self.ngram_hash[key])
      string_rep += '\t\n'.join(f'{key}: {counts_hash[key]}' for key in sorted(counts_hash, key = counts_hash.get, reverse = True))
      string_rep += '\n'
    return string_rep
    
  def insert(self, the_read):
    # track the reads where we've seen this UMI/well_id combo
    if the_read.umi_well_seq in self.umi_well_seq_hash:
      if (self.store_reads):
        self.umi_well_seq_hash[the_read.umi_well_seq].append(the_read)
      else:
        self.umi_well_seq_hash[the_read.umi_well_seq] += 1
    else:
      if (self.store_reads):
        self.umi_well_seq_hash[the_read.umi_well_seq] = [the_read]
      else:
        self.umi_well_seq_hash[the_read.umi_well_seq] = 1
      # insert ngrams
      ngram_histogram = {}
      for offset in range(self.seq_length - self.ngram_length):
        ngram = the_read.umi_well_seq[offset:(offset + self.ngram_length)]
        if ngram in self.ngram_hash:
          self.ngram_hash[ngram].append((the_read.umi_well_seq, offset))
        else:
          self.ngram_hash[ngram] = [(the_read.umi_well_seq, offset)]
      if (self.build_ngram_histogram_cache):
        # Insert the umi_well_seq into the ngram_histogram_cache at build time, so it is saved with the hash before any queries
        self.get_ngram_histogram(the_read.umi_well_seq)
      
  def num_reads(self, umi_well_seq):
    if (self.store_reads):
      return(len(self.umi_well_seq_hash[umi_well_seq]))
    else:
      return(self.umi_well_seq_hash[umi_well_seq])

  def get_ngram_histogram(self, umi_well_seq):
    if umi_well_seq not in self.ngram_histogram_cache:
      ngram_histogram = {}
      for offset in range(self.seq_length - self.ngram_length):
        ngram = umi_well_seq[offset:(offset + self.ngram_length)]
        if ngram in ngram_histogram:
          ngram_histogram[ngram] += 1
        else:
          ngram_histogram[ngram] = 1
      self.ngram_histogram_cache[umi_well_seq] = ngram_histogram
    return(self.ngram_histogram_cache[umi_well_seq])
    
  def umi_well_seq_query(self, query_umi_well_seq, min_hist_int = 0.8, sort_result = False):
    min_similarity = 0.8*(self.seq_length - self.ngram_length)
    # build a hash of matches
    match_umi_well_seq_counts = {}
    for offset in range(self.seq_length - self.ngram_length):
      ngram = query_umi_well_seq[offset:(offset + self.ngram_length)]
      for match_umi_well_seq, offset in self.ngram_hash[ngram]:
        if match_umi_well_seq in match_umi_well_seq_counts:
          match_umi_well_seq_counts[match_umi_well_seq] += 1
        else:
          match_umi_well_seq_counts[match_umi_well_seq] = 1
    # compute similarities to query
    match_similarities = [(match_umi_well_seq, \
      histogram_intersection( \
      self.get_ngram_histogram(query_umi_well_seq), \
      self.get_ngram_histogram(match_umi_well_seq)) \
      ) for match_umi_well_seq in match_umi_well_seq_counts if match_umi_well_seq_counts[match_umi_well_seq] > min_similarity]
    final_matches = [(match_umi_well_seq, similarity) for (match_umi_well_seq, similarity) in match_similarities if similarity > min_similarity]
    if sort_result:
      final_matches = sorted(final_matches, key = lambda match : match[1], reverse = True)
    return(final_matches)
          
  def seq_ngrams_query(self, query_seq, max_mismatch = 2):
    '''
    Return a list of tuples of mqtches, each with a reference to a read that contains the
    query sequence query_seq, and the postion at which it was found in that read's umi_well_seq
    '''
    if len(query_seq) < self.ngram_length:
      sys.stderr.write(f'Cannot query a ReadNgramHash with a query string shorter than the ngram length ({self.ngram_length})\n')
    match_umi_well_seqs = {}
    num_ngrams = len(query_seq) - self.ngram_length + 1
    for offset in range(num_ngrams):
      ngram = query_seq[offset:(offset + self.ngram_length)]
      if ngram in self.ngram_hash:
        new_matches = self.ngram_hash[ngram]
      else:
        new_matches = []
      for (match_umi_well_seq, match_offset) in new_matches:
        if match_umi_well_seq in match_the_reads:
          match_the_reads[match_umi_well_seq].update({match_offset: offset})
        else:
          match_the_reads[match_umi_well_seq] = {match_offset: offset}
    final_matches = []
    for match_umi_well_seq, offsets in match_umi_well_seqs.items():
      if (num_ngrams - len(offsets) <= max_mismatch):
        sorted_match_offsets = sorted(offsets.keys())
        # print(query_seq)
        # print(the_read.umi_well_seq)
        # print('\n'.join([f'{match_offset}: {offsets[match_offset]}' for match_offset in sorted_match_offsets]))
        best_run = _longest_consecutive_run(sorted_match_offsets)
        if best_run == None:
          continue
        # print(best_run)
        # print(offsets[sorted_match_offsets[_longest_consecutive_run(sorted_match_offsets)]])
        best_match_start = sorted_match_offsets[best_run] - offsets[sorted_match_offsets[_longest_consecutive_run(sorted_match_offsets)]]
        # print(f'best_match_start: {best_match_start}')
        if best_match_start > self.seq_length - len(query_seq) or best_match_start < 0:
          continue
        best_match_dist = Levenshtein.hamming(query_seq, match_umi_well_seq[best_match_start:best_match_start + len(query_seq)])
        # print(f'the_read: {the_read}\nbest_match_start: {best_match_start}\nbest_match_dist: {best_match_dist}')
        final_matches.append((match_umi_well_seq, best_match_start, best_match_dist))
        # print('#########')
    return(final_matches)
    
def histogram_intersection(hist1, hist2):
  similarity = 0
  for entry in hist1.keys():
    if entry in hist2:
      similarity += min(hist1[entry], hist2[entry])
  return similarity

def _longest_consecutive_run(offsets):
  best_run_start = None
  best_run_length = 0
  current_run_length = 0
  current_run_start = 0
  for i in range(1, len(offsets)):
    if (offsets[i - 1] == offsets[i] - 1):
      current_run_length += 1
      if current_run_length > best_run_length:
        best_run_start = current_run_start
        best_run_length = current_run_length
    else:
      current_run_start = i
      current_run_length = 0
  return(best_run_start) 

def seq_target_query(query_seq, target_seq, expected_pos, max_pos_miss = 6, dist_measure = Levenshtein.hamming):
  '''Exhaustive search for best match for query_seq in in target_seq'''
  best_pos = None
  best_dist = len(query_seq)
  if max_pos_miss > 0:
    positions = range(max(0, expected_pos - max_pos_miss), min(len(target_seq) - len(query_seq), expected_pos + max_pos_miss))
  else:
    positions = [expected_pos]
  for pos in positions:
    dist = dist_measure(query_seq, target_seq[pos:pos + len(query_seq)])
    if dist < best_dist:
      best_pos = pos
      best_dist = dist
      if dist == 0: # no need to continue if we've found a perfect match
        break
  return (best_pos, best_dist)

