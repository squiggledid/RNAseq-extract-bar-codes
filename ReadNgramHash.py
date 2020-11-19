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
    # testing
    self.hist_match_diff = 0
    self.num_hist_ints = 0

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
      for offset in range(self.seq_length - self.ngram_length):
        ngram = the_read.umi_well_seq[offset:(offset + self.ngram_length)]
        if not (ngram in self.ngram_hash):
          self.ngram_hash[ngram] = {the_read.umi_well_seq : [offset]}
        else:
          if the_read.umi_well_seq in self.ngram_hash[ngram]:
            self.ngram_hash[ngram][the_read.umi_well_seq].append(offset)
          else:
            self.ngram_hash[ngram].update({the_read.umi_well_seq : [offset]})

      if (self.build_ngram_histogram_cache):
        # Insert the umi_well_seq into the ngram_histogram_cache at build time, so it is saved with the hash before any queries
        self.get_ngram_histogram(the_read.umi_well_seq)
  
  def delete(self, umi_well_seq):
    del self.umi_well_seq_hash[umi_well_seq] # NB the key umi_well_seq *must* have been inserted before this is called
    for offset in range(self.seq_length - self.ngram_length):
      ngram = umi_well_seq[offset:(offset + self.ngram_length)]
      if umi_well_seq in self.ngram_hash[ngram]: # if there were multiple entries for the ngram for the same umi_well_seq, it could already have been deleted
        del self.ngram_hash[ngram][umi_well_seq]
        if len(self.ngram_hash[ngram]) == 0:
          self.ngram_hash[ngram]
    
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
    
  def umi_well_seq_query(self, query_umi_well_seq, min_similarity_fraction = 0.8, sort_result = False):
    min_similarity = min_similarity_fraction*(self.seq_length - self.ngram_length)
    # build a hash of matches
    match_umi_well_seq_counts = {}
    for offset in range(self.seq_length - self.ngram_length):
      ngram = query_umi_well_seq[offset:(offset + self.ngram_length)]
      for match_umi_well_seq in self.ngram_hash[ngram]:
        if match_umi_well_seq in match_umi_well_seq_counts:
          match_umi_well_seq_counts[match_umi_well_seq] += len(self.ngram_hash[ngram][match_umi_well_seq]) #TODO if we're never going to use the offsets, just make the hash-hash entry a count?
        else:
          match_umi_well_seq_counts[match_umi_well_seq] = len(self.ngram_hash[ngram][match_umi_well_seq])
    # compute similarities to query
    matches = [(match_umi_well_seq, match_umi_well_seq_counts[match_umi_well_seq]) \
      for match_umi_well_seq in match_umi_well_seq_counts if match_umi_well_seq_counts[match_umi_well_seq] > min_similarity]
    if sort_result:
      matches = sorted(matches, key = lambda match : match[1], reverse = True)
    return(matches)

  def umi_well_seq_query_with_histogram_intersection(self, query_umi_well_seq, min_similarity_fraction = 0.8, sort_result = False):
    min_similarity = min_similarity_fraction*(self.seq_length - self.ngram_length)
    preliminary_matches = self.umi_well_seq_query(query_umi_well_seq, min_similarity_fraction = min_similarity_fraction, sort_result = False) # Don't want to sort preliminary matches, no matter what.
    match_similarities = [(match_umi_well_seq, \
      histogram_intersection( \
      self.get_ngram_histogram(query_umi_well_seq), \
      self.get_ngram_histogram(match_umi_well_seq)) \
      ) for match_umi_well_seq, match_umi_well_seq_count in preliminary_matches]
    final_matches = [(match_umi_well_seq, similarity) for (match_umi_well_seq, similarity) in match_similarities if similarity > min_similarity]
    self.num_hist_ints += len(preliminary_matches)
    self.hist_match_diff += len(preliminary_matches) - len(final_matches)
    if sort_result:
      final_matches = sorted(final_matches, key = lambda match : match[1], reverse = True)
    return(final_matches)
    
  def seq_ngram_query(self, query_seq, min_similarity_fraction = 0.8, sort_result = False):
    # TODO replace umi_well_seq_query with this - at almost zero cost
    min_similarity = min_similarity_fraction*(len(query_seq) - self.ngram_length)
    # build a hash of matches
    match_umi_well_seq_counts = {}
    for offset in range(len(query_seq) - self.ngram_length):
      ngram = query_seq[offset:(offset + self.ngram_length)]
      for match_umi_well_seq in self.ngram_hash[ngram]:
        if match_umi_well_seq in match_umi_well_seq_counts:
          match_umi_well_seq_counts[match_umi_well_seq] += len(self.ngram_hash[ngram][match_umi_well_seq]) #TODO if we're never going to use the offsets, just make the hash-hash entry a count?
        else:
          match_umi_well_seq_counts[match_umi_well_seq] = len(self.ngram_hash[ngram][match_umi_well_seq])
    # compute similarities to query
    matches = [(match_umi_well_seq, match_umi_well_seq_counts[match_umi_well_seq]) \
      for match_umi_well_seq in match_umi_well_seq_counts if match_umi_well_seq_counts[match_umi_well_seq] > min_similarity]
    if sort_result:
      matches = sorted(matches, key = lambda match : match[1], reverse = True)
    return(matches)
    
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
      if dist == 0: # no need to continue if we've found a perfect match TODO consider weighting dist_measure vs miss in pos?
        break
  return (best_pos, best_dist)

