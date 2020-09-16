import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.insert(0, project_dir)

from Trie import Trie

class HashTrie(Trie):
  '''
  A HashTrie uses a dictionary to store refences to objects stored in the trie.
  This allows constant time exact match retrieval. The trie is only traversed when doing
  inexact matching.
  '''
  
  def __init__(self):
    super().__init__()
    self.entries = {}
    
  def insert(self, word, trie_object = None):
    '''
    Insert an object into the HashTrie
    '''
    if trie_object == None:
      trie_object = word # if nothing to store was passed, store the word itself
    
    if word in self.entries:
      pass # it's already stored, so do nothing
    else:
      self.entries[word] = trie_object
      super().insert(word, trie_object)

  def full_query(self, word):
    if word in self.entries:
      return {word: {
                'trie_object':self.entries[word],
                'num_errors': 0
                }
              }
    else:
      return {} # it's not in the trie, return an empty dictonary
    
  def full_query_with_errors(self, word, max_errors = 0):
    if max_errors == 0:
      return self.full_query(word)
    else:
      return super().full_query_with_errors(word, max_errors)

  def contains(self, word):
    return word in self.entries
  
