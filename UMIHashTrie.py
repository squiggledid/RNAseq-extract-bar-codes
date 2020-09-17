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

    def _destructive_full_query_dfs_with_errors(self, node, query_word, partial_result, num_errors, max_errors, level):
        """Do depth first searches until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        When a word is matched according to the criterion, it is removed from the trie,
        and this method returns True, signalling to the caller to remove the token that
        led to the match from the children array
        """
        #print('\t'*level, 'query_word: ', query_word)
        #print('\t'*level, 'node.token: ', node.token)
        #print('\t'*level, 'partial_result: ', partial_result)
        #print('\t'*level, 'num_errors: ', num_errors)
        #print('\t'*level, 'max_errors: ', max_errors, '\n')
        
        if num_errors > max_errors:
            #print('\t'*level, 'Too many errors\n')
            return False
        
        if node.is_end:
            #print('\t'*level, 'Found end.\n')
            self.result[partial_result] = {
                'trie_object': node.trie_object,
                'num_errors': num_errors
            }
            return True
            
        if len(query_word) < 1:
            #print('\t'*level, 'Empty word.\n')
            return False

        query_token = query_word[0]
        tokens_to_destroy = [] # note, with inexact matching, there can be more than one path that leads to a match
        for token in node.children:
            #print('\t'*level, 'Level: ', level, 'token: ', token, '\n')
            if token != query_token:
                error_increment = 1
            else:
                error_increment = 0
            destroy_token = self._destructive_full_query_dfs_with_errors(
                node = node.children[token],
                query_word = query_word[1:],
                partial_result = partial_result + token,
                num_errors = num_errors + error_increment,
                max_errors = max_errors,
                level = level + 1)
            if destroy_token:
                tokens_to_destroy.append(token)
        #print('\t'*level, 'Tokens to destroy:', tokens_to_destroy)
        for token in tokens_to_destroy:
          del node.children[token]

    def full_query_with_errors(self, word, max_errors, destructive = False):
        """Given an input (a word), retrieve all words that match that word from the trie,
        with max_errors allowed.
        When a word is matched according to the criterion, it is removed from the trie
        """
        
        #print('word: ', word, '\n')
        
        query_token = word[0]
        num_errors = 0
        self.result = {}
        level = 0
        full_query_dfs = self._destructive_full_query_dfs_with_errors if destructive else self._full_query_dfs_with_errors
        for token in self.root.children:
            #print('Level: ', level, 'token: ', token, '\n')
            if token != query_token:
                error_increment = 1
            else:
                error_increment = 0
            full_query_dfs(
                node = self.root.children[token],
                query_word = word[1:],
                partial_result = token,
                num_errors = num_errors + error_increment,
                max_errors = max_errors,
                level = level + 1)
        return self.result
