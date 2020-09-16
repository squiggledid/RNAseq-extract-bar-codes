import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.append(project_dir)

from TrieNode import TrieNode # TODO learn the proper Python way to do this

class Trie(object):
    """A trie class"""

    def __init__(self):
        """
        The trie needs to have a root node.
        The root node does not store any token
        """
        self.root = TrieNode("")
    
    def insert(self, word, trie_object = None):
        """Insert a word into the trie"""
        node = self.root
        
        if trie_object == None:
            trie_object = word # if nothing to store was passed, store the word itself
            
        # Loop through each token in the word (e.g. characters, bytes)
        # Check if there is no child containing the token, create a new child for the current node
        for token in word:
            if token in node.children:
                node = node.children[token]
            else:
                # If a token is not found,
                # create a new node in the trie
                new_node = TrieNode(token)
                node.children[token] = new_node
                node = new_node
        
        node.is_end = True # Mark the end of a word
        node.trie_object = trie_object # store a reference to the object represented by this path through the trie
        
    def _depth_first_search(self, node, prefix):
        """Depth-first traversal of the trie
        
        Args:
            - node: the node to start with
            - prefix: the current prefix, for tracing a
                word while traversing the trie
        """
        if node.is_end:
            self.output.append((prefix + node.token, node.trie_object))
        
        for child in node.children.values():
            self._depth_first_search(child, prefix + node.token)
        
    def prefix_query(self, prefix):
        """Given an input (a prefix), retrieve all words stored in
        the trie with that prefix, sort the words by the number of 
        times they have been inserted
        """
        # Use a variable within the class to keep all possible outputs
        # As there can be more than one word with such prefix
        self.output = []
        node = self.root
        
        # Check if the prefix is in the trie
        for token in prefix:
            if token in node.children:
                node = node.children[token]
            else:
                # cannot find the prefix, return empty list
                return []
        
        # Traverse the trie to get all candidates
        self._depth_first_search(node, prefix[:-1])

        # Sort the results in reverse order and return
        return sorted(self.output, key=lambda prefix: prefix[1], reverse=True)

    def _full_query_dfs(self, node, word, partial_result):
        """Do a depth first search until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        """
            
        print('word: ', word)
        print('node.token: ', node.token)
        print('partial_result: ', partial_result, '\n')

        if node.is_end:
            print('Found end.\n')
            return (partial_result + node.token, node.trie_object) # Success!
        if len(word) < 1:
            print('Empty word.\n')
            return[] # ran out of characters before finding a match

        token = word[0]
        if token in node.children:
            return self._full_query_dfs(node.children[token], word[1:], partial_result + node.token)
        else:
            # print('Not in trie.\n')
            return [] # if is not in trie

    def _full_query_dfs(self, node, word, partial_result):
        """Do depth first search until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        """
            
        # print('word: ', word)
        # print('node.token: ', node.token)
        # print('partial_result: ', partial_result)

        if node.is_end:
            # print('Found end.\n')
            return (partial_result + node.token, node.trie_object) # Success!
        if len(word) < 1:
            # print('Empty word.\n')
            return[] # ran out of characters before finding a match

        token = word[0]
        if token in node.children:
            return self._full_query_dfs(node.children[token], word[1:], partial_result + node.token)
        else:
            # print('Not in trie.\n')
            return [] # if is not in trie

    def full_query(self, word):
        """Given an input (a word), retrieve that word from the trie 
        """
        
        # print('word: ', word)
        
        node = self.root
        token = word[0]
        if token in node.children:
            node = node.children[token]
        else:
            # cannot find the first character of id, return empty list
            return []
        
        # Traverse the trie to find the id
        return self._full_query_dfs(node, word[1:], '')
        
    def _full_query_dfs_with_errors(self, node, query_word, partial_result, num_errors, max_errors, level):
        """Do depth first searches until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        """
        # print('\t'*level, 'query_word: ', query_word)
        # print('\t'*level, 'node.token: ', node.token)
        # print('\t'*level, 'partial_result: ', partial_result)
        # print('\t'*level, 'num_errors: ', num_errors)
        # print('\t'*level, 'max_errors: ', max_errors, '\n')
        
        if num_errors > max_errors:
            # print('\t'*level, 'Too many errors\n')
            return
        
        if node.is_end:
            # print('\t'*level, 'Found end.\n')
            self.result[partial_result] = {
                'trie_object': node.trie_object,
                'num_errors': num_errors
            }
            return
            
        if len(query_word) < 1:
            # print('\t'*level, 'Empty word.\n')
            return

        query_token = query_word[0]
        for token in node.children:
            # print('\t'*level, 'Level: ', level, 'token: ', token, '\n')
            if token != query_token:
                error_increment = 1
            else:
                error_increment = 0
            self._full_query_dfs_with_errors(
                node = node.children[token],
                query_word = query_word[1:],
                partial_result = partial_result + token,
                num_errors = num_errors + error_increment,
                max_errors = max_errors,
                level = level + 1)

    def full_query_with_errors(self, word, max_errors):
        """Given an input (a word), retrieve all words that match that word from the trie,
        with max_errors allowed
        """
        
        # print('word: ', word, '\n')
        
        query_token = word[0]
        num_errors = 0
        self.result = {}
        level = 0
        for token in self.root.children:
            # print('Level: ', level, 'token: ', token, '\n')
            if token != query_token:
                error_increment = 1
            else:
                error_increment = 0
            self._full_query_dfs_with_errors(
                node = self.root.children[token],
                query_word = word[1:],
                partial_result = token,
                num_errors = num_errors + error_increment,
                max_errors = max_errors,
                level = level + 1)
        return self.result
