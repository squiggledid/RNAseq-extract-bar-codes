import sys

project_dir = "/Users/davids/git/WEHI_CoViD_RNAseq/RNAseq-extract-bar-codes"
sys.path.append(project_dir)

from ID_Trie_Node import ID_Trie_Node # TODO learn the proper Python way to do this

class ID_Trie(object):
    """The trie object"""

    def __init__(self):
        """
        The trie has at least the root node.
        The root node does not store any character
        """
        self.root = ID_Trie_Node("")
    
    def insert(self, word):
        """Insert a word into the trie"""
        node = self.root
        
        # Loop through each character in the word
        # Check if there is no child containing the character, create a new child for the current node
        for char in word:
            if char in node.children:
                node = node.children[char]
            else:
                # If a character is not found,
                # create a new node in the trie
                new_node = ID_Trie_Node(char)
                node.children[char] = new_node
                node = new_node
        
        # Mark the end of a word
        node.is_end = True

        # Increment the counter to indicate that we have seen this word once more
        node.counter += 1
        
    def depth_first_search(self, node, prefix):
        """Depth-first traversal of the trie
        
        Args:
            - node: the node to start with
            - prefix: the current prefix, for tracing a
                word while traversing the trie
        """
        if node.is_end:
            self.output.append((prefix + node.char, node.counter))
        
        for child in node.children.values():
            self.depth_first_search(child, prefix + node.char)
        
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
        for char in prefix:
            if char in node.children:
                node = node.children[char]
            else:
                # cannot find the prefix, return empty list
                return []
        
        # Traverse the trie to get all candidates
        self.depth_first_search(node, prefix[:-1])

        # Sort the results in reverse order and return
        return sorted(self.output, key=lambda prefix: prefix[1], reverse=True)

    def full_query_dfs(self, node, word, partial_result):
        """Do a depth first search until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        """
            
        print('word: ', word)
        print('node.char: ', node.char)
        print('partial_result: ', partial_result, '\n')

        if node.is_end:
            print('Found end.\n')
            return (partial_result + node.char, node.counter) # Success!
        if len(word) < 1:
            print('Empty word.\n')
            return[] # ran out of characters before finding a match

        char = word[0]
        if char in node.children:
            return self.full_query_dfs(node.children[char], word[1:], partial_result + node.char)
        else:
            # print('Not in trie.\n')
            return [] # if is not in trie

    def full_query_dfs(self, node, word, partial_result):
        """Do depth first search until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        """
            
        # print('word: ', word)
        # print('node.char: ', node.char)
        # print('partial_result: ', partial_result)

        if node.is_end:
            # print('Found end.\n')
            return (partial_result + node.char, node.counter) # Success!
        if len(word) < 1:
            # print('Empty word.\n')
            return[] # ran out of characters before finding a match

        char = word[0]
        if char in node.children:
            return self.full_query_dfs(node.children[char], word[1:], partial_result + node.char)
        else:
            # print('Not in trie.\n')
            return [] # if is not in trie

    def full_query(self, word):
        """Given an input (a word), retrieve that word from the trie 
        """
        
        # print('word: ', word)
        
        node = self.root
        char = word[0]
        if char in node.children:
            node = node.children[char]
        else:
            # cannot find the first character of id, return empty list
            return []
        
        # Traverse the trie to find the id
        return self.full_query_dfs(node, word[1:], '')
        
    def full_query_dfs_with_errors(self, node, query_word, partial_result, num_errors, max_errors, level):
        """Do depth first searches until we either reach an end node (meaning
        have a match), or cannot find the next node (meaning the initial query
        word is not in the trie)
        """
        # print('\t'*level, 'query_word: ', query_word)
        # print('\t'*level, 'node.char: ', node.char)
        # print('\t'*level, 'partial_result: ', partial_result)
        # print('\t'*level, 'num_errors: ', num_errors)
        # print('\t'*level, 'max_errors: ', max_errors, '\n')
        
        if num_errors > max_errors:
            # print('\t'*level, 'Too many errors\n')
            return
        
        if node.is_end:
            # print('\t'*level, 'Found end.\n')
            self.result[partial_result] = {
                'count': node.counter,
                'num_errors': num_errors
            }
            return
            
        if len(query_word) < 1:
            # print('\t'*level, 'Empty word.\n')
            return

        query_char = query_word[0]
        for char in node.children:
            # print('\t'*level, 'Level: ', level, 'char: ', char, '\n')
            if char != query_char:
                error_increment = 1
            else:
                error_increment = 0
            self.full_query_dfs_with_errors(
                node = node.children[char],
                query_word = query_word[1:],
                partial_result = partial_result + char,
                num_errors = num_errors + error_increment,
                max_errors = max_errors,
                level = level + 1)

    def full_query_with_errors(self, word, max_errors):
        """Given an input (a word), retrieve all words that match that word from the trie,
        with max_errors allowed
        """
        
        # print('word: ', word, '\n')
        
        query_char = word[0]
        num_errors = 0
        self.result = {}
        level = 0
        for char in self.root.children:
            # print('Level: ', level, 'char: ', char, '\n')
            if char != query_char:
                error_increment = 1
            else:
                error_increment = 0
            self.full_query_dfs_with_errors(
                node = self.root.children[char],
                query_word = word[1:],
                partial_result = char,
                num_errors = num_errors + error_increment,
                max_errors = max_errors,
                level = level + 1)
        return self.result