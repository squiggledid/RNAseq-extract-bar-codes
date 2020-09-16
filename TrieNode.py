class TrieNode:
    """A node in the trie structure"""

    def __init__(self, token):
        # the token represented by this node in the path through the trie. Could be a character, could be a byte
        self.token = token

        # flag indicating that this is a leaf node
        self.is_end = False

        # An object name, optional stored in the leaf nodes of the treioe
        self.trie_object = None

        # a dictionary of child nodes
        # keys are tokens; values are nodes
        self.children = {}
