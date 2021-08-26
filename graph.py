"""
Find connected components.
"""


class Node:
    def __init__(self, value, parent):
        self.value = value
        self.parent = parent

    def __repr__(self):
        return "Node(value={}, parent={})".format(self.value, self.parent)


class ComponentFinder:

    def __init__(self, values):
        self.nodes = {x: Node(x, None) for x in values}

    def merge(self, x, y):
        assert x != y
        x_root = self._find_node(x)
        y_root = self._find_node(y)

        if x_root is y_root:
            return

        # Merge while making sure that the node with the smaller value is the
        # new parent.
        if x_root.value < y_root.value:
            y_root.parent = x_root
        else:
            x_root.parent = y_root

    def _find_node(self, x):
        node = root = self.nodes[x]
        while root.parent is not None:
            root = root.parent

        # compression path
        while node.parent is not None:
            node.parent, node = root, node.parent
        return root

    def find(self, x):
        """
        Return which component x belongs to, identified by the smallest value.
        """
        return self._find_node(x).value

    def print(self):
        for x in sorted(self.nodes):
            print(x, ":", self.nodes[x], "is represented by", self._find_node(x))
