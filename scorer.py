def nice_print(debug, *args):
    if debug:
        print(*args)

class ScoreVector:
    def __init__(self, x, y, z):
        self.multiloops = x
        self.unpaired = y
        self.branches = z
    
    def __iter__(self):
        yield self.multiloops
        yield self.unpaired
        yield self.branches

    def __str__(self) -> str:
        return str(tuple(self))

    def __add__(self, other):
        return ScoreVector(*(x1+x2 for x1, x2 in zip(self, other)))


class IntervalTreeNode:
    """
    Node in an interval search tree

    Each node is an interval in some range, with the requirement
    that each node's interval is a strict subinterval of its parent's
    """
    def __init__(self, start=0, end=0):
        self.start = start
        self.end = end
        self.children = []

    def insert(self, start, end, sort=False):
        """
        Insert a new node, either as a child here or by handing off to a child
        """
        if start <= self.start or end >= self.end:
            raise ValueError("Could not insert node")

        if start >= end:
            raise ValueError("Invalid interval")

        # Check whether [start, end] is a subinterval of any child and, if so, insert it there
        for child in self.children:
            if child.start < start and end < child.end:
                child.insert(start, end, sort=sort)
                return

        # If the interval was not contained in a child, we insert it here
        child = IntervalTreeNode(start, end)
        self.children.append(child)

        if sort:
            self.sort()

    def __lt__(self, b):
        return self.start < b.start or (self.start == b.start and self.end < b.end)

    def sort(self):
        """
        Sort the children in increasing order of start base position
        """
        self.children.sort()

        # Recurse! Recurse!
        for child in self.children:
            child.sort()

    def valency(self):
        """
        Return the number of children of this node
        """
        return len(self.children)

class RNAStructure():
    def __init__(self, structure, debug=False):
        self.debug = debug

        self.structure = structure
        self.root = IntervalTreeNode(-1, len(structure))
        pairs = RNAStructure.get_pairs(structure)[::-1]
        for pair in pairs:
            self.root.insert(pair[0], pair[1], sort=True)
    
    @classmethod
    def get_pairs(self, structure):
        open_helices = []
        pairs = []
        for i, c in enumerate(structure):
            if c == "(":
                open_helices.append(i)
            if c == ")":
                start = open_helices.pop()
                pairs.append((start, i))
        
        if len(open_helices):
            raise Exception(f"Helices not closed, still have {len(open_helices)} remaining.") 
        
        return pairs

    def score_structure(self):
        nice_print(self.debug, "Starting structure scoring.")
        score = ScoreVector(0,0,0)

        for child in self.root.children:
            score += self.score_internal_node_recursively(child)

        return score

    def score_internal_node_recursively(self, node):
        score = ScoreVector(0,0,0)

        i = node.start
        j = node.end
        if node.valency() >= 2:
            score += self.multiloop(i, j, node, score)

        for child in node.children:
            score += self.score_internal_node_recursively(child)

        return score
    
    def multiloop(self, i, j, node, score):
        loop = self.scoreM(node)
        nice_print(self.debug, f"Multiloop ({i}, {j})")
        return loop

    def scoreM(self, node):
        score = ScoreVector(0,0,0)
        score.multiloops += 1
        score.branches += (node.valency() + 1)

        if node.valency() > 0:
            score += self.scoreMUnpairedRegion(node.start, node.end, node.children[0].start, node.children[0].end)
            score += self.scoreMUnpairedRegion(node.children[-1].start, node.children[-1].end, node.start, node.end)
            for child, neighbor in zip(node.children, node.children[1:]):
                score += self.scoreMUnpairedRegion(child.start, child.end, neighbor.start, neighbor.end)

        return score

    def scoreMUnpairedRegion(self, i1, j1, i2, j2):
        start, end = self.determine_unpaired_region(i1, j1, i2, j2)
        score = ScoreVector(0,(end-start-1),0)

        return score

    def determine_unpaired_region(self, i1, j1, i2, j2):
        if i1 < i2 and i2 < j2 and j2 < j1:
            # (i1, j1) is the initiating pair of the loop, so
            # region is between i1 and i2
            start, end = i1, i2
        elif i2 < i1 and i1 < j1 and j1 < j2:
            # (i2, j2) is the initiating pair of the loop, so
            # region is between j1 and j2
            start, end = j1, j2
        elif i1 < j1 and j1 < i2 and i2 < j2:
            # neither (i1, j1) and (i2, j2) is the initiating pair, so
            # region is between j1 and i2
            start, end = j1, i2
        else:
            raise ValueError("Invalid multiloop nesting")

        return start, end

if __name__ == "__main__":
    # s1 = ScoreVector(1, 2, 3)
    # s2 = ScoreVector(4,3,7)
    # s3 = s1+s2
    # print(s3)
    # s1.multiloops = 4
    # print(s3)

    struct = RNAStructure("((...)(...)).(...)(...).", debug=False)
    print(struct.score_structure())
    # for i in range(200000):
    #     struct.score_structure()