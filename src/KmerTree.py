from experimental.LCSHash import Block


class Node:

    def __init__(self, parent, value=None, a=None, c=None, g=None, t=None, block=None):

        self.value = value
        self.a = a
        self.c = c
        self.g = g
        self.t = t
        self.parent = parent
        self.block = block

    def __iter__(self):

        if self.a:
            yield from self.a.__iter__()

        if self.c:
            yield from self.c.__iter__()

        if self.g:
            yield from self.g.__iter__()

        if self.t:
            yield from self.t.__iter__()

    def add_block(self, block):
        if self.block:
            raise Exception("The node already has a block")
        else:
            self.block = block

    def append(self, acc):
        if not self.block:
            self.block = Block(acc)
        else:
            self.block.append(acc)

    def has_block(self)->bool:
        return bool(self.block)

    def has_value(self)->bool:
        return bool(self.value)

    def has_children(self)->bool:
        return bool(self.get_children_count())

    def get_children_count(self) -> int:
        if not (self.a or self.c or self.g or self.t):
            return 0
        return sum([int(bool(self.a)), int(bool(self.c)), int(bool(self.g)), int(bool(self.t))])


class KmerTree:

    def __init__(self):
        self.count = 0
        self.root = Node(parent=None)
        self.kmer_count = 0

    def __iter__(self):
        yield from self.root.__iter__()

    def __repr__(self):
        return "Number of Nodes: {count}, k-mers: {kmer}".format(count=self.count, kmer=self.kmer_count)

    def add(self, kmer:str):
        try:
            self.kmer_count += 1

            cs = kmer[1:-1]
            acc = kmer[0] + kmer[-1]
            parent, node, val = self._find_parent(cs)
            if node is None:
                # assume that the node that is returned has a value, if not then we're adding to a non-MU
                # TODO: after fixing MU, assert that nodes have value
                if parent.has_value():
                    compare()
                assert parent.has_block
                parent.append(acc)
                return
            new_node = Node(parent=parent, block=Block(acc))
            if node == 'A':
                parent.a = new_node
            elif node == 'C':
                parent.c = new_node
            elif node == 'G':
                parent.g = new_node
            elif node == 'T':
                parent.t = new_node
            else:
                raise Exception("Received " + str(node) +  " as char")

            self.count += 1

        except AssertionError as e:
            raise e

    # def _compare(self, node_a, node_b, parent):

    def _find_parent(self, kmer):

        try:
            assert len(kmer) > 2

            def _inner(parent, i, kmer):
                if not parent:
                    raise AssertionError("Parent is nonetype. Kmer: " + kmer)

                if i > len(kmer) - 1:
                    return parent, None, None

                elif kmer[i] == 'A':
                    if not parent.a:
                        return parent, 'A', kmer[i:] if i >= len(kmer) - 1 else None
                    return _inner(parent.a, i+1, kmer)

                elif kmer[i] == 'C':
                    if not parent.c:
                        return parent, 'C', kmer[i:] if i >= len(kmer) - 1 else None
                    return _inner(parent.c, i+1, kmer)

                elif kmer[i] == 'G':
                    if not parent.g:
                        return parent, 'G', kmer[i:] if i >= len(kmer) - 1 else None
                    return _inner(parent.g, i+1, kmer)

                elif kmer[i] == 'T':
                    if not parent.t:
                        return parent, 'T', kmer[i:] if i >= len(kmer) - 1 else None
                    return _inner(parent.t, i+1, kmer)

                else:
                    raise Exception("unknown char: " + str(kmer[i]) + " found.")

            return _inner(self.root, 0, kmer)

        except AssertionError:
            raise

        except Exception:
            raise

