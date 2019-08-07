"""
    LCSTree:
        Hash table that has the lowest common substring (combination of lowest common prefix and suffix) as keys
        and has an array of chars (i = 0 and i = n)

        hopefully this can reduce the memory usage by about a 16th in the best case

"""
from collections import deque


class Block:

    def __init__(self, val:str = None):
        self.vals = [val] if val else []

    def __eq__(self, blk):
        return self.vals == blk.vals

    def __repr__(self):
        yield from self.vals

    def __iter__(self):
        return self.vals

    def append(self, val:str):
        assert len(val) <= 2
        if val in self.vals:
            raise Exception("K-mer already exists in Block")
        self.vals.append(val)
        self.vals.sort()

    def exists(self, val:str):
        return val in self.vals

    def get_different(self, blk):
        if not self == blk:

            diff = set(self.vals).symmetric_difference(set(blk.values))
            for char in diff:

                if 'Z' in char:
                    if char[0] == "Z":
                        yield char[1]
                    elif char[1] == "Z":
                        yield char[0]
                    else:
                        raise Exception("An unknown char pair found in Block.values: ", char)

                else:
                    yield char[0], char[1]

    def values(self):
        return self.vals


class LCSHash:
    """
        base is a dictionary where keys: common substring, values = Block
        in case that a SNP-inserted substring is unique, change values to (Block, address of SNP)

    """

    def __init__(self):
        self.base = {}
        self.kmer_count = 0

    def __len__(self):
        return len(self.base)

    def __repr__(self):
        return "Number of CS: {length}, k-mers: {kmer} ".format(length=len(self.base), kmer=self.kmer_count)

    def insert(self, kmer:str, is_snp=False):
        # ensure no ambs
        assert 'N' not in kmer
        cs = kmer[1:-1]
        acc = kmer[0] + kmer[-1]

        if cs not in self.base:
            if is_snp:
                self.base[cs] = (Block(), [acc])
            else:
                self.base[cs] = Block(acc)

        elif type(self.base[cs]) == Block and not self.base[cs].exists(acc):
            if is_snp:
                self.base[cs] = (self.base[cs], [acc])
            else:
                self.base[cs].append(acc)

        elif type(self.base[cs]) == tuple and not self.base[cs][0].exists(acc) and not acc in self.base[cs][1]:
            if is_snp:
                self.base[cs][1].append(acc)
            else:
                self.base[cs][0].append(acc)

        else:
            return

        self.kmer_count += 1

    def contains(self, kmer:str)->bool:
        """
            True: The SNP-inserted k-mer is not unique
            False: The SNP-inserted k-mer is unique
        :param kmer:
        :return:
        """
        # cs = common substring
        cs = kmer[1:-1]
        # acc = accessories (i.e. first and last chars)
        acc = kmer[0] + kmer[-1]

        if cs in self.base:
            # first check if a SNP has been inserted for that cs:
            if type(self.base[cs]) == Block:

                # kmers = list(self.base[kmer[1:-1]].get_different())
                #if self.base[cs].exists(acc):
                    # check if the kmer is not unique
                #    return True # SNP-inserted kmer is not unique
                #else:
                    # SNP-inserted kmer is unique
                    # self.base[cs].append(acc)
                    # self.base[cs] = (self.base[cs], [acc])
                #    return False

                return True if self.base[cs].exists(acc) else False

            # else, unique SNP-mer has been included
            elif type(self.base[cs]) == tuple:
                #if acc in self.base[cs][1]:

                #    return True
                #elif self.base[cs][0].exists(acc):
                #    return True
                #else:
                    # self.base[cs][1].append(acc)

                return True if (acc in self.base[cs][1] or self.base[cs][0].exists(acc)) else False

            else:
                raise Exception("Unknown type found for base[cs]")

        # else, cs was not in self.base. is unique
        else:
            # self.base[cs] = (Block(), acc)
            return False


def save_LCSTree(lcs:LCSHash, filename:str):

    base = lcs.base
    with open(filename, 'w') as file:

        for key in base:
            file.write('{ CS: ' + key + ' acc: [')
            vals = base[key].values() if type(base[key]) == Block else base[key][0].values()
            for v in vals:
                file.write(' ' + v)

            file.write(' ]')
            if type(base[key]) == tuple:
                file.write(' SNP: [')
                snp = base[key][1]
                for s in snp:
                    file.write(' ' + s)
                file.write(' ]')

            file.write(' } ')




