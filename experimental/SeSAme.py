"""
    Open se SA me
        SA guide to quick indexing of strings
        let guide store up to 4^16 (4,294,967,296) indices
"""

alpha = {
    'A':0,
    'C':1,
    'G':2,
    'T':3
}

from tqdm import tqdm
from collections import deque

class SAGuide:

    def __len__(self):
        return len(self.guide)

    def __repr__(self):
        return '{guide}'.format(guide=self.guide)

    def __init__(self, max=13, s_arr=None, geno=None):

        self.max = max
        self.guide = [0 for _ in range(4**max)]
        if not (s_arr is None or geno is None):
            self.guide = list(self._guide_make(s_arr=s_arr, geno=geno))

    def _guide_make(self, s_arr, geno:str):
        end = len(s_arr)

        # assumes that text is unam DNA chars (A, C, G, T)
        # will run 4^max times

        compare_str = 'A' * self.max
        end_str = 'T' * self.max

        # runs = 4**self.max
        runs = len(s_arr)
        s_arr = deque(s_arr)

        pbar = tqdm(total=runs, desc="building guide")

        while s_arr:
            pbar.update(1)
            s = s_arr.popleft()
            if s + self.max > end:
                continue

            elif geno[s:s+self.max] == compare_str:
                yield s
                if compare_str == end_str:
                    return

            elif geno[s:s+self.max] > compare_str:

                while geno[s:s+self.max] > compare_str:
                    compare_str = self._next_string(compare_str)
                    yield None

                if geno[s:self.max] == compare_str:
                    yield s

                    if compare_str == end_str:
                        return

            elif geno[s:s+self.max] < compare_str:
                continue

            compare_str = self._next_string(compare_str)


    def _next_string(self, instr:str):
        assert len(instr) > 1

        def _inner_next(instr:str):
            if len(instr) == 1 and instr == 'T':
                return 'AA'

            if instr[-1] == 'T':
                return _inner_next(instr[:-1]) + 'A'

            elif instr[-1] == 'A':
                return instr[:-1] + 'C'

            elif instr[-1] == 'C':
                return instr[:-1] + 'G'

            elif instr[-1] == 'G':
                return instr[:-1] + 'T'

        return _inner_next(instr)


    def get_index(self, instr:str):
        """
            from inputted string, compute the corresponding index in the guide that the prefix should reside
        :param instr: inpuutted string
        :return:
        """

        if len(instr) < self.max:
            raise Exception("inputted string shorter than the SA guide entries")

        index = 0
        for i in range(self.max):
            index += alpha[instr[-(i+1)]]





