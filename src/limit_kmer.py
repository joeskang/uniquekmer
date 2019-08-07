"""
    inputs: SA, lcp array, bt2 index
    purpose: uses bt2 index to find prohibited address for
        start of unique kmers

    'ambs' = ambiguous character (ie. 'N')
    'unam' = unambiguous character (ie. 'A', 'G', 'C', 'T')

    functions:
            ambs_and_unam: creates corresponding ambs and unam files
            within_range: only returns sa,lcp pairs that meet the top and bot requirements for lengths
            __within_distance: only returns sa,lcp pairs that meet the distance requirement
"""
from collections import deque, defaultdict, OrderedDict
import sys
import os
from tqdm import tqdm, trange
from exception_manager import print_memory_usage

sys.path.append(os.getcwd().split('uniquekmer')[0] + '/uniquekmer/src')
from exception_manager import InvalidStart, SamePositionError


def sort_sa_lcp(lcp:deque, s_arr:deque):
    try:
        # +++ 04-08-2019:
        # choosing to go without unique_starts
        # sort sa and lcp according to sorted sa
        assert len(lcp) == len(s_arr)

        if type(lcp) != deque:
            lcp = deque(lcp)
        if type(s_arr) != deque:
            s_arr = deque(s_arr)

        #d = {}
        d = [0] * len(s_arr)
        past = 0
        runs = len(lcp)
        #pbar = tqdm(total=runs, desc="Creating Dict to Sort")

        while lcp and s_arr:
            #pbar.update(1)
            l = lcp.popleft()
            d[s_arr.popleft()] = past if past > l else l
            past = l

        print("sorting suffix array")
        s_arr = deque(sorted(d.keys()))

        for sa in tqdm(s_arr,desc="reconstructing lcp array"):
            lcp.append(d[sa])

        return lcp, s_arr

    except MemoryError:
        print_memory_usage()
        raise

    except Exception as e:
        raise e

def within_distance(trues: list, top: int, distance=300):
    trues = deque(trues)
    past = 0
    flag = False
    try:
        for _ in trange(len(trues), desc='Finding starting positions within '+str(distance)+' bp: '):
            if past:
                curr = trues.popleft()
                if curr <= past + distance:
                    if not flag:
                        trues.append(past)

                    flag = True
                    trues.append(curr)

                else:
                    flag = False

                past = curr

    except Exception:
        raise Exception



def within_distance(in_dict: dict, top: int, distance: int):
    """
        s = start of current k-mer
        e = end of current k-mer

        n = start of next k-mer
        m = end of next k-mer

        y = top + sa + distance

    """

    print("Keys of dictionary sorted. Keys: ", len(in_dict.keys()))
    keys = deque()
    if type(in_dict) == dict:
        keys = deque(sorted(in_dict))
    else:
        keys = deque(in_dict)

    assert type(top) == int

    assert type(distance) == int

    # all the lcps remaining in the dict are <= distance after within_range()
    # TODO: assumes that distance > top. write for case distance < top
    try:
        pbar = tqdm(total=len(keys),desc="Finding k-mer starts within "+str(distance)+" bp: ")

        # initialize
        e = 0  # the end point of the top-mer + distance
        past = 0  # the last start pos that was seen

        while keys:
            pbar.update(1)
            sa = int(keys.popleft())

            assert type(sa) == int

            assert type(sa) == int

            # if the start of the next address is greater than end of current + distance
            #   yet is also less than the
            if e < sa:
                if sa < past + distance + top:  # if the sa is still less than end of past top-mer + distance
                    # set the new e to this top-mer starting at sa + distance
                    e = sa + top + distance
                # else, check if the distance between current sa and next start pos is greater than distance
                elif keys[0] > sa + top + distance:
                    # if so, delete the key,value pair from the dictionary
                    del in_dict[sa]

            past = sa

        return in_dict

    except IndexError as e:
        return in_dict

    except Exception as e:
        raise e


if __name__ == '__main__':
    print("Please do not run this file directly.")
