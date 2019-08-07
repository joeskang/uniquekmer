"""
    compute lcp array from suffix array based on the Kasai Algorithm
"""

import sys
import os
import time
import numpy as np
from tqdm import trange, tqdm

sys.path.append(os.getcwd().split('uniquekmer')[0]+'/uniquekmer/src')
from exception_manager import get_time


def kasai(genome, s_array, inv_suff=None, verbose=False, print=print) -> tuple:
    """
        returns an unordered dictionary with SA as keys and LCP as values
    :param genome:
    :param s_array:
    :return:


    consumes too much memory using dict for sa_lcp_dict
    """

    past = time.time()

    if verbose:
        print('Truncating sequence to A/C/G/T')

    # (1): truncate genome to only A/C/G/T
    if 'N' in genome:
        print('N found in genome for kasai(). Parsing for unambiguous only')

        splits = genome.split('N')
        genome = ''
        for part in splits:
            genome += part

    size = len(s_array)
    gen_size = len(genome)

    if verbose:
        print('Truncation done.')
        past = get_time(past,print)
        print("SA size: %s", str(size))
        print("genome size: %s", str(gen_size))

    lcp_arr = [0] * size

    if verbose:
        print("Empty List created for LCP")
        print("Creating Empty Lists for Inverse Suffix Array")

    if not inv_suff:
        inv_suff = inverse1(s_array=s_array)
        inv_suff = inv_suff.astype(int)

    inv_suff = inv_suff.tolist()

    if verbose:
        # ____________________________________ #
        print("Inverse Suffix Array Completed")
        past = get_time(past, print)
        # ____________________________________ #

    # length of previous lcp
    k = 0

    s_array = s_array.tolist()

    for i in trange(size, desc="Computing LCP"):

        sa = inv_suff[i]
        if sa == size - 1:
            k = 0
            # sa_lcp_dict[sa] = k
            continue

        if type(sa) != int:
            sa = int(sa)

        j = s_array[sa + 1]

        while i + k < gen_size and j + k < gen_size and genome[i + k] == genome[j + k]:
            k += 1
        # sa_lcp_dict[sa] = k
        lcp_arr[sa] = k

        if k > 0:
            k -= 1

    if verbose:
        print("LCP Array Completed")
        get_time(past)

    # return sa_lcp_dict
    return inv_suff, lcp_arr


def inverse1(s_array, verbose=False):

    size = len(s_array)
    b = np.array([_ for _ in trange(size, desc='Create empty list for inverse suffix array')])

    inv_suff = np.asarray([0] * size)

    inv_suff[s_array] = b

    return inv_suff


def inverse2(s_array, verbose=False):
    size = len(s_array)
    inv_suff = [0] * size

    for i in range(size):
        inv_suff[s_array[i]] = i

    return inv_suff


def inverse3(s_array):
    size = len(s_array)
    from collections import deque
    s_array = deque(s_array)
    inv_suff = [0] * size

    for i in range(size):
        inv_suff[s_array.popleft()] = i

    return inv_suff

def naive_SA(string: str, verbose=False):
    """
		inputs a string and outputs deque for the suffix array and string
	:param string:
	:return:
	"""

    str_dict = {}
    sa_array = deque()

    # generate suffixes
    for i in trange(len(string), desc='Creating suffixes: '):
        temp = string[i:]

        # store the position
        str_dict[temp] = i

    # generate deque with actual suffixes
    suffixes = deque(sorted(str_dict.keys()))

    # generate sa_array with position values to corresponding suffixes
    # also, generate last column, L
    L = ""
    pbar = tqdm(total=len(suffixes), desc='Generating suffix array: ')
    while suffixes:
        index = str_dict[suffixes.popleft()]
        sa_array.append(index)
        L += string[index - 1]
        pbar.update(1)

    if verbose:
        print("Suffix array: ", sa_array)
        print("The last column is: ", L)

    return sa_array, L

if __name__ == '__main__':
    from collections import deque
    str = 'banana'

    verbose = True
    sa_array, L = naive_SA(str, verbose=False)


    size = len(sa_array)
    sa_array = np.asarray(sa_array).ravel()
    print(type(sa_array))
    if verbose:
        print("Suffix Array : \n", sa_array)

    inv, lcp = kasai(str, sa_array)

    if verbose:
        print("LCP Array : \n", lcp)

    d = {}
    for i in range(len(sa_array)):
        d[sa_array[i]] = lcp[i]
