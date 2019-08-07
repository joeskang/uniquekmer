import sys
import os

from tqdm import tqdm

sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')
from file_manager import append_file_name, read_unambiguous, read_byte_numpy
from kasai import kasai

TWO_WAY_GENO = append_file_name('data/2_way_genome.fa')
S_ARRAY = append_file_name('data/2_way_genome.sa')

# todo: does building a 2-way genome and computing SA and LCP intermingle the forward and reverse?


class OddSuffixArrayLength(Exception):
    """the length of forward + reverse suffix array is not even"""
    raise Exception("Length of suffix array is not even!")


def driver():

    geno = read_unambiguous(TWO_WAY_GENO)
    s_arr = read_byte_numpy(S_ARRAY)


def mu_2way_internal(geno:str, s_arr):
    try:
        # check if the length of s_arr is even. (forward length + reverse length)
        if not len(s_arr) % 2 == 0:
            raise OddSuffixArrayLength

        length = int(len(s_arr)/2)

        inv_suff, lcp = kasai(genome=geno,s_array=s_arr)
        uniques = list(get_uniques(lcp))


    except OddSuffixArrayLength as e:
        raise e


def get_uniques(lcp):
    """
        as kasai builds lcp array of the 2-way genome in the forward direction,
        no need to reverse the flow in which the unique is computed for the reverse
        simply apply "yield l if l > past else past" to the whole 2-way

    :param lcp:
    :return:
    """
    past = 0
    for l in tqdm(lcp, desc="getting uniques: "):
        yield l if l > past else past
        past = l

def compare()

if __name__ == "__main__":
    driver()
