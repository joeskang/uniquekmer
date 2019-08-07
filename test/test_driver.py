"""
    test the different components of driver

    part 0: read genome and s_array files
    part 1: compute lcp array with kasai algorithm
    part 2: find groups of unambiguous and ambiguous chars in genome (unam/ambs)
    part 3: with provided length and distance, use unam/ambs to find valid start addresses of k-mers

"""

import sys
sys.path.append('../src')
import driver
#import argparse
import time
import numpy as np
from collections import OrderedDict, deque
import file_manager as fm
from SA import naive_SA

class Args:

    genome = str
    SA = str
    outfile= str
    length = int
    low = int
    quiet = bool
    distance = int
    uniquesfile = str
    ambs=str


    def __init__(self, genome='', SA='', outfile='', length=100, low=20, distance=300, quiet=False, uniquesfile='', ambs=''):
        self.genome=genome
        self.SA=SA
        self.outfile=outfile
        self.length=length
        self.low=low
        self.distance=distance
        self.quiet=quiet
        self.uniquesfile=uniquesfile
        self.ambs=ambs






def main():

    #args = read_args()
    args = Args()

    args.genome = fm.append_file_name('test/fake_genome')

    genome, past, s_array, start = _test_part_0(args)

    past, sa_uniques = _test_part_1(genome=genome, s_array=s_array, args=args, past=past)

    a_u_dict, past, unam = _test_part_2(args=args, past=past)

    _test_part_3(args=args, past=past, sa_uniques=sa_uniques, unam=unam, a_u_dict=a_u_dict)

def _test_part_0(args:Args):

    if not args.SA:
        sequence = fm.reads(args.genome)
        s_array, L = naive_SA(string=sequence)
        fm.write_array_to_byte(byte_arr = s_array, filename=fm.append_file_name('test/fake_SA'))
        args.SA = fm.append_file_name('test/fake_SA')

    genome, past, s_array, start = driver._part_0(args=args)

    # genome is a string of the sequence
    assert type(genome) is str and genome

    # s_array is a numpy array
    assert type(s_array) is np.ndarray

    # start is a time.time object
    assert type(start) is float

    # past is also a time.time object
    assert type(past) is float

    return genome, past, s_array, start

def _test_part_1(genome:str, s_array, args, past=time.time):

    past, sa_uniques = driver._part_1(genome=genome, past=past, s_array=s_array, args=args)

    # past is a time.time object
    assert type(past) is float

    # sa_uniques is OrderedDict
    assert type(sa_uniques) is OrderedDict and sa_uniques

    return past, sa_uniques

def _test_part_2(args, past=time.time):

    a_u_dict, past, unam = driver._part_2(args=args, past=past)

    # a_u_dict is dict
    assert type(a_u_dict) is dict and a_u_dict

    # past is a time.time object
    assert type(past) is float

    # unam is a deque
    assert unam is deque and unam

    return a_u_dict, past, unam

def _test_part_3(args, sa_uniques, unam, a_u_dict, past=time.time):

    return

"""
def read_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('-g', dest='genome')
    parser.add_argument('-s', dest='SA')
    parser.add_argument('-o', dest='outfile')
    parser.add_argument('-m', dest='length')
    parser.add_argument('-l', dest='low')
    parser.add_argument('-q', dest='quiet')
    parser.add_argument('-d', dest='distance')

    args = parser.parse_args()

    args.length = 100 if not args.length else int(args.length)
    args.low = 20 if not args.low else int(args.low)
    if not args.genome:
        args.genome = fm.append_file_name('test/fake_genome')

    args.distance = 300 if not args.distance else int(args.distance)

    args.quiet = False if (args.quiet == 'f' or args.quiet == 'F' or args.quiet == 'false' or args.quiet == 'False') else True

    return args
"""

if __name__ == '__main__':
    args = Args()
    args.genome = fm.append_file_name('data/genome.fa')
    past = time.time()
    driver._part_2(past=past, args=args)

