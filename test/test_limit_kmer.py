# test limit_kmer.py modules with pytest
from collections import deque
import random
import sys
import os
sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')
from limit_kmer import  within_distance
from exception_manager import SamePositionError

def test_valid_address():
    try:
        # parameters for valid_address: unam(rb_tree), sa_lcp_dict(dict), ambs_unam_dict(dict), top(int)

        # fill lcp with random numbers from 0-100, with replacement
        lcp = random.choices(population=range(100),k=1000)

        # fill s_array with random numbers from 0-100, without replacement
        s_array = deque(random.sample(population=range(1000), k=1000))


        # set top as 10 for now
        top = 10

        # set sa_lcp_dict
        d2 = dict()
        # setting array with random 20 ints 1~100. sorted
        myarr = sorted(random.sample(range(1,100),20))
        for entry in myarr:
            d2[entry] = random.randint(1,10)


        unam = deque(myarr)

        result = valid_address(unam=unam, lcp=lcp, au_dict=d2, top=top, s_arr=s_array)

        # check is not None
        assert result is not None

        # check if dict
        assert type(result) is dict

        # check is not empty dict
        #assert len(result) is not 0

    except SamePositionError:
        # with random values, to be expected, try again
        pass




def test_within_distance():
    # make random parameters to check very basic requirements
    # parameters: in_dict, top, distance
    d = dict()
    myarr = sorted(random.sample(range(1,100),20))

    for item in myarr:
        d[item] = random.randint(1,10)


    # set top as 10
    top = 10

    # set distance as 3
    distance = 3


    result = within_distance(in_dict=d, top=top, distance=distance)

    # check if not none
    assert result is not None

    # check if dict
    assert type(result) is dict

    # check if not empty dict
    assert len(result) is not 0
