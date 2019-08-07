from collections import deque, OrderedDict
import numpy as np
import sys
sys.path.insert(0,'../src')
#from rb_tree import RedBlackTree
from exception_manager import NotSameLengthError
from tqdm import tqdm


def sort_dict(in_dict: dict):
    '''
        simply sorts the keys and returns new dict with sorted keys

    :param in_dict:
    :return:
    '''
    return_d = OrderedDict()
    sorted_keys = sorted(in_dict.keys())
    for key in tqdm(sorted_keys, desc="Sorting dictionary: "):
        return_d[key] = in_dict[key]

    return return_d


def ambs_unam_dict(ambs: deque, unam: deque):
    """
        return dictionary with unam as key, ambs as value
    :param ambs:
    :param unam:
    :return:

    ******** DEPRECATED ************
    (ambs_unam_split in ambiguous_split should handle this)
    """
    ambs_unam_d = OrderedDict()
    try:
        if len(ambs) != len(unam):
            raise NotSameLengthError

        else:
            while True:
                ambs_unam_d[unam.popleft()] = ambs.popleft()

    except IndexError:
        pass
    except NotSameLengthError:
        print("The length of the AMBS and UNAM deques were different. aborting")

    finally:
        return ambs_unam_d


def sort_dict_to_numpy(in_dict: dict):
    '''
        returns sorted keys and corresponding values
            returns as a numpy array for plotting

    :param in_dict:
    :return:
    '''

    sorted_values = deque()
    sorted_keys = sorted(in_dict.keys())
    for key in sorted_keys:
        sorted_values.append(in_dict[key])

    return np.array(sorted_keys), np.array(sorted_values)

'''
def set_rb_tree(in_deque: deque):
    rb = RedBlackTree()

    try:
        while True:
            rb.add(in_deque.popleft())
    except IndexError as e:
        raise
    finally:
        return rb
'''

def count_dict(in_dict: dict, top: int):
    """
        returns a dict with:
            - keys: unique values from in_dict
            - values: counts of the values in in_dict

    :param in_dict:
    :return count:
    """

    myq = deque(list(in_dict.values()))
    count = {}
    # keys = lengths of lcps
    # values = counts of lcps
    while myq:
        key = myq.pop()
        count[key] = count[key] + 1 if key in count else 1
    return count


def set_Ordered_sa_uniques(key_d: deque, value_d: deque, top:int)->OrderedDict:
    """
        unique_start: the starting k length from a given starting address from s_array
            where uniqueness begins
        take sa and lcp arrays and return a dict with keys: sa, values: unique_starts
        also get valid ranges


        consuming too much memory
        TODO: deprecate and use other functions
    :param key_d:
    :param value_d:
    :return:
    """

    d = OrderedDict()
    past = 0

    try:
        # check and see that there are same number of lcp values as sa addresses
        if len(key_d) != len(value_d):
            raise NotSameLengthError
        else:
            runs = len(key_d)
        pbar = tqdm(total=runs, desc="Creating OrderedDict with s_array as key and valid unique start address as values")
        # start unique_start check
        while key_d and value_d:
            current_lcp = value_d.popleft()

            pbar.update(1)
            # ensure that unique_start is within top limit
            if current_lcp > top or past > top:
                past = current_lcp
                continue

            # unique_start depends will be the bigger of the current or past lcp values
            # unique_start = past_lcp if past_lcp > current_lcp else current_lcp
            d[key_d.popleft()] = past if past > current_lcp else current_lcp
            past = current_lcp



        # when returning, key_d and value_d will be empty
        return d

    except NotSameLengthError as e:
        print(e)
    except IndexError:
        print("Index Error caught at set_Ordered_sa_uniques()")
        return d


def count_kmer(in_dict: dict, bot: int, top: int):
    """
        returns a dict that counts all the unique kmers
        in_dict = sa_lcp_dict
            keys = sa address
            values = lcp values

    :param in_dict:
    :param top:
    :return:
    """

    myq = deque((list(in_dict.values())))
    count = {}
    # keys = lengths of kmers
    # values = counts of kmers

    # loop through the dictionary and increase corresponding lcp values
    while myq:
        value = myq.pop()
        # if the value is anything less than the minimum limit, then set value to the minimum
        if value < bot:
            value = bot
        # if the count is already in dictionary, increase by 1, else set to 1
        count[value] = 1 if value not in count else count[value] + 1

    # initialize all the empty spots in count dictionary that did not have corresponding lcp values
    for i in range(top - bot + 1):
        if bot + i not in count:
            count[bot + i] = 0

    # without having to loop through each time the lcp values were set in the first while loop,
    # run for loop and just increase sequentially
    sorted_keys = sorted(count)
    counter = 0
    for key in sorted_keys:
        counter += count[key]
        count[key] = counter

    return count

    """ 
        # naive implementation:
        # go from lowest to top, counting all unique kmers
        for i in range(top-value + 1):

            count[value + i] = 1 if value + i not in count else count[value + i] + 1
            
    """


if __name__ == '__main__':
    # testing sort_dict()
    mydict = {5: 'ab', 1: 'ij', 8: 'ih', 4: 'qu'}
    keys, values = sort_dict_to_numpy(mydict)
    print(keys)
    print(values)

    # test count_kmer
    ckmer = {123: 6, 124: 7, 125: 8, 126: 3, 127: 6}
    top = 10
    bot = 4
    count = count_kmer(in_dict=ckmer, bot=bot, top=top)
    print("result of count_kmer = ", sorted(count.items()))
