"""
    begin anew
"""
from collections import deque, OrderedDict
import re
from limit_kmer import sort_sa_lcp
from tqdm import tqdm, trange
from rb_tree import RedBlackTree


def unique_starts(lcp:list):
    past = 0
    runs = len(lcp)
    for l in tqdm(lcp, desc="calculating unique starts: "):
        yield l if l > past else l
        past = l


def just_true_address(au:RedBlackTree, sorted_lcp: list):
    assert type(au) == RedBlackTree

    au_dict = {}
    for item in list(au):
        au_dict[item[0]] = item[1]

    del au

    u_ceil = list(au_dict)[0]
    u_floor = 0

    a_offset = au_dict[u_ceil]
    for sa in trange(len(sorted_lcp), desc="Calculating True Addresses: "):
        if sa < u_floor:
            raise Exception("SA is less than u_floor. Possible that s_arr not sorted correctly?")

        if sa > u_ceil and len(au_dict) > 1:

            u_floor = u_ceil
            del au_dict[u_ceil]
            u_ceil = list(au_dict)[0]
            a_offset = au_dict[u_ceil]

        elif len(au_dict) <= 1:
            print("not au_dict reached")
            break

        yield (sa + a_offset)



def true_address_with_sort(lcp: list, au: RedBlackTree, inv_suff:list, top=100, bot=20, distance=300):
    # as suffix array will be sorted anyways, unnecessary to retrieve suffix array
    # dict to be returned

    assert type(au) == RedBlackTree
    # realize that sa becomes sorted, use simple list instead of rb_tree
    au_dict = {}
    for item in list(au):
        au_dict[item[0]] = item[1]
    del au

    # no need for s_arr to be deque, use list
    if type(lcp) != list:
        lcp = list(lcp)
    assert len(lcp) == len(inv_suff)
    true_addresses = []
    past = 0
    tops = []

    lcp = list(unique_starts(lcp))

    #######
    first_counter = 0
    second_counter = 0

    # TODO: sa is sorted, is it necessary to use RedBlackTree?(04/11/2019)
    u_ceil = list(au_dict)[0]
    u_floor = 0
    a_offset = au_dict[u_ceil]
    sorted_lcp = []
    for sa in trange(len(inv_suff), desc='Calculating True Addresses: '):
        # a sorted Suffix Array should be a list from 0 --> length - 1, thus use i instead of retrieving from list

        # get next suffix address and corresponding lcp (unique start) value
        lcp_curr = int(lcp[int(inv_suff[sa])])

        if sa < u_floor:
            raise Exception("SA is less than u_floor. Possible that s_arr not sorted correctly?")

        if sa > u_ceil and len(au_dict) > 1:

            # myl = list from RedBlackTree
            # [unam ambs]

            # delete and save as much memory
            u_floor = u_ceil
            del au_dict[u_ceil]
            u_ceil = list(au_dict)[0]
            a_offset = au_dict[u_ceil]

        elif len(au_dict) <= 1:
            print("not au_dict reached")
            break

        if (lcp_curr + sa) > u_ceil:
            first_counter +=1
            continue
        elif lcp_curr > top:
            second_counter += 1
            continue

        #difference = u_ceil - sa # do we need to add 1 here?

        # try yielding instead
        # true_addresses.append(sa + a_offset)
        # tops.append(difference if difference < top else top)
        # sorted_lcp.append(lcp_curr)
        yield (sa+a_offset, lcp_curr)



def true_address_no_sort(lcp: list, au: RedBlackTree, s_arr: deque, top=100, bot=20, distance=300):

    # dict to be returned
    d = OrderedDict()

    assert type(au) == RedBlackTree

    if type(lcp) != deque:
        lcp = deque(lcp)

    if type(s_arr) != deque:
        s_arr = deque(s_arr)

    assert len(lcp) == len(s_arr)
    d = {}
    past = 0
    runs = len(lcp)
    for _ in trange(runs,desc="calculating unique starts: "):
        l = lcp.popleft()
        lcp.append(past if past > l else l)
        past = l
    runs = len(s_arr)
    pbar = tqdm(total=runs, desc="Finding true addresses: ")

    u_ceil = 0
    u_floor = 0
    myl = []
    while s_arr:
        sa = s_arr.popleft()
        pbar.update(1)
        # get next suffix address and corresponding lcp (unique start) value
        lcp_curr = lcp.popleft()

        if type(sa) != int and type(lcp) != int:
            sa = int(sa)
            lcp_curr = int(lcp_curr)

        if sa > u_ceil or sa < u_floor:

            # myl = list from RedBlackTree
            # [unam ambs]
            myl = au.floor([sa])
            u_floor = 0 if not myl else myl[0]

            myl = au.ceil([sa])
            if not myl:
                continue
            u_ceil = myl[0]

        if (lcp_curr + sa) > u_ceil:
            continue

        else:
            difference = u_ceil - sa # do we need to add 1 here?
# todo: is the u_floor going too low?
        d[sa + myl[1]] = (lcp_curr, (difference if difference < top else top))

    return d


def true_address_dict(indict:dict, au:RedBlackTree):

    # ! ensure that indict values are sorted numerically
    keys = deque(indict.keys())
    values = deque(indict.values())
    del indict

    assert type(au) == RedBlackTree
    # realize that sa becomes sorted, use simple list instead of rb_tree
    au_dict = {}
    for item in list(au):
        au_dict[item[0]] = item[1]
    del au

    # TODO: sa is sorted, is it necessary to use RedBlackTree?(04/11/2019)
    u_ceil = list(au_dict)[0]
    u_floor = 0
    a_offset = au_dict[u_ceil]

    pbar = tqdm(total=len(keys), desc="Calculating True Addresses: ")
    # for _ in trange(len(keys), desc='Calculating True Addresses: '):
    while keys and values:
        pbar.update(1)
        # a sorted Suffix Array should be a list from 0 --> length - 1, thus use i instead of retrieving from list

        # get next suffix address and corresponding lcp (unique start) value
        adr = keys.popleft()
        end = values.popleft()

        if adr < u_floor:
            raise Exception("SA is less than u_floor. Possible that s_arr not sorted correctly? SA: ", adr, " u_floor ", u_floor)

        if adr > u_ceil and len(au_dict) > 1:

            # myl = list from RedBlackTree
            # [unam ambs]

            # delete and save as much memory
            u_floor = u_ceil
            del au_dict[u_ceil]
            u_ceil = list(au_dict)[0]
            a_offset = au_dict[u_ceil]

        elif len(au_dict) < 1:
            print("not au_dict reached")
            break

        if end > u_ceil:
            # first_counter += 1
            continue

        # try yielding instead
        # true_addresses.append(sa + a_offset)
        # tops.append(difference if difference < top else top)
        # sorted_lcp.append(lcp_curr)
        yield adr + a_offset, end + a_offset


def true_address(lcp: list, au_dict: dict, top: int, s_arr: deque, bot:int, distance:int):

    # dict to be returned
    d = {}

    a_offset = 0

    assert au_dict

    unam = deque(au_dict.keys())

    # first entry is 0 for unam
    # pop off and increase a_offset
    u_offset = unam.popleft()
    a_past = au_dict[u_offset]
    a_offset = 0
    u_offset -= 1

    if type(s_arr) != deque:
        s_arr = deque(s_arr)
    if type(lcp) != deque:
        lcp = deque(lcp)

    lcp, s_arr = sort_sa_lcp(lcp=lcp, s_arr=s_arr)

    while s_arr and lcp:
        # get next suffix address and corresponding lcp (unique start) value
        sa = s_arr.popleft()
        lcp_curr = lcp.popleft()

        if not (type(sa) == int and  type(u_offset) == int and type(lcp) == int):
            sa = int(sa)
            u_offset = int(u_offset)
            lcp_curr = int(lcp_curr)

        if sa > u_offset:
            u_curr = unam.popleft()
            a_offset += a_past
            a_past = au_dict[u_curr]
            u_offset += u_curr

        if (lcp_curr + sa) > u_offset:
            continue

        else:
            difference = u_offset - sa + 1

        d[sa + a_offset] = (lcp_curr, (difference if difference < top else top))

    return d


def forbiddens(inv_suff:list, lcp:list, au):

    """
        as the number of forbidden addresses are actually relatively low compared
        to the total number of trues, instead of converting and storing all trues
        return list of forbidden addresses instead.

        mainly for use in MU with forwards backwards
    :param invsuff:
    :param lcp:
    :return:
    """

    try:
        assert type(au) == RedBlackTree
        # realize that sa becomes sorted, use simple list instead of rb_tree
        au_dict = {}
        for item in list(au):
            au_dict[item[0]] = item[1]
        del au

        # no need for s_arr to be deque, use list
        if type(lcp) != list:
            lcp = list(lcp)
        assert len(lcp) == len(inv_suff)
        true_addresses = []
        past = 0
        tops = []

        lcp = list(unique_starts(lcp))

        #######
        first_counter = 0
        second_counter = 0

        u_ceil = list(au_dict)[0]
        u_floor = 0
        a_offset = au_dict[u_ceil]
        sorted_lcp = []
        for sa in trange(len(inv_suff), desc='Calculating Forbidden Addresses: '):
            # a sorted Suffix Array should be a list from 0 --> length - 1, thus use i instead of retrieving from list

            # get next suffix address and corresponding lcp (unique start) value
            lcp_curr = int(lcp[int(inv_suff[sa])])

            if sa < u_floor:
                raise Exception("SA is less than u_floor. Possible that s_arr not sorted correctly?")

            if sa > u_ceil and len(au_dict) > 1:

                # myl = list from RedBlackTree
                # [unam ambs]

                # delete and save as much memory
                u_floor = u_ceil
                del au_dict[u_ceil]
                u_ceil = list(au_dict)[0]
                a_offset = au_dict[u_ceil]

            elif len(au_dict) <= 1:
                print("not au_dict reached")
                break

            if (lcp_curr + sa) > u_ceil:
                first_counter +=1
                yield sa

    except Exception as e:
        print(e)
        breakpoint()
