from sys import path
from tqdm import tqdm, trange
import os
import numpy as np
import traceback
from collections import deque, OrderedDict
from ambiguous_split import split_sequence, rb_tree_ambs
from file_manager import read_byte_numpy, reads, msgunpack_dict, append_file_name, json_it, unjson_it, chr_splits, write_array_to_byte, read_byte_array, read_unambiguous
from true_address import true_address_with_sort, true_address_no_sort, true_address_dict, forbiddens, just_true_address

# TODO: will hard code directory paths for now 05/21
if "bioinformatics" in os.getcwd():
    path.append('/work/bioinformatics/s187520/minimal_unique/src')
else:
    path.append(os.getcwd().split('Documents')[0] + 'PycharmProjects/minimal_uniquemer/src')

from MU_internal import mu, compare, sort_mu
from contiguous import find_contigs, contig_seqs, find_perfect_contigs, within_distance
PATH = append_file_name('data/')
#ORIGINAL = '22.fa'
#FLIPPED = 'flip22'
ORIGINAL = 'genome.fa'
FLIPPED = 'flippedGeno'

bt2_SA1 = '22'
bt2_SA2 = 'flip22'

files = {
    'MU_RESULT':'/output/c22_mu_result',
    'PERFECT_CONTIGS':'/output/c22_perfect_contigs',
    'PERFECT_CONTIGS_WITH_DISTANCE':'/output/c22_perfect_contigs_with_distance'

}


def mu_driver():

    """
        similar function as driver.py, except include minimal uniques instead of finding 20-100 uniquemers
    :return:
    """

    try:

        # gitignore()
        print('reading original genome: ', end='')
        chrs, geno = chr_splits(filename=PATH + ORIGINAL)
        json_it(chrs, append_file_name("json_chrs"))
        del chrs
        print('done.\nreading original SA...: ', end='')
        s_arr = read_byte_numpy(append_file_name('data/genome.sa'))

        lcp1 = kasai(geno, s_arr)[1]
        d1 = OrderedDict(mu(SA=s_arr, LCP=lcp1))
        del lcp1
        del s_arr

        au = _part_2(genome_file_name=PATH + ORIGINAL)
        print("au list: ", list(au))

        # *************************
        # (2) flipped
        # *************************

        print("performing flips: ")
        geno2 = read_unambiguous(PATH + FLIPPED)

        s_arr2 = read_byte_numpy(append_file_name('data/flippedGeno.sa'))

        lcp2 = kasai(geno2, s_arr2)[1]
        del geno2

        mu_result = dict(compare(d=d1, SA=s_arr2, LCP=lcp2))
        del lcp2
        mu_result = OrderedDict(sort_mu(mu_result))

        mu_result = OrderedDict(true_address_dict(mu_result, au))

        json_it(mu_result, append_file_name(files['MU_RESULT']))

        #contigs = list(find_contigs(d=old_mu_result_without_true_addresses, bot=20, top=100))
        contigs = OrderedDict(find_perfect_contigs(d=mu_result, bot=20, top=100))

        json_it(contigs, append_file_name(files['PERFECT_CONTIGS']))

        contigs = list(within_distance(d=contigs, distance=300))

        json_it(contigs, append_file_name(files['PERFECT_CONTIGS_WITH_DISTANCE']))

        print("number of contigs: ", len(contigs))

        print("done")

    except Exception as e:
        raise


if 'bioinformatics' in os.getcwd():
    path.append(os.getcwd().split('uniquekmer')[0] + 'minimal_unique/src')
else:
    path.append(os.getcwd().split('Documents')[0] + 'PycharmProjects/minimal_uniquemer/src')


from file_manager import chr_splits, append_file_name, read_byte_numpy, read_unambiguous
from kasai import kasai

PATH = append_file_name('data/')
FORWARD = '22.fa'
FLIPPED = 'flip22'


COMMENTS = deque([
    'READING ORIGINAL GENOME: ',
    'DONE.\nREADING FORWARD SA: ',
    'AU LIST: ',
    'PERFORMING FLIPS: ',
    'NUMBER OF STITCHED: ',
    'NUMBER OF CONTIGS: ',
    'DONE'
])


def efficient_mu_driver():
    """
        NOTES:
            07/05: You MUST run get_uniques first before sorting the lcp

    :return:
    """


    try:
        # comment()
        geno = reads(filename=PATH + FORWARD)
        geno_length = len(geno)
        # comment()
        s_arr = read_byte_numpy(append_file_name('data/22.sa'))

        inv_suff, lcp = kasai(geno, s_arr)
        lcp = kasai(geno, s_arr)[1]
        del geno, s_arr


        # comment()
        au = _part_2(genome_file_name=PATH + FORWARD)

        lcp = list(get_uniques(lcp))

        trues0 = list(sort_lcp(lcp=lcp, inv_suff=inv_suff))

        del lcp

        bad_address = forbiddens(inv_suff=inv_suff, lcp=trues0, au=au)
        del inv_suff

        geno = read_unambiguous(filename=PATH + FLIPPED)
        s_arr = read_byte_numpy(append_file_name('data/f22.sa'))

        inv_2, lcp = kasai(geno, s_arr)
        lcp = kasai(geno, s_arr)[1]
        del geno, s_arr


        lcp = list(get_uniques(lcp))

        trues1 = list(sort_lcp(lcp=lcp, inv_suff=inv_2))
        del lcp, inv_2


        # mu_s, mu_e = list(compare(inv0=inv_suff, trues0=trues0, inv1=inv_2, trues1=trues, bad_address=bad_address))
        mu_s = []
        mu_e = []

        au_dict = {}
        for item in list(au):
            au_dict[item[0]] = item[1]

        del au

        u_ceil = list(au_dict)[0]
        u_floor = 0
        a_offset = au_dict[u_ceil]

        # mu_s, mu_e = list(compare(trues0=trues0, trues1=trues1, bad_address=bad_address))
        for tup in compare_no_inv_suff(trues0=trues0, trues1=trues1, bad_address=bad_address, geno_length=geno_length):
            sa = tup[0]

            if sa < u_floor:
                raise Exception("SA is less than u_floor. Possible that s_arr not sorted correctly?")

            if sa > u_ceil and len(au_dict) > 1:
                u_floor = u_ceil
                del au_dict[u_ceil]
                u_ceil = list(au_dict)[0]
                a_offset = au_dict[u_ceil]

            elif len(au_dict) < 1:
                print("not au_dict reached")
                break

            # mu_s.append(tup[0])
            mu_s.append(sa + a_offset)
            mu_e.append(tup[1])

        # TODO: 07/05 made the line below return a dict as well as accept geno
        #   to return to before, do not input geno and output two lists
        # myd = dict(compare(trues0 = trues0, trues1 = trues1, bad_address=bad_address))
        # json_it(data=myd, filename="c22_mu")
        assert len(mu_s) == len(mu_e)
        just_dump(myl=mu_s, fn="c22_mu_starts_0709", print_length=True)
        just_dump(myl=mu_e, fn="c22_mu_ends_0709", print_length=True)
        # 07/08: changed get_uniques so that it doesn't yield past or lcp + 1


        # json_it(mu_s, "efficient_mu_starts")
        # json_it(mu_e, "efficient_mu_ends")

        # stitched = list(stitch(starts=mu_s, uniques=mu_e))
        # json_it(stitched, "stitched")

        # print("Number of stitched: " + str(len(stitched)))

        # print("Number of MU: " + str(len(mu_s)))
        # findmean(mys = mu_s, mye=mu_e)

    except IndexError:
        pass
    except Exception:
        print(traceback.format_exc())
        breakpoint()


def temp_forward_unique_check():

    geno = read_unambiguous(filename=PATH + FORWARD)
    # comment()
    s_arr = read_byte_numpy(append_file_name('data/22.sa'))

    inv_suff, lcp = kasai(geno, s_arr)

    myd = {}
    for num in range(len(s_arr)):
        myd[s_arr[num]] = lcp[num]

    #trues0 = list(get_uniques(lcp))
    json_it(trues0, "c22_forward_uniques")


def sort_lcp(lcp:list, inv_suff:list):
    assert len(lcp) == len(inv_suff)
    for s in trange(len(lcp), desc="sorting lcp"):
        yield lcp[inv_suff[s]]


def compare(inv0: list, trues0: list, inv1: list, trues1: list, bad_address: list):
    """
        match (s,e) to (e,s)
    :param inv1:
    :param trues1:
    :param inv0:
    :param trues0:
    :param bad_address:
    :return:
    """

    try:
        mys = []
        mye = []
        forbidden_count = 0
        assert len(inv1) == len(inv0) == len(trues1) == len(trues0)
        length = len(inv0)
        for s0 in trange(length, desc="finding minimal uniques: "):

            if s0 in bad_address:
                forbidden_count += 1
                continue

            # forward: start, end (s0, e0)
            # reverse: end, start (e1, s1)
            # if lcp_0 == lcp_1
            unique_0 = trues0[inv0[s0]]
            e0 = s0 + unique_0
            s1 = length - 1 - e0
            unique_1 = trues1[inv1[s1]]

            # breakpoint()
            if unique_0 == unique_1:
                mys.append(s0)
                mye.append(e0)


        print("number of forbidden addresses: " + str(forbidden_count))

        return mys, mye

    except Exception as e:
        print(traceback.format_exc())
        breakpoint()


def compare_with_starts(true_start_0: list, true_start_1: list, trues0: list, trues1: list, bad_address: list, geno_length: int):
    try:
        forbidden_count = 0

        length = len(true_start_0)
        for _ in trange(length, desc="finding minimal uniques: "):
            s0 = true_start_0[_]

            if s0 in bad_address:
                forbidden_count += 1
                continue

            unique_0 = trues0[_]
            e0 = s0 + unique_0
            s1 = geno_length - 1 - e0
            unique_1 = trues1[s1]
            breakpoint()

            if unique_0 == unique_1:
                yield (s0, e0)

        print("number of forbidden addresses: " + str(forbidden_count))

    except Exception as e:
        print(traceback.format_exc())
        breakpoint()


def compare_no_inv_suff(trues0: list, trues1: list, bad_address: list, geno_length: int):
    """
        match (s,e) to (e,s)
    :param inv1:
    :param trues1:
    :param inv0:
    :param trues0:
    :param bad_address:
    :return:
    """

    try:
        forbidden_count = 0
        assert len(trues1) == len(trues0)
        length = len(trues0)
        for s0 in trange(length, desc="finding minimal uniques: "):

            if s0 in bad_address:
                forbidden_count += 1
                continue

            # forward: start, end (s0, e0)
            # reverse: end, start (e1, s1)
            # if lcp_0 == lcp_1
            unique_0 = trues0[s0]
            e0 = s0 + unique_0
            s1 = length - 1 - e0
            unique_1 = trues1[s1]

            # breakpoint()
            if unique_0 == unique_1:
                # assert(forward[s0:e0] == flipped[s1:s1+unique_1][::-1])
                # mys.append(s0)
                # mye.append(e0)
                yield (s0, unique_0)

        print("number of forbidden addresses: " + str(forbidden_count))
        # return mys, mye

    except Exception as e:
        print(traceback.format_exc())
        breakpoint()

def temp():
    c22 = unjson_it(filename="c22_mu")
    l0 = list(c22.keys())
    l1 = list(c22.values())
    just_dump(l0, l1, fn = "c22_mu_just_dump")


def just_dump(l0:list,l1:list, fn:str):
    # just dump a list into the file
    assert len(l0) == len(l1)
    with open(fn, "w") as file:
        for _ in trange(len(l0)):
            if _ != 0:
                file.write(", ")
            file.write(str(l0[_]) + " " + str(l1[_]))

    return


def just_dump (myl:list, fn:str, print_length: bool = False):
    with open(fn, "w") as file:

        for _ in range(len(myl)):
            if _ == 0 and print_length:
                file.write("length: " + str(len(myl)) + "\n")
            file.write(str(myl[_]) + '\n')




def findmean(mys, mye):
    myl = []
    assert len(mys) == len(mye)
    for _ in range(len(mys)):
        myl.append(mye[_] - mys[_] + 1)

    myl = np.asarray(myl)
    print("mean: ", np.mean(myl))
    print("std dev: ", np.std(myl))
    print("median: ", np.median(myl))


def comment():
    try:
        print(COMMENTS.popleft())
    except IndexError:
        return


def _part_2(genome_file_name: str):

    try:
        # _________________________________
        print('\n_____________________________________')
        print('PART 2: FIND UNAM AND AMBS')
        print('_____________________________________\n')
        # _________________________________
        au = None

        ambs, unam = split_sequence(filename=genome_file_name)
        au = rb_tree_ambs(ambs, unam)

        assert au
        return au

    except Exception as e:
        breakpoint()


def gitignore():
    with open('../.gitignore','a') as file:
        for f in files:
            if not os.path.exists(f):
                file.write('\n')
                file.write(f)


def get_uniques(lcp):
    past = 0
    for l in tqdm(lcp, desc="getting uniques: "):
        yield l+1 if l > past else past+1
        past = l


def stitch(starts:list, uniques:list):
    assert len(starts) == len(uniques)
    cluster_start = starts[0]
    past_unique = uniques[0]

    for _ in trange(len(uniques), desc="stitching MU's: "):
        start = starts[_]
        un = uniques[_]

        if start > past_unique:
            yield cluster_start, past_unique
            if _ == len(uniques) - 1:
                return start, un
            cluster_start = start

        elif _ == len(uniques) - 1:
            return cluster_start, un

        past_unique = un


if __name__ == '__main__':
    efficient_mu_driver()
    #temp_forward_unique_check()
    # temp()
