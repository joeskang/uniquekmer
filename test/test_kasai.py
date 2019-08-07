import sys
import os
from tqdm import trange, tqdm
sys.path.append(os.getcwd().split('uniquekmer')[0] + '/uniquekmer/test')
sys.path.append(os.getcwd().split('uniquekmer')[0] + '/uniquekmer/src')
from SA import naive_SA
from kasai import kasai, inverse1
import numpy as np
import file_manager as fm
import os
import time
from unittest import TestCase


def main_lcp():
    loud = False

    eman_proc = fm.append_file_name('test/emancipation_proclamation')

    if os.path.isfile(eman_proc):
        kasai, naive = find_sa_lcp(eman_proc, loud=loud)
        print('passed emancipation proclamation lcp comparison')
        print('kasai mean time: ', np.mean(kasai), ' with std dev: ', np.std(kasai))
        print('naive mean time: ', np.mean(naive), ' with std dev: ', np.std(naive))

    share_wealth = fm.append_file_name('test/share_our_wealth')

    if os.path.isfile(share_wealth):
        kasai, naive = find_sa_lcp(share_wealth, loud=loud)
        print('passed share our wealth lcp comparison')
        print('kasai mean time: ', np.mean(kasai), ' with std dev: ', np.std(kasai))
        print('naive mean time: ', np.mean(naive), ' with std dev: ', np.std(naive))

    fake_geno = fm.append_file_name('test/fake_genome')

    if not os.path.isfile(fake_geno):
        fm.fake_genome(filename=fake_geno, lines=1000, alpha={'N': 0, 'A': 1, 'C': 2, 'G': 3, 'T': 4})

    kasai, naive = find_sa_lcp(filename=fake_geno, loud=loud)
    print('passed fake genome lcp comparison')
    print('kasai mean time: ', np.mean(kasai), ' with std dev: ', np.std(kasai))
    print('naive mean time: ', np.mean(naive), ' with std dev: ', np.std(naive))


def naive_lcp(s_array, T:str):
    # naive approach to calculating lcp
    size = len(s_array)
    lcp_array = [None for _ in tqdm(s_array,desc='creating empty list for lcp')]
    lcp = 0
    for i in trange(size, desc='computing lcp naively: '):
        index = s_array[i]
        if i == size - 1:
            lcp_array[i] = 0
            continue

        suffix = T[index:]

        next_index = s_array[i + 1]
        next_suffix = T[next_index:]
        lcp = 0

        for j in range(len(suffix)):
            if suffix[j] != next_suffix[j]:
                break
            lcp += 1

        lcp_array[i] = lcp

    lcp_array[size - 1] = 0
    return lcp_array


def find_sa_lcp(filename: str, loud: bool):
    kasai_ = list()
    naive_ = list()
    with open(filename, 'r') as file:

        for str in file:
            if str:
                sa_array, L = naive_SA(str,verbose=loud)
                size = len(sa_array)
                sa_array = np.asarray(sa_array)

                # print("Suffix Array : \n", sa_array)
                past = time.time()
                kasai_lcp = kasai(str, sa_array, verbose=loud)
                kasai_.append(time.time() - past)


                if loud:
                    print("Kasai LCP Array : \n", kasai_lcp)
                past = time.time()
                naive_lcp_array = naive_lcp(T=str, s_array=sa_array)
                naive_.append(time.time() - past)

                if loud:
                    print("Naive LCP Array : \n", naive_lcp_array)

                assert kasai_lcp and naive_lcp_array
                assert kasai_lcp == naive_lcp_array # Todo: why is this broken?
    return kasai_, naive_


def find_best_inverse_macro(filename: str, loud: bool):
    ut = TestCase()
    sequence = fm.simple_string_read(filename=filename)
    sa_array, L = naive_SA(sequence, verbose=loud)

    sa_array = np.asarray(sa_array)
    past = time.time()

    inv1 = inverse1(sa_array)

    time_i_1 = time.time() - past

    past = time.time()

    inv2 = inverse2(sa_array)

    time_i_2 = time.time() - past

    assert type(inv1) is type(inv2)

    ut.assertCountEqual(first=inv1, second=inv2)

    return time_i_1, time_i_2


def find_best_inverse_micro(filename: str, loud: bool):
    ut = TestCase()

    inv1_time = list()
    inv2_time = list()
    with open(filename, 'r') as file:

        for str in file:
            if str:
                sa_array, L = naive_SA(str, verbose=loud)
                size = len(sa_array)

                sa_array = np.asarray(sa_array)

                past = time.time()

                inv1 = inverse1(sa_array)

                inv1_time.append(time.time() - past)

                past = time.time()

                inv2 = inverse2(sa_array)

                inv2_time.append(time.time() - past)

                assert type(inv1) is type(inv2)

                ut.assertCountEqual(first=inv1, second=inv2)

    return inv1_time, inv2_time


def inverse2(sa_array):

    size = len(sa_array)
    inv_suff = [None for _ in range(size)]
    for i in range(size):
        inv_suff[sa_array[i]] = i

    return np.asarray(inv_suff)


def main_inverse():
    loud = False

    inv1_ = list()
    inv2_ = list()

    eman_proc = fm.append_file_name('test/emancipation_proclamation')

    if os.path.isfile(eman_proc):
        inv1_, inv2_ = find_best_inverse_micro(eman_proc, loud=loud)
        print('passed emancipation proclamation line-by-line lcp comparison\n')
        print('mean of inverse1 method: ', np.mean(inv1_), ' with std dev: ', np.std(inv1_))
        print('\nmean of inverse2 method: ', np.mean(inv2_), ' with std dev: ', np.std(inv2_))
        print()

        inv1_, inv2_ = find_best_inverse_macro(eman_proc, loud=loud)

        print('passed emancipation proclamation whole lcp comparison\n')
        print('time of inverse1 method: ', inv1_)
        print('\ntime of inverse2 method: ', inv2_)
        print()

    fake_geno = fm.append_file_name('test/fake_genome')

    if not os.path.isfile(fake_geno):
        fm.fake_genome(filename=fake_geno, lines=2000, alpha={0:'N', 1:'A', 2:'C', 3:'G', 4:'T'})

    inv1_, inv2_ = find_best_inverse_micro(fake_geno, loud=loud)
    print('passed fake genome line-by-line lcp comparison\n')
    print('mean of inverse1 method: ', np.mean(inv1_), ' with std dev: ', np.std(inv1_))
    print('\nmean of inverse2 method: ', np.mean(inv2_), ' with std dev: ', np.std(inv2_))
    print()

    inv1_, inv2_ = find_best_inverse_macro(fake_geno, loud=loud)
    print('passed fake genome whole lcp comparison\n')
    print('time of inverse1 method: ', np.mean(inv1_))
    print('\ntime of inverse2 method: ', np.mean(inv2_))
    print()

    print('End.')


# test lcp already made
def kasai_truth_test():
    print("reading suffix array:")
    s_arr = fm.read_byte_numpy(filename=fm.append_file_name('data/22.sa'))
    print("suffix array read.\nreading genome:")
    genome = fm.reads(filename=fm.append_file_name('data/22.fa'))
    print("genome read.\nsplitting genome:")
    splits = genome.split('N')
    gen = ''
    for part in splits:
        gen += part

    print("genome split.\nreading lcp file:")

    lcp = fm.unjson_it(filename=fm.append_file_name('output/0408_c22_1552_json_lcp'))

    past = ''
    current = ''
    past_sa = 0

    print('lcp file read.\ntesting uniqueness')

    assert len(s_arr) == len(lcp)

    for i in trange(len(s_arr)):
        sa = s_arr[i]
        l = lcp[i]
        current = gen[sa:sa+l+1]
        if past == current:
            raise(Exception("Not unique sequences found: "+past+" at address " + str(past_sa)+" and "+current+" at address " + str(sa)))

        else:
            past = current
            past_sa = sa


def test_inverse1():

    dummy_suffix_array = [5,3,1,2,4,0]
    inv = inverse1(dummy_suffix_array)
    inv = inv.tolist()
    actual_inverse = [5,2,3,1,4,0]
    assert inv == actual_inverse


if __name__ == '__main__':
    #main_lcp()
    #kasai_truth_test()
    test_inverse1()
