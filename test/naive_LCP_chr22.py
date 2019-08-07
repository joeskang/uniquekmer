"""
    get naive lcp of chr22

"""
import sys
import os

filepath = os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src'
sys.path.append(filepath)

print(filepath)

import file_manager as fm
import test_kasai

def naive_lcp_22():
    s_arr = fm.read_byte_numpy(filename=fm.append_file_name('data/22.sa'))
    lcp = test_kasai.naive_lcp(s_array=s_arr,T=simple_genome())
    fm.json_it(data=lcp,filename=fm.append_file_name('output/naive_lcp_22'))

def simple_genome():
    seq = ''
    with open(fm.append_file_name('data/22.fa'),'r') as file:
        for line in file:
            if line[0] != '>':
                seq += line.split()[0]
    return seq

if __name__ == '__main__':
    naive_lcp_22()
