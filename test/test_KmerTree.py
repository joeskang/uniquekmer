import sys
import os
sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')

from KmerTree import KmerTree

def test_kt():
    kt = KmerTree()
    kt.add("ACGCGT")
    kt.add("CGGCGA")
    kt.add("TCGCGA")
    kt.add("ACTCTCTATCTGCGCGTCG")
    print(kt)


if __name__ == '__main__':
    test_kt()
