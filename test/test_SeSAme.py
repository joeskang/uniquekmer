import sys
import os

sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')
from file_manager import append_file_name, read_byte_numpy, read_unambiguous
sys.path.append(append_file_name('experimental'))
from SeSAme import SAGuide

def sesame_seed():
    se = SAGuide()
    test_str = "A" * 10

    a = se._next_string(test_str)
    print(a)

    test_str = "ATTTTT"
    print(se._next_string(test_str))

    print(se._next_string("TTTTT"))


def sesame_plant():
    s_arr = read_byte_numpy(append_file_name('data/22.sa'))
    geno = read_unambiguous(append_file_name('data/22.fa'))

    se = SAGuide(s_arr=s_arr, geno=geno)
    print(se)



if __name__ == "__main__":
    sesame_plant()
