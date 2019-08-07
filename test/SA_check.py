import sys
import os
from tqdm import trange,tqdm

sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')
from file_manager import append_file_name, read_byte_numpy, reads, read_byte_array, read_unambiguous


def main():
    print('Reading SA: ')
    s_array = read_byte_array(append_file_name('data/22.sa'))
    print('SA read!\nReading genome: ')
    # !! reads returns ambs
    genome = reads(filename=append_file_name('data/22.fa'))
    gen_list = genome.split('N')
    genome = ''
    for part in gen_list:
        genome += part
    #genome = read_unambiguous(filename=append_file_name('data/22.fa'))
    print('Genome read!')

    length = 30
    s_len = len(s_array)
    for i in trange(s_len, desc='Checking validity of suffix array: '):
        sa = s_array[i]
        if sa+length+1 < s_len:
            s0 = genome[sa:sa+length+1]
            s1 = genome[s_array[i+1]:s_array[i+1]+length+1]

            print('s0: ',s0)
            print('s1: ',s1)

            assert s0 <= s1

        else:
            pass


if __name__ == '__main__':
    main()
