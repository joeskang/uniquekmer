import struct
import numpy as np
import array
import random
import msgpack
import Bio.Seq
import time
from sys import path
import os
path.append(os.getcwd().split('uniquekmer')[0]+'/uniquekmer/src')
from file_manager import append_file_name, read_byte_array, read_unambiguous, read_byte_numpy, reads, msgunpack_dict, chr_splits, flip_text


def test_chr_splits():
    # ASSUMES THAT reads() HAS BEEN THOROUGHLY TESTED AND VALID
    s0 = reads(filename='../data/22.fa')
    chrs, s1 = chr_splits('../data/22.fa')

    assert s0 == s1


def flip_check():

    # (1) checkpoint 1
    dummy_text = "banana"
    new_text = flip_text(intext=dummy_text)
    assert new_text == "ananab"

    # (2) checkpoint 2
    flip_text(intext=dummy_text, outfile="flip_check_dummy_text")

    # (3) checkpoint 3
    new_text = flip_text(infile="flip_check_dummy_text")
    assert new_text == "banana"


def genome_reads_main():

    filename = 'file_manager_test.txt'
    create_test_txt_file(filename)
    #prefix = os.getcwd().split('uniquekmer')[0]
    #filename = prefix + 'uniquekmer/test/test_file_manager.txt'
    byte_read_array_numpy(filename=filename)

    genome_reads_test(filename='data/genome.fa')


def byte_read_array_numpy(filename:str):
    array_array = read_byte_array(filename)
    numpy_array = read_byte_numpy(filename)

    assert type(numpy_array) is np.ndarray
    assert type(array_array) is array.ArrayType
    assert list(numpy_array) == list(array_array)


def msg_pack_ambs_unam(filename:str):

    a_u_dict = msgunpack_dict(filename)
    with open('a_u_dict.txt', 'w') as file:
        for unam in a_u_dict:
            file.write('unam: ')
            file.write(str(unam))
            file.write('\n')
            file.write('ambs ')
            file.write(str(a_u_dict[unam]))
            file.write('\n')

    return


def create_test_txt_file(filename):
    out_list = random.choices(range(1000000), k=50)

    with open(filename,'wb') as file:
        file.write(struct.pack('I', 50))
        for item in out_list:
            file.write(struct.pack('I', item))


def create_msgpack_test_file(filename):
    out_list = random.choices(range(1000000), k=50)
    print(out_list)
    with open(filename, 'wb') as file:
        msgpack.pack(out_list,file)


def genome_reads_test(filename):
    # testing reads() from file manager
    filename = append_file_name(filename)
    past = time.time()
    read_bioseq = read_unambiguous(filename=filename)
    current = time.time()
    print('genome read with Bio.Seq. Time elapsed: ',current - past )
    past = current


    read_reads = reads(filename=filename)
    current = time.time()
    print('genome read with Reads.py. Time elapsed: ', current - past)

    assert type(read_bioseq) is Bio.Seq.Seq

    assert type(read_reads) is str

    assert read_reads == str(read_bioseq)


def append_file_test():

    #appended = append_file_name('test')

    #print('appended: ', appended)
    #print('os.getcwd(): ', os.getcwd())

    #assert appended == os.getcwd()

    appended = append_file_name('test/append_test.txt')

    with open(appended, 'w') as file:
            for i in range(100):
                    file.write(str(i))
                    file.write('\n')


if __name__ == '__main__':
    #filename = append_file_name('src/chr22msgpack_ambs_unam')
    #msg_pack_ambs_unam(filename)

    #test_chr_splits()
    flip_check()



