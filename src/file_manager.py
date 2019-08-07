"""
    handles byte file io
    will read or write in 4 bytes

"""
import json
import array
from tqdm import tqdm, trange
import time
import struct
from collections import deque, OrderedDict
import pickle
from Bio.Seq import Seq
from Bio import SeqIO
import os
import numpy as np
import msgpack
import random
import re
from sys import path, platform
from math import ceil

path.append(os.getcwd().split('uniquekmer')[0] + '/uniquekmer/src')

from reads import get_fstream


# overriding reads() from reads.py
def reads(filename: str):
    fp, fstream = get_fstream(filename)
    sequence = ""
    for read in tqdm(fstream, desc='Reading genome file'):
        sequence += read.seq
    fp.close()
    return sequence


def pickle_dict(indict, filename: str):
    with open(filename, 'wb') as file:
        pickle.dump(indict, file, protocol=pickle.HIGHEST_PROTOCOL)
    return


def unpickle_dict(filename: str):
    with open(filename, 'rb') as file:
        return pickle.load(file)


def msgpack_dict(indict, filename: str):
    with open(filename, 'wb') as file:
        msgpack.pack(o=indict, stream=file)
    return


def msgunpack_dict(filename: str):
    with open(filename, 'rb') as file:
        return msgpack.unpack(stream=file)


def change_dir():
    # changes directory to ../data
    prefix = os.getcwd().split('uniquekmer')[0]
    return prefix + '/uniquekmer/data'


def append_file_name(filename: str) -> str:
    if not os.path.isfile(filename):
        filename = filename.split('uniquekmer')[0]

    return os.getcwd().split('uniquekmer')[0] + 'uniquekmer/' + filename


def read_byte_array(filename: str) -> array.array:
    """
        reads byte file in groups of 4 bytes
    """

    # the code below throws OSError Errno 22 when file is large
    with open(filename, 'rb') as inf:
        buf = inf.read(4)
        length = struct.unpack('I', buf)[0]
        read_array = array.array('I')
        read_array.fromfile(inf, length)
        return read_array


def read_byte_numpy(filename: str) -> np.ndarray:
    with open(filename, 'rb', buffering=0) as file:
        buf = file.read(4)
        count = struct.unpack('I', buf)[0]
        return np.fromfile(file, dtype=np.uint32, count=count).ravel()


def write_array_to_byte(byte_arr, filename):
    """
        reads values from bool_q in 4 bytes and writes a byte file
    :param byte_arr:
    :param out_file:
    :return: none:
    """

    out_file = open(filename, 'wb')
    for bytes in tqdm(byte_arr, desc='Writing array to bytes: '):
        out_file.write(struct.pack('I', bytes))

    return


def read_byte_to_queue(filename: str) -> list:
    return_queue = []
    try:
        with open(filename, 'rb') as f:
            while True:
                buf = f.read(4)
                if not buf: break
                return_queue.append(struct.unpack('I', buf)[0])
    except struct.error:
        raise
    finally:
        return return_queue


def read_big_endian(filename):
    return_queue = deque()
    with open(filename, 'rb') as f:
        while True:
            buf = f.read(4)
            if not buf: break
            return_queue.append(buf)

    return return_queue


def write_queue_to_byte(filename, in_deque):
    try:
        if type(in_deque) == deque:
            with open(filename, 'wb') as f:
                while in_deque:
                    f.write(in_deque.popleft())
        else:
            with open(filename, 'wb') as f:
                for i in range(len(in_deque)):
                    f.write(in_deque[i])

    except IndexError:
        return


def write_queue_to_int(filename, in_deque):
    try:
        with open(filename, 'w') as f:
            while in_deque:
                f.write(str(in_deque.popleft()) + '\n')

    except IndexError:
        return


def read_unambiguous(filename):
    """
        read from file and truncate to only A, C, G, or T
    """
    full = reads(filename=filename)
    seq_split = full.split('N')
    seq = ''
    for part in seq_split:
        seq += part

    return seq


def bioseq_read(filename):
    sequence = Seq("")
    for seq_record in tqdm(SeqIO.parse(filename, "fasta"), desc='Reading FASTA file'):
        sequence += seq_record.seq

    return sequence


def json_it(data, filename: str):
    # writing to output file
    try:
        json.dump(data, open(filename, 'w'))

    except IndexError:
        raise


def unjson_it(filename: str):
    # reading from json file
    try:
        return json.load(open(filename, 'r'))

    except IndexError:
        raise


def read_anything(filename: str):
    # just read anything and everything from file
    try:
        r_str = ''
        with open(filename, 'r') as file:
            for line in file:
                r_str = line.strip()

        return r_str

    except Exception:
        raise


def flip_fasta(filename=str, outfile=str):
    # flip a fasta file
    try:
        text = ''
        with open(filename, 'r') as file:
            for line in file:
                if line[0] == '>':
                    pass

                else:
                    text += line.strip()

            text = flip_text(intext=text)

        with open(outfile,'r') as file:
            file.write(text)

    except Exception:
        raise


def flip_text(intext=None, infile=None, outfile=None):
    if not intext and not infile and not outfile:
        raise Exception("No values passed")

    elif intext and infile:
        raise Exception("Both input file and text received")

    try:
        text = read_anything(infile) if infile else intext
        seq = ''
        for _ in range(len(text)):
            seq += text[-_ - 1]

        if outfile:
            with open(outfile, 'w') as file:
                file.write(seq)
        else:
            return seq

    except Exception:
        raise


def read_kmer_custom_file(filename: str):
    # read the custom made k-mer table file
    seqs = []
    with open(filename, 'r') as file:
        for line in tqdm(file, desc='reading file'):

            # check for N's as well

            if line[0:10] == 'sequence: ':
                seq = line.split('sequence: ')[1]
                if 'N' in seq:
                    raise Exception('ambiguous char found in sequence: ' + seq)
                # seqs.append(seq)
                if not seq in seqs:
                    # seqs[seq]=address
                    seqs.append(seq.split()[0])
                else:
                    raise Exception('NOT UNIQUE SEQUENCES FOUND: ' + seq)

            elif line[0:9] == 'address: ':
                pass

    return seqs


def fake_genome(filename: str, lines: int, alpha=None):
    # first check if the file is already made
    # if not os.path.isfile(filename):

    if not alpha:
        alpha = {
            0: 'N',
            1: 'A',
            2: 'C',
            3: 'G',
            4: 'T'
        }
    flag = True
    with open(filename, 'w') as file:
        # initialize first chr
        chr_no = 1
        file.write('>chr')
        file.write(str(chr_no))
        file.write('\n')

        # print the lines of A,C,G,T,N
        for i in range(lines):
            # choose whether to start a new chr
            print_chr = random.choice(range(50))
            if print_chr == 7 and i > 0:
                chr_no += 1
                file.write('>chr')
                file.write(str(chr_no))
                file.write('\n')
                flag = True

            # flag determines whether we are beginning of chr,
            # if so, output ambs

            string = ''
            # make array of 60 random ints from 0-4
            myarr = random.choices(range(5), k=60)

            if flag:
                flag = False
                how_many = random.choice(range(3, 10))
                for i in range(how_many):
                    myarr[i] = 0
            for item in myarr:
                # use the alphabet dictionary to translate the ints into chars
                string += alpha[item]

            file.write(string)
            file.write('\n')

    return


def truncate(infile: str, outfile: str):
    # truncates the input sequence file and outputs to outfile
    # 60 chars/line

    sequence = reads(infile)
    unam = re.split(r'N+', sequence)
    unam = deque(filter(None, unam))

    with open(outfile, 'w') as file:
        while unam:
            count = 0
            item = unam.popleft()
            for char in item:
                count += 1

                file.write(char)

                if count == 60:
                    file.write('\n')
                    count = 0

    return


def simple_string_read(filename: str):
    seq = ''
    with open(filename, 'r') as file:
        for line in file:
            seq += line

    return seq


def chr_splits(filename: str):
    chrs = OrderedDict()
    sequence = ''

    with open(filename, 'r') as file:
        chr_info = ''
        seq = ''
        for line in file:

            if line[0] == '>':
                if chr_info and seq:
                    chrs[chr_info.split()[0]] = len(seq)

                elif not seq and chr_info:
                    raise Exception("Found Empty Chromosome!")

                elif not chr_info and seq:
                    raise Exception("Found Chromosome Without Chromosome Info!")

                # reset
                sequence += seq
                seq = ''
                chr_info = line

            else:
                seq += line.split()[0]
        chrs[chr_info.split()[0]] = len(seq)
        sequence += seq

    return chrs, sequence


def write_2way_genome(infile:str, outfile:str):
    geno = reads(infile)
    flipped = flip_text(geno)

    with open(outfile,'w') as file:

        file.write('>' + outfile + '\n')

        lines = ceil(len(geno)/60)

        # write the forward genome first
        for i in range(lines):
            file.write(geno[i*60:i*60+60] + '\n')

        # write the flipped genome next
        file.write('>flipped genome\n')
        for i in range(lines):
            file.write(flipped[i*60:i*60+60]+'\n')


def json_to_text(infile:str, outfile:str):
    data = unjson_it(infile)

    with open(outfile, 'w') as file:
        if type(data) == list and type(data[0]) == list:
            for myl in data:
                file.write()
                # TODO: continue if necessary


if __name__ == '__main__':
    # make fake genome
    # filename = append_file_name('test/fake_genome')
    # fake_genome(filename=filename, lines=2000)

    # print("writing 2way genome")
    # write_2way_genome(infile=append_file_name('data/genome.fa'), outfile=append_file_name('data/2way_genome.fa'))
    # write_2way_genome(infile=append_file_name('data/dummy'), outfile=append_file_name('data/2way_dummy'))
    print("do not run this file directly")
