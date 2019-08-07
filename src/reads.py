#!/usr/bin/env python

#
# Copyright 2019, Chanhee Park <parkchanhee@gmail.com> and Daehwan Kim <infphilo@gmail.com>
#
# This file is part of HISAT 2.
#
# HISAT 2 is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# HISAT 2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with HISAT 2.  If not, see <http://www.gnu.org/licenses/>.
#

# small change made to print statements to make compatible with Python3


import os, sys, gzip, bz2
from argparse import ArgumentParser, FileType

"""
"""
COMPRESSION_NON   = 0
COMPRESSION_GZIP  = 1
COMPRESSION_BZIP2 = 2

SEQUENCE_FASTA    = 0
SEQUENCE_FASTQ    = 1

"""
"""
class Read:
    
    def __init__(self):
        self.init()
        return

    def __init__(self, id, seq):
        self.init()
        self.id = id
        self.seq = seq
        return

    def init(self):
        self.id = None
        self.seq = None
        return

    def __repr__(self):
        return '[{}, {}]'.format(self.id, self.seq)



"""
"""
def parser_FQ(fp):
    # skip empty line
    while True:
        line = fp.readline()

        if line == "":
            # end of file
            return

        if line[0] == '@':
            break;

    while True:
        id = line[1:].split()[0]
        seq = ""

        line = fp.readline()
        if line == "":
            return

        seq = line.strip()
        yield Read(id, seq)

        line = fp.readline() # '+'
        line = fp.readline() # quality
        line = fp.readline() # next ID
        if line == "":
            return

"""
"""
def parser_FA(fp):
    # skip empty line
    while True:
        line = fp.readline()

        if line == "":
            # end of file
            return

        if line[0] == '>':
            break;

    while True:
        id = line[1:].split()[0]
        seq = ""

        while True:
            line = fp.readline()
            if line == "":
                break

            if line[0] == '>':
                break

            seq += line.strip()

        yield Read(id, seq)

        if line == "":
            return

"""
"""
def parse_type(fname):
    compression_type = COMPRESSION_NON
    sequence_type = SEQUENCE_FASTA

    ff = fname.split('.')

    ext = ff[-1]
    if ext.lower() == "gz":
        compression_type = COMPRESSION_GZIP
        ext = ff[-2]
    elif ext.lower() == "bz2":
        compression_type = COMPRESSION_BZIP2
        ext = ff[-2]

    if ext.lower() == "fq" or \
        ext.lower() == "fastq":
        sequence_type = SEQUENCE_FASTQ

    return sequence_type, compression_type


"""
"""
def get_fstream(fname):
    sequence_type, compression_type = parse_type(fname)

    if compression_type == COMPRESSION_GZIP:
        fp = gzip.open(fname, 'r')
    elif compression_type == COMPRESSION_BZIP2:
        fp = bz2.BZ2File(fname, 'r')
    else:
        assert (compression_type == COMPRESSION_NON)
        fp = open(fname, 'r')

    if sequence_type == SEQUENCE_FASTA:
        fstream = parser_FA(fp)
    else:
        assert (sequence_type == SEQUENCE_FASTQ)
        fstream = parser_FQ(fp)

    return fp, fstream


"""
"""
def reads(read_fname, read_count = 0):
    fp, fstream = get_fstream(read_fname)
    for read in fstream:
        print(read)

    fp.close()
    return



if __name__ == '__main__':
    parser = ArgumentParser(
            description='')

    parser.add_argument('-f',
                        dest='read_fname',
                        action='store',
                        type=str,
                        help='read filename')

    args = parser.parse_args()

    if not args.read_fname:
        parser.print_help()
        exit(1)

    reads(args.read_fname)

