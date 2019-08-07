import argparse


def sa_lcp_read_args():
    """
        Parses arguments from command line.
    """

    parser = argparse.ArgumentParser()
    parser.add_argument('-g', dest='genome', metavar='genome_file', help='genome file')
    parser.add_argument('-s', dest='SA', metavar='SA_file', help='suffix array file')
    parser.add_argument('-o', dest='outfile', action='store', help='output file')
    parser.add_argument('-m', dest='length', action='store',
                        help='maximum length of k-mer (default = 100)')
    parser.add_argument('-l', dest='low', action='store',
                        help='minimum length of k-mer (default = 20)')
    parser.add_argument('-t', dest='threads', action='store',
                        help='number of threads (default = 1)')
    parser.add_argument('-b', dest='bt2', action='store', help='bt2 index file')
    parser.add_argument('-a', dest='ambs',metavar='ambs_infile',help='File name for ambs input')
    parser.add_argument('-lcp', dest='lcpfile',metavar='lcp_file',help='File name for lcp array (if it exists)')
    parser.add_argument('-d', dest='distance',metavar='distance',help='Maximum distance between k-mers')
    parser.add_argument('-v', dest='verbose',metavar='verbose',help='Print comments (T/F)')
    parser.add_argument('-p', dest='pickle',metavar='pickle',help='Name of pickle file')
    parser.add_argument('-inv', dest='inverse',metavar='inverse',help='File name for inverse suffix array')
    args = parser.parse_args()


    if not args.length:
        args.length = 100
    else:
        args.length = int(args.length)

    if not args.low:
        args.low = 16
    else:
        args.low = int(args.low)

    if not args.distance:
        args.distance=300
    else:
        args.distance = int(args.distance)

    if not args.threads:
        args.threads = 1
    else:
        args.threads = int(args.threads)

    if not args.verbose:
        args.verbose = True

    elif args.verbose.lower() == 'f' or args.verbose.lower() == 'false' or args.verbose.lower() == 'n' or args.verbose.lower() == 'no':
        args.verbose = False

    return args





def specific_unique_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('-o', dest='outfile', metavar='output_file', help="name of output file")
    parser.add_argument('-low',dest='low',metavar='low_end_of_range', help='length of the shortest k-mer')
    parser.add_argument('-high',dest='high',metavar='high_end_of_range', help='length of the longest k-mer')
    parser.add_argument('-true', dest='trues', metavar='trues_file', help='filename of trues file')
    parser.add_argument('-top', dest='tops', metavar='tops_file', help='filename of tops file')
    parser.add_argument('-unq', dest='unqs', metavar='uniques_file', help='filename of unique starts')
    parser.add_argument('-chr', dest='chrs', metavar='chrs_file', help='filename of chrs info')

    args = parser.parse_args()

    args.low = int(args.low) if args.low else 20

    args.high = int(args.high) if args.high else 100

    #if not args.top and args.unq and args.true and args.chr:
    #   raise Exception('Not enough input files!')

    if not (args.unqs and args.tops and args.chrs and args.chrs):
        raise FileNotFoundError("Not enough input files!")

    return args

def ambs_unam():

    parser = argparse.ArgumentParser()
    parser.add_argument('-a', dest='ambsout',metavar='ambs_outfile',help='File name for ambs output')
    parser.add_argument('-u', dest='unamout',metavar='unam_outfile',help='File name for unam output')
    parser.add_argument('-b', dest='bt2file',metavar='bt2_infile',help='File name for bt2 input file')
    parser.add_argument('-t', dest='top',metavar='top',help='Maximum length of unique k-mers')
    parser.add_argument('-g', dest='genome', metavar='genome_file', help='File name for genome (FASTA)')

    args = parser.parse_args()


    # set default maximum length of unique kmer
    if not args.top:
        args.top = 100



