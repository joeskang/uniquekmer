"""
    get lcp_sa_dict and output k-mers of specific length k
        likely required to get default dict with unique starts instead

"""
import sys
import os
from tqdm import tqdm, trange
from collections import deque
import pdb
from itertools import repeat
from multiprocessing import Pool, Process

sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')
from file_manager import append_file_name, reads, unpickle_dict, unjson_it, msgunpack_dict, read_byte_numpy, read_byte_to_queue, json_it
from exception_manager import InsufficientArguments
from arg_reader import specific_unique_args
#from plot_manager import draw_histo_list, draw_bar_list


def specific_k(k: int, d: dict, seq_file='', sequence='', high=0, outfile=append_file_name('k_mer')):
    """
        d: { key( true address of genome ): [ unique start, top ] }
    :param k:
    :param d:
    :return:
    """
    # if the sequence has been passed, then use sequence.
    # else, if only seq_file was passed, then read sequence
    if not sequence and seq_file:
        sequence = reads(seq_file)
        print('length of sequence: ' + str(len(sequence)))

    elif sequence:
        pass

    else:
        raise InsufficientArguments

    valids = {}
    # keys: seq
    # values: sa
    try:
        special_end = 0
        for sa in tqdm(d, desc='finding ' + str(k) + '-mers: '):
            # first find all the k-mers
            if d[sa][1] != high:
                special_end += 1
            if d[sa][0] + 1 < k < d[sa][1]:
                seq = sequence[int(sa):int(sa) + k]
                if seq in valids:
                    del valids[seq]
                    continue
                else:
                    valids[seq] = sa

    except IndexError:
        pass

    return valids


def file_write(valids: dict, outfile: str, k, d: dict):
    with open(outfile, 'w') as file:
        # next, output the total number of valid k-mers
        first_line = 'Number of ' + str(k) + '-mers: ' + str(len(valids))
        file.write('____________________________\n')
        file.write(first_line)
        file.write('\n____________________________\n\n\n')

        myd = deque()

        pbar = tqdm(total=len(valids), desc='Ensuring uniqueness')
        for seq in valids:
            pbar.update(1)
            address = valids[seq]
            if len(seq) != k:
                continue
            file.write('address: ' + str(address))
            file.write('\nsequence: ')
            if 'N' in seq:
                raise (Exception('N found at: ' + str(address) + ': ' + seq + '\nWith top: ' + str(d[address][1])))
            file.write(seq + '\n\n')


def read_dict(filename: str, filetype='pickle'):
    if filetype.lower() == 'pickle':
        return unpickle_dict(filename)

    elif filetype.lower() == 'json':
        return unjson_it(filename)

    elif filetype.lower() == 'msgpack':
        return msgunpack_dict(filename=filename)

    else:
        raise Exception('unknown file type')


def conseq_count(trues, tops):
    count = []
    clusters=[]
    u_tops=[]
    # key: start address
    # value: size of cluster (number of adjacent k-mers), top

    n = 0
    cluster = 0
    start = trues[0]
    hold_tops = []
    u_top = int()

    max_cluster_size = 0
    for i in range(len(trues)):


        curr = trues[i]
        if curr == n + 1:
            cluster += 1
            hold_tops.append(tops[i])

        elif n and hold_tops:
            u_top = hold_tops.pop()

            while hold_tops:
                # work backwards and get the tops
                # u_top = top that ensures that the cluster is unique internally
                # if curr top > curr past, then set the
                top_p = hold_tops.pop()
                u_top = u_top if u_top <= top_p else top_p

            max_cluster_size = max_cluster_size if max_cluster_size >= cluster else cluster
            #count[start] = (cluster, u_top)
            count.append(start)
            clusters.append(cluster)
            tops.append(u_top)
            cluster = 0
            start = curr

        n = curr
    print("max size of cluster: ", max_cluster_size)
    print("number of cluster: ", len(count))
    return count, clusters, u_tops



def bytes_unique(trues, tops, lcps, full_seq):
    counts = 0

    # begin with first address and see how many k-mers have consecutive addresses


def read_data(true_file: str, tops_file: str, lcps_file:str, chrs_file:str):
    print("Reading data: ")
    print("Reading True Address file: ", end="")
    trues = read_byte_to_queue(true_file)
    print("done.\nReading Tops file: ", end="")
    tops = read_byte_to_queue(tops_file)
    print("done.\nReading Unique Start file: ", end='')
    lcps = read_byte_to_queue(lcps_file)
    print("done.\nReading chrs file: ", end='')
    chrs_dict = unjson_it(chrs_file)
    print("done.")
    return trues, tops, lcps, chrs_dict

def read_trues_chr(true_file: str, chrs_file:str):

    print("Reading data: ")
    print("Reading True Address file: ", end="")
    trues = read_byte_to_queue(true_file)

    print("done.\nReading chrs file: ", end='')
    chrs_dict = unjson_it(chrs_file)
    print("done.")
    return trues, chrs_dict


def bytes_main():
    args = specific_unique_args()
    #trues, tops, lcps, chrs_d = read_data(chrs_file= append_file_name(args.chrs), true_file=append_file_name(args.trues), tops_file=append_file_name(args.tops), lcps_file=append_file_name(args.unqs))
    trues, chrs_d = read_trues_chr(chrs_file=append_file_name(args.chrs), true_file=append_file_name(args.trues))

    #counts, chr_counts = _run_with_trues(trues=trues, lcps=lcps, tops=tops, chr_d=chrs_d, bot=args.low, top=args.high)
    chr_counts = _count_trues(trues=trues, chr_d=chrs_d)
    #print("counts of uniques: ", end='')
    #print(counts)

    print('chrs: ')
    print(chr_counts)




def _run_through_internal_no_trues(tops:list, lcps:list, bot=20, top=100):
    counts = [0 for _ in range(top - bot + 1)]
    assert len(lcps) == len(tops)
    for j in trange(len(lcps)):
        for i in range(top-bot+1):
            if lcps[j] < i+bot <= tops[j]:
                counts[i] += 1

    return counts

def _run_with_trues(trues:list, lcps:list, tops:list, chr_d:dict, bot=20, top=100):
    chrs = list(chr_d)
    counts = [0 for _ in range(top-bot+1)]
    assert len(lcps) == len(tops) == len(trues)
    chr_counts = {}


    for i in range(top - bot + 1):
        text = "finding " + str(i+bot) + "-mers:"

        chr_ = 0
        curr_chr = chrs[chr_]
        chr_len = chr_d[curr_chr]
        chr_ += 1

        count = 0

        #for t in trange(len(trues),desc=text):
        print(text)
        for t in range(len(trues)):
            if lcps[t] < i+bot <= tops[t]:
                count += 1

                if trues[t] <= chr_len:
                   chr_counts[curr_chr] = 1 if not curr_chr in chr_counts else chr_counts[curr_chr] + 1

                else:
                    curr_chr = chrs[chr_]
                    chr_len += chr_d[curr_chr]
                    chr_ += 1

        counts[i] = count

    return counts, chr_counts


def _count_trues(trues:list, chr_d:dict):
    chrs = list(chr_d)
    chr_counts = {}

    chr_ = 0
    curr_chr = chrs[chr_]
    chr_len = chr_d[curr_chr]
    chr_ += 1

    for t in trues:
        if t <= chr_len:
            chr_counts[curr_chr] = 1 if not curr_chr in chr_counts else chr_counts[curr_chr] + 1

        else:
            curr_chr = chrs[chr_]
            chr_len += chr_d[curr_chr]
            chr_ += 1

    return chr_counts


def SNP_uniquemer(trues:list, tops:list,filename:str, lcps:list, full_seq:str, snp_addr:list, SNPs:list, bot=20, top=100):
    with open(filename,'w') as file:
        t_adr = len(trues)
        top = tops[-1]
        while snp_addr and SNPs:

            # get an address from SNP list
            adr = snp_addr.pop()
            snp = SNPs.pop()

            # ensure that trues has the address:
            # pop until less than or equal to SNP list address
            while t_adr > adr-top and trues and tops:
                t_adr = trues.pop()
                top = top.pop()
            # TODO: complete
            #    if t_adr




    # if adr = i, then check for each address from i-top+1 to i

def test_run(d: dict, full_seq: str, s_arr, bot=20, top=23):
    counts = {}
    runs = len(s_arr)
    for _ in range(top - bot):
        i = bot + _
        count = 0
        past = ''
        text = "finding " + str(i) + "-mers:"
        for sa in tqdm(total=runs, iterable=s_arr, desc=text):
            if d[str(sa)][0] < i < d[str(sa)][1]:
                seq = full_seq[sa:sa + i]
                assert 'N' not in seq and seq != past
                past = seq
                count += 1
        print(count)
        counts[i] = count

    return counts


def main():
    args = specific_unique_args()
    print('reading dict:')
    d = read_dict(filename=args.dfile, filetype=args.filetype)
    print('dict read!')

    output = '' if not args.outfile else args.outfile

    # get the number of times to loop
    cycle = 0
    if args.low and args.high:
        if args.low == args.high:
            cycle = 1
        else:
            assert args.high > args.low
            cycle = args.high - args.low + 1
    output = output + '_' + str(args.low) + '-mer'
    for i in range(cycle):
        if i > 0:
            spl = str(args.low + i - 1) + "-mer"
            output = output.split(spl)[0] + str(args.low + i) + "-mer"
        valids = specific_k(k=int(args.low + i), d=d, seq_file=args.genome, outfile=output, high=args.high)
        file_write(valids=valids, outfile=output, k=args.low + 1, d=d)


if __name__ == '__main__':
    """
    print('reading dict: ')
    d = read_dict(filename='../0409_c22_1115_json_default_dict', filetype='json')
    print('dict read!')

    full_seq = reads('../data/22.fa')
    s_arr = read_byte_numpy(filename=append_file_name('data/22.sa'))
    # print(test_run(d=d,full_seq=full_seq, s_arr=s_arr))
    print(run_through_pool(d=d, full_seq=full_seq, out_file=''))
    """
    bytes_main()
