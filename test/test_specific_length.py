import sys
import os
from tqdm import tqdm
from collections import deque

sys.path.append(os.getcwd().split('uniquekmer')[0]+'uniquekmer/src')
import specific_length as sl
from file_manager import append_file_name
from arg_reader import specific_unique_args


def naive_kmer_uniqueness_test(filename:str):

    seqs = {}
    #seqs = []
    with open(filename, 'r') as file:
        # ignore lines that do not begin with 'sequence: '
        print('reading file: ')
        address = 0
        ambs = 0
        for line in tqdm(file, desc='reading file'):

            # check for N's as well

            if line[0:10] == 'sequence: ':
                seq = line.split('sequence: ')[1]
                if 'N' in seq:
                    raise Exception('ambiguous char found in sequence at address: ' + str(address) + ', sequence: ' + seq)
                    ambs += 1
                #seqs.append(seq)
                if not seq in seqs:
                    seqs[seq]=address
                else:
                    raise Exception('NOT UNIQUE SEQUENCES FOUND: ' + seq + ' at addresses: ' + str(address) + ' and ' + str(seqs[seq]))

            elif line[0:9] == 'address: ':
                address = int(line.split('address: ')[1])

        if ambs:
            print('Ambiguous chars found in ', ambs, ' k-mers')

        print('sorting sequences: ')

    seqs_d = deque(sorted(seqs.values()))
    # seqs_d = deque(sorted(seqs))
    # with open('delete_me_next.txt','w') as file:
    #    for item in tqdm(seqs, desc='writing to file: '):
    #        file.write(item + '\n')
    total = len(seqs)
    pbar = tqdm(total=total, desc='ensuring uniqueness: ')
    past = ''
    not_unique = 0
    while seqs_d:
        pbar.update(1)
        seq = seqs_d.popleft()
        print(seq)
        if seq == past:
            raise Exception('NOT UNIQUE SEQUENCES FOUND: ' + past + ' and ' + seq + '\naddress1: ' + str(seqs[past]) + ' and address2: ' + str(seqs[seq]))
            #raise Exception('NOT UNIQUE SEQUENCES FOUND: ' + past + ' and ' + seq)
        past = seq

    print('Number of not unique pairs: ', not_unique, '.\nPercent not unique: ', str(not_unique/total))
    return


def naive_address_uniqueness_test(filename:str):

    print('reading dict: ')
    d = sl.read_dict(filename=filename, filetype='json')
    print('dict read!\ndequeing keys of dict: ')
    addresses = deque(sorted(d))
    print('dequed!')

    past = 0
    pbar = tqdm(total=len(addresses), desc='checking uniqueness of addresses: ')
    while addresses:
        pbar.update(1)
        adr = addresses.popleft()
        if adr == past:
            raise Exception('NOT UNIQUE ADDRESSES FOUND: ' + str(past) + ' and ' + str(adr))

        past = adr
    return


def with_args():
    args = specific_unique_args()

    print('reading dict:')
    d = sl.read_dict(filename=args.dfile, filetype=args.filetype)
    print('dict read!')

    output = '' if not args.outfile else append_file_name('output/' + args.outfile)
    valids, sequence, special_end = sl.specific_k(k=args.length, d=d, seq_file=args.genome, outfile=output)
    sl.file_write(valids=valids,outfile=output,sequence=sequence,k=args.length, special_end=special_end,d=d)
    naive_address_uniqueness_test(args.dfile)
    naive_kmer_uniqueness_test(output)

if __name__ == '__main__':
    naive_kmer_uniqueness_test(filename=append_file_name('output/0403_c22_30mer'))
