"""
    run through chr22 and get the line numbers where unam and ambs group ends.


"""
from collections import deque, OrderedDict
from tqdm import tqdm
import file_manager as fm
import time
import re
from rb_tree import RedBlackTree

def main():
    ambs, unam = ambs_unam_split(True, '../data/22.fa')


    ##########
    print(len(ambs))
    print(len(unam))

def ambs_unam_split(filename:str, chr_option=False):
    """
        chr_option: take into consideration lines for chr labels, not taken into consideration by default
        
    :param chr_option:
    :param filename: should be the *.fa file and not the bt2 index file
    :return:


        Deprecated
    """

    # TODO: outputting weird data. fix
    try:
        with open(filename, 'r') as inf:
            # flag that tells if we're in ambs or unam group
            # True = ambs
            # False = unam
            ambs_flag = True

            #ambs_end = deque()
            #unam_end = deque()
            ambs = 0
            line_num = 0
            unam = 0
            past = 0
            unam_offset = 0
            unam_line = 0
            a_u_dict = OrderedDict() 
            # reading from input file
            for line in tqdm(inf, desc="Finding AMBS/UNAM groups"):
                # implement counter that will save the line number
                column = 0

                if line[0] == '>':
                    # preserve the line numbers if the chr info lines are to be counted
                    if chr_option:
                        line_num += 1
                    continue

                # switch to unam from ambs
                elif ambs_flag and ('A' in line or 'G' in line or 'C' in line or 'T' in line):
                    ACGT = [None]*4
                    ACGT[0] = line.find('A')
                    ACGT[1] = line.find('C')
                    ACGT[2] = line.find('G')
                    ACGT[3] = line.find('T')

                    ACGT = sorted(ACGT)

                    ambs_flag = False
                    # need to subtract 1 from the result of find()
                    #   as find returns the position of the first char in the new group
                    # we need the position of end of the previous group

                    # ambs < unam
                    ambs = line_num * 60 + ACGT[0] - 1 - unam_line + ambs

                # switch to ambs from unam
                elif not ambs_flag and 'N' in line:
                    ambs_flag = True
                    """
                    Is the below implementation necessary?
                    
                    # need to subtract 1 from the result of find()
                    #   as find returns the position of the first char in the new group
                    # we need the position of end of the previous group
                    current = line_num*60 + line.find('N') - 1
                    unam_offset += current - past
                    past = current
                    #unam_end.append(unam_offset)
                    a_u_dict[unam_offset] = ambs
                    """

                    # ambs < unam
                    unam_line = line_num*60 + line.find('N') - 1
                    unam = past - ambs
                    # ambs should have been initialized by now, thus no need for 'past'
                    a_u_dict[unam] = ambs

                line_num += 1

            # Keys: unam
            # Values: ambs
            return a_u_dict

    except FileNotFoundError as e:
        raise SystemExit(e)

    except IndexError as e:
        # there should be no index errors with given algorithm
        raise(e)

    except KeyError:
        raise('ambs key does not exist for unam value')


def split_sequence(sequence=None, filename=None):
    """
        a more concise implementation of above with regex

    :param sequence:
    :param filename:
    :return:
    """

    if not sequence and not filename:
        raise(Exception('not enough values passed to split_sequence()'))

    if filename and not sequence:
        sequence = ''
        with open(filename, 'r') as file:
            for line in file:
                if line[0] != '>':
                    sequence += line.strip()

    # split the sequence to unam groups
    unam_list = re.split(r'N+', sequence)
    unam_list = deque(filter(None, unam_list))
    ambs_list = re.split(r'[ACGT]+', sequence)
    ambs_list = deque(filter(None, ambs_list))

    print('len(unam_list): ', str(len(unam_list)))
    print('len(ambs_list): ', str(len(ambs_list)))

    ambs_list = deque([len(ambs_list[_]) for _ in range(len(ambs_list))])

    unam_list = deque([len(unam_list[_]) for _ in range(len(unam_list))])

    return ambs_list, unam_list


def rb_tree_ambs(ambs:deque, unam:deque):
    rb = RedBlackTree()
    u_off = -1
    a_off = 0
    # assuming that genome starts with ambs
    while ambs and unam:
        a_off += ambs.popleft()
        u_off += unam.popleft()
        pair = [u_off, a_off]
        rb.add(pair)

    return rb


def test_rb_ambs(rb:RedBlackTree, sequence:str):
    myl = list(rb)
    true_adr = 0
    past = 0
    for pair in myl:
        if not true_adr:
            true_adr = pair[1] - 1
        else:
            true_adr = pair[1] + past
        assert sequence[true_adr] == 'N'
        true_adr = pair[0] + pair[1]
        past = pair[0]
        assert sequence[true_adr] != 'N'

    return


def test_validity(d:dict, sequence:str):
    print("testing validity of split")
    a_off = -1
    u_off = 0
    for unam in d:
        a_off += (d[unam] + unam)
        curr = sequence[a_off]
        surround = sequence[a_off-5:a_off+5]
        dum = sequence[a_off:a_off+5]
        assert curr == 'N'
        if u_off:
            u_off = a_off - d[unam]
            curr = sequence[u_off]
            surround = sequence[u_off-5:u_off+5]
            dum = sequence[u_off:u_off+5]
            assert curr != 'N'
        else:
            u_off = a_off - d[unam]

    return


def with_args():

    seq = fm.reads(input('Enter file name of sequence file: '))
    d = split_sequence(sequence=seq)
    print('writing to file: ')
    with open(fm.append_file_name('ambs_unam'),'w') as file:
        unam = list(d.keys())
        ambs = list(d.values())
        file.write('unam\n')
        for u in unam:
            file.write(str(u) + '\n')

        file.write('\n\nambs\n')
        for a in ambs:
            file.write(str(a) + '\n')


if __name__ == '__main__':
    ambs, unam = split_sequence(filename='../data/22.fa')
    rb = rb_tree_ambs(ambs, unam)
    print(list(rb))
    seq = fm.reads('../data/22.fa')
    test_rb_ambs(rb,sequence=seq)



