# test the module ambiguous_split.py
# contains: ambs_unam_split()
import sys
import os

sys.path.append(os.getcwd().split('test')[0]+'src')

from ambiguous_split import ambs_unam_split, split_sequence
import os.path
import file_manager as fm
import random
from csv import DictWriter, writer

alpha = {
    0:'N',
    1:'A',
    2:'C',
    3:'G',
    4:'T'
}

def test_ambs_unam_split():
    filename = fm.append_file_name('test/fake_genome')
    # make a fake genome sequence with file name test_a_u_split.txt
    if not os.path.isfile(filename):
        fm.fake_genome(filename=filename, lines=1000,alpha=alpha)
    try:

        ambs, unam = ambs_unam_split(chr_option = False, filename=filename)

        assert ambs is not None and unam is not None

        assert len(ambs) is not 0 and len(unam) is not 0

        # test to see that the entries in each deque are increasing
        past_amb = ambs.popleft()
        past_una = unam.popleft()
        while ambs:
            cur_amb = ambs.popleft()
            assert cur_amb > past_amb
            past_amb = cur_amb

            cur_una = unam.popleft()
            assert cur_una > past_una
            past_una = cur_una



    except IndexError as e:
        print(e)

def create_fake(filename:str):
    sequence = ''

    # create fake genome with known ambs unam groups
    # thus, unable to simply to file_manager's fake_genome()
    # create about 40 ambs/unam groups:

    # lists of the lengths of ambs/unam groups
    ambs = random.choices(population=range(1, 50), k=40)
    unam = random.choices(population=range(50, 100), k=40)

    file1 = fm.append_file_name('test/'+filename+'ambs_unam')

    # write the ambs/unam lists to a file
    with open(file1, 'w') as file:
        file.write('ambs: ')
        for a in ambs:
            file.write(str(a) + ' ')

        file.write('\nunam: ')
        for u in unam:
            file.write(str(u) + ' ')

    file2 = fm.append_file_name('test/'+filename)

    # start with ambs
    for i in range(len(ambs)):
        a = ambs[i]
        for _ in range(a):
            sequence += 'N'
        u = unam[i]
        for _ in range(u):
            sequence += random.choice(['A', 'C', 'G', 'T'])


    # split string into 60 chars, as the fasta file is
    lines = [sequence[_:_+60] for _ in range(0, len(sequence), 60)]

    with open(file2,'w') as file:
        for line in lines:
            file.write(line + '\n')

def compare(filename:str):
    d = ambs_unam_split(filename=filename)
    print(d)
    ambs = list(d.values())
    unam = list(d.keys())
    print('ambs: ',end='')
    print(ambs)
    print('unam: ', end='')
    print(unam)

    with open('ambs_unam','w') as file:
        file.write('ambs\n')
        for i in range(len(ambs)):
            file.write(str(ambs[i])+'\n')

        file.write('\n\n')
        file.write('unam\n')
        for i in range(len(unam)):
            file.write(str(unam[i])+'\n')

def compare2(filename:str):
    d = split_sequence(filename=fm.append_file_name('test/'+filename))
    ambs = list(d.values())
    unam = list(d.keys())

    naive_a, naive_u = naive_ambiguous_split(filename=filename)

    assert len(ambs) == len(naive_a)
    assert len(unam) == len(naive_u)

    for i in range(len(ambs)):
        assert ambs[i] == naive_a[i]

    for i in range(len(unam)):
        assert unam[i] == naive_u[i]

    with open('ambs_unam_2','w') as file:
        file.write('ambs\n')
        for i in range(len(ambs)):
            file.write(str(ambs[i])+'\n')

        file.write('\n\n')
        file.write('unam\n')
        for i in range(len(unam)):
            file.write(str(unam[i])+'\n')

def naive_ambiguous_split(filename:str):
    ambs = []
    unam = []
    is_amb= False
    first = True

    with open(fm.append_file_name('test/'+filename),'r') as file:
        a_count = -1
        u_count = 0

        for line in file:
            for char in line.strip():
                if char == 'N':

                    if u_count and not is_amb:
                        first = False
                        is_amb = True
                        unam.append(u_count)

                    a_count += 1


                elif char == 'A' or char == 'C' or char == 'G' or char == 'T':

                    if (a_count and is_amb) or first:
                        is_amb = False
                        first = False
                        ambs.append(a_count)

                    u_count += 1

        if a_count not in ambs:
            ambs.append(a_count)

        if u_count not in unam:
            unam.append(u_count)

    return ambs, unam


def read_json_file(filename:str):
    a_u_dict = fm.unjson_it(filename)
    with open('ambs', 'w') as ambs:
        with open('unam', 'w') as unam:
            for key in a_u_dict:
                unam.write(str(key)+'\n')
                ambs.write(str(a_u_dict[key])+'\n')

def crawler(filename:str, geno_file:str):
    # goes through ambs/unam json file to ensure that the validity of ambs/unam groups
    a_u_dict = fm.unjson_it(filename)
    seq = fm.reads(geno_file)
    max = len(seq)
    u_offset = 0
    a_offset = 0 # important to note that ambs measures length and not position


    for una in a_u_dict:
        # key = unam
        # value = ambs



        u_offset += int(una)
        a_offset += (int( a_u_dict[una] ) + int(una) - 1)


        if u_offset != 0:
            assert seq[u_offset] in ['A', 'C', 'G', 'T']

            if u_offset < max:
                assert seq[u_offset + 1] == 'N'

        assert seq[a_offset] == 'N'

        if a_offset < max:
            assert seq[a_offset + 1] in ['A', 'C', 'G', 'T']

        u_offset += int(a_u_dict[una])

    return


if __name__ == '__main__':
    #lines = 1000
    #test_ambs_unam_split()
    #filename = 'test_a_u_split_2'
    #create_fake(filename=filename)
    #compare(filename)
    #compare2(filename)

    #read_json_file(fm.append_file_name('0402_c22_1449_json_ambs_unam'))

    crawler(filename='../0402_c22_1544_json_ambs_unam', geno_file='../data/22.fa')