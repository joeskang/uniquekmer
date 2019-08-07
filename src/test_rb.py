"""
 read the trues, uniques, and tops files

 insert into rb_tree
"""

from rb_tree import RedBlackTree


import file_manager as fm
from tqdm import tqdm, trange


alpha = {
    "A":1,
    "C":2,
    "G":3,
    "T":4,
    "N":0
}

def main(filename):
    kmers = get_kmers(filename)
    rb = make_tree(kmers)
    return


def convert_geno_to_num(filename):
    geno = fm.reads(filename)
    return [alpha[_] for _ in tqdm(geno, desc="converting genome to numbers: ")]

def make_tree(kmers):
    rb = RedBlackTree()
    for k in tqdm(kmers, desc="adding to tree: "):
        if not rb.contains(k):
            rb.add(k)

    return rb

def get_kmers(filename):

    geno = convert_geno_to_num(filename)
    print("reading true addresses: ", end='')
    trues = fm.read_byte_to_queue(fm.append_file_name('0415_c22_1132_true_addresses'))
    print("done!\nreading unique start addresses: ", end='')
    uniques = fm.read_byte_to_queue(fm.append_file_name('0415_c22_1132_unique_starts'))
    print("done!")
    kmer_list = []

    for i in trange(len(trues), desc="finding kmers: "):
        t = trues[i]

        # kmer = filter(str.isdigit, repr(geno[t:t+uniques[i]]))
        kmer = int(''.join(map(str, geno[t:t+uniques[i]])))
        kmer_list.append(kmer)

    return kmer_list


if __name__ == "__main__":
    main(fm.append_file_name("data/22.fa"))
