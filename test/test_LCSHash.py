import sys
import os
sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')

from experimental.LCSHash import LCSHash, Block


def test_LCSHash():
    lcs = LCSHash()

    common_substring = "ACAGTCGGAT"
    kmer1 = "T" + common_substring + "A"
    kmer2 = "A" + common_substring + "A"
    kmer3 = "T" + common_substring + "C"
    kmer4 = "C" + common_substring + "G"

    lcs.insert(kmer1)
    lcs.insert(kmer2)
    lcs.insert(kmer3)
    lcs.insert(kmer4)

    lcs.insert('TCGCGTAACGCGT')
    # lcs.insert('T' + common_substring + 'A')

    lcs.insert("ACGCGTAACGCGT", is_snp=True)

    lcs.insert("ACGTTTTTTCA", is_snp=True)

    print(lcs.contains("TCGCGTAACGCGT"))

    for key in lcs.base:
        # print(type(lcs.base[key]))

        if type(lcs.base[key]) == Block:
            print('Key: ', key, ' acc: ', lcs.base[key].values())

        elif type(lcs.base[key]) == tuple:
            print('Key: ', key, ' acc: ', lcs.base[key][0].values(), ' SNP: ', lcs.base[key][1])

    print(lcs)


if __name__ == "__main__":
    test_LCSHash()
