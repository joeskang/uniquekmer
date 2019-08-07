"""

"""


from tqdm import trange, tqdm
from file_manager import append_file_name, reads, unjson_it, json_it
from collections import deque, OrderedDict
from experimental.LCSHash import LCSHash
from KmerTree import KmerTree


def parse_SNP(filename:str):
    # read from dbSNP file
    with open(filename,'r') as file:
        for line in file:
            myl = line.split()

            # check if SNP
            if myl[1] == 'single':
                # column 4 is the address
                # column 5 is the SNP
                yield int(myl[3]), myl[4]


def print_simple_file(myl:list, filename:str):
    with open(filename, 'w') as file:
        for item in myl:
            file.write(str(item) + "\n")


def sort_SNP(input:dict):
    keys = list(input.keys())
    keys.sort()

    for k in tqdm(keys, desc="Sorting the dict: "):
        yield k, input[k]


def check_numerical(input:dict):
    past = 0
    for k in tqdm(input, desc="checking if sorted in numerical order: "):
        assert k < past
        past = k


def validate_addr(addr:list, mu:dict, snp:list):
    # run through addr list from SNP file and see which match to uniquemers found
    # assumes that addr is sorted
    addr = deque(addr)
    snp = deque(snp)
    a = addr.popleft()
    s = snp.popleft()

    for key in tqdm(mu.keys(), desc="validating SNPs: "):

        if a < int(key):
            continue
            # print(a)
            # raise Exception("Error, SNP address less than uniquemer address. Possible that SNP addresses not sorted?")

        while int(key) <= a <= mu[key]:
            yield a, s
            a = addr.popleft()
            s = snp.popleft()


def output_reads(filename:str,snp:list, addr:list, geno:str, trues=[], uniques=[], tops=[], top=100, bot=20):
    assert len(addr) == len(snp)
    with open(filename, 'w') as file:
        # for each address found in SNP file
        for a in trange(len(addr)):
            end = addr[a] + top - 1
            start = addr[a] - (top - 1)
            # get sequence that corresponds to i - 99 to i + 99
            # insert snp
            seq = geno[start:start+top-1]+snp[a]+geno[start+top:end]
            for t in range(top):
                for j in range(top-bot+1):
                    kmer = seq[t:bot+j]
                    if "N" not in kmer:
                        file.write(kmer)

    return


####################
# if minimal unique was the short k-mer such that any substring of the k-mer is not a unique k-mer,
# then longest common un-uniques are substrings of such k-mers that are n-1 length.
# theoretically, two LCU can be produced per MU
# using LCU can compress


def use_LCSTree(mu:dict, geno:str):
    lcs = LCSHash()

    for key in tqdm(mu, desc="building LCS table: "):
        kmer = geno[int(key):int(mu[key]) + 1]
        if len(kmer) > 2:
            if 'N' in kmer:
                raise Exception("N found in kmer: ", kmer, " at address: ", key, " to: ", mu[key])
            lcs.insert(kmer)

        else:
            # ignore "unique" k-mers that have less than length 3 for now
            pass
    return lcs


def SNP_compare(lcs:LCSHash, snp:dict, mu:dict, geno:str):
    count = 0
    for a in snp:
        for start in mu:
            end = int(mu[start])
            start = int(start)
            if a < start:
                break
            elif start <= a <= end:
                new_mer = geno[start:a] + snp[a] + geno[a:end + 1]
                if lcs.contains(new_mer):
                    continue
                    # TODO: *** need to catch new-mers that have been inserted multiple times ***
                else:
                    lcs.insert(new_mer, is_snp=True)
                    count += 1

            elif end < a:
                break

    return lcs, count


def use_KmerTree(mu, geno):
    kt = KmerTree()

    for key in tqdm(mu, desc="building KmerTree: "):
        kmer = geno[int(key):int(mu[key]) + 1]
        if len(kmer) > 2:
            if 'N' in kmer:
                raise Exception("N found in kmer: " + kmer + " at address: " + str(key) + " to: " + str(mu[key]))
            kt.add(kmer)

        else:
            pass

    return kt

def driver():

    addr = []
    snp = []

    for tup in parse_SNP(filename=append_file_name("data/22.snp")):
        addr.append(tup[0])
        snp.append(tup[1])

    with open("../output/ez_c22_SNP_pos", "w") as file:
        with open ("../output/ez_c22_SNP_char", 'w') as file_2:
            for tup in parse_SNP(filename=append_file_name("data/22.snp")):
                file.write(str(tup[0]) + '\n')
                file_2.write(str(tup[1]) + '\n')


if __name__ == '__main__':
    driver()

