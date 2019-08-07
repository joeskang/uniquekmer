
from collections import deque
from tqdm import trange, tqdm

from file_manager import write_array_to_byte, reads, read_unambiguous, append_file_name

'''
	suffix array generator for relatively small strings
	use for debugging
'''
def naive_main():

    # simple script to test suffix array

    in_str = input("Enter a string ")
    naive_SA(in_str)
    return

def naive_SA(string: str, verbose=False):
    """
		inputs a string and outputs deque for the suffix array and string
	:param string:
	:return:
	"""

    str_dict = {}
    sa_array = deque()

    # generate suffixes
    for i in trange(len(string), desc='Creating suffixes: '):
        temp = string[i:]

        # store the position
        str_dict[temp] = i

    # generate deque with actual suffixes
    suffixes = deque(sorted(str_dict.keys()))

    # generate sa_array with position values to corresponding suffixes
    # also, generate last column, L
    L = ""
    pbar = tqdm(total=len(suffixes), desc='Generating suffix array: ')
    while suffixes:
        index = str_dict[suffixes.popleft()]
        sa_array.append(index)
        L += string[index - 1]
        pbar.update(1)

    if verbose:
        print("Suffix array: ", sa_array)
        print("The last column is: ", L)

    return sa_array, L


def SA_file(filename:str):
    print('Reading file: ')
    sequence = read_unambiguous(filename)
    print('File read!\nCreating Suffix Array Naively:')

    s_array, _ = naive_SA(sequence)

    length = len(s_array)
    s_array.insert(0,length)

    print('Suffix Array Created!\nWriting to file: ')

    write_array_to_byte(filename='fake_genome_sa',byte_arr=s_array)

    return



if __name__ == '__main__':
    SA_file(append_file_name('test/fake_genome'))
