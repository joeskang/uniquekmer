# truncate sequence to only include nucleotides "A","C","G", or "T"
"""
	read from file and truncate to only A, C, G, or T
	return list with the strings
"""

def trun(myfile: str, start: int, end: int):

    myfile.seek(start)
    return_str = ''

    for i in range(end - start):
        line = myfile.readline()
        if line[0:1] == ">":
            pass
        else:

            newline = ""

            for char in line:
                if char != "A" and char != "C" and char != "G" and char != "T":
                    pass
                else:
                    newline = newline + char
        return_str += newline

    return return_str


# outfile.write(line[:-1])#get rid of \n
def toBinary(string:str):
    # given that the alphabet = { 'A', 'C', 'G', 'T' }
    alphabet = {

        'A': '00',
        'C': '01',
        'G': '10',
        'T': '11'
    }

    new_str = ''
    for char in string:
        new_str += alphabet[char]


def only_nucleotide(infile: str, outfile: str):
    """
        dispose of chr info
    :param infile:
    :return:
    """
    inf = open(infile, 'r')
    outf = open(outfile, 'w')
    for line in inf:
        if line[0] == '>':
            continue
        else:
            outf.write(line)


if __name__ == '__main__':
    inputfilename = input("Enter the name of file to be read: ")
    outputfilename = input("Enter the name of file to be written: ")

    only_nucleotide(infile=inputfilename, outfile=outputfilename)
