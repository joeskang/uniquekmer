import sys
import os
from tqdm import tqdm


sys.path.append(os.getcwd().split('uniquekmer')[0] + 'uniquekmer/src')

import file_manager as fm


filename = '../22_json_default_dict'
print('reading jsoned dict')
d = fm.unjson_it(filename)
print('read!')

print('reading genome')
genome = fm.reads('../data/22.fa')
print('read!')

ambs=0
tops=0

# checking how many tops are less than 100
for sa in tqdm(d, desc='checking dict'):
    top = int(d[sa][1])
    sa = int(sa)
    string = genome[sa:sa+top]
    if 'N' in string:
        ambs+=1

    if top < 100:
        tops += 1


print('number of sequences with ambs: ',ambs)
print('number of tops less than 100: ', tops)
