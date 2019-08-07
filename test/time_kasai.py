import sys
sys.path.append('../src')
from kasai import inverse1, inverse2, inverse3
import random
import time
from exception_manager import get_time
from tqdm import trange


def main():
    size = 1000000000

    s_array = random.sample(population=range(size), k=size)

    past = time.time()

    inverse1(s_array=s_array)
    print('inverse1 completed.')
    past = get_time(past=past)
    print()

    inverse2(s_array=s_array)
    print('inverse2 completed.')
    past = get_time(past=past)
    print()

    inverse3(s_array)
    print('inverse3 completed')
    past = get_time(past=past)
    print()


if __name__ == '__main__':
    main()