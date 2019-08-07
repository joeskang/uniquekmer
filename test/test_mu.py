import sys
import os
sys.path.append(os.getcwd().split('uniquekmer')[0] + '/uniquekmer/test')
sys.path.append(os.getcwd().split('uniquekmer')[0] + '/uniquekmer/src')
from kasai import kasai, inverse1
import numpy as np
import file_manager as fm
from tqdm import tqdm, trange



def test_validity():
    # assumes that mu's have already been json'ed
    mu = fm.unjson_it("../src/c22_mu")
    geno = fm.read_unambiguous("../data/22.fa")

    myd = {}

    for key in tqdm(mu, desc="checking uniqueness"):
        seq = geno[int(key):mu[key] + 1]
        assert seq not in myd
        myd[seq] = 1
        if key == "700":
            myd["ACGT"]


if __name__ == "__main__":
    test_validity()

