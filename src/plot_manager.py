from collections import deque
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import freedman_bin_width
from math import ceil
import sys
from Bio import motifs
from Bio.Seq import Seq
sys.path.append('../src')
from dict_manager import sort_dict_to_numpy, count_dict, count_kmer
from file_manager import read_kmer_custom_file
def draw_lcp_bar(in_dict: dict, top:int):


    """
        draws bar graph with in_dict values
        in_dict: sa_lcp_dict
    :param in_dict:
    :return:
    """

    count = count_dict(in_dict,top)

    keys, values = sort_dict_to_numpy(count)

    #find smallest and greatest counts
    new_v = sorted(values)


    #############
    print(keys)
    #############


    fig, ax = plt.subplots()

    # adjust bar width as necessary
    bar_width = 0.7

    opacity = 0.4

    error_config = {'ecolor':'0.3'}

    rects = ax.bar(x=np.array(keys), height=values, width=bar_width, alpha=opacity, color='b')

    x_max = keys[-1]

    ax.set_xlabel('Length of LCP')
    ax.set_ylabel('Frequency')
    #ax.set_xticks(keys)
    #ax.set_xlim(0,x_max+1)
    #ax.set_xticklabels(keys)
    ax.set_title('Bar Graph of LCP')
    x_ax = ax.xaxis


    ###############
    print("x_max = ", str(x_max))
    ###############

    # set axes: must include min and max of both x and y
    #ax.axis([0, x_max+1, new_v[0], new_v[-1] + 1 ])
    x_ax.get_ticklines(minor=True)


    # turn on grids
    plt.grid(True)

    # prevent clipping of y-label
    fig.tight_layout()

    plt.show()

    return

def draw_histo_list(myl: list):
    """
        draws histogram with in_dict values
        in_dict: full dictionary with sa address and corresponding lcp value
    :param in_dict:
    :return:
    """

    full = deque(sorted(myl))
    '''
    for key in in_dict.keys():
        for i in range(in_dict[key]):
            full.append(key)
    '''
    bin_width = freedman_bin_width(data=full)


    ############
    print("width of bars: ", str(bin_width))
    ############


    x_min = full[0]
    x_max = full[-1]

    ############
    print("Greatest lcp = ", str(x_max))
    print("Lowest lcp = ", str(x_min))
    ############

    n_bins = ceil( (x_max - x_min) / bin_width)

    ###########
    print("Number of bins: ", str(n_bins))
    ###########

    fig, ax = plt.subplots()
    x_ax = ax.xaxis

    # draw histogram
    n, bins, patches = plt.hist(x=full, bins=n_bins,density=1)


    #plt.axis([0, x_max+1, 0, 1])
    ax.set_xlabel('LCP')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram of LCP')

    plt.grid(True)
    plt.tight_layout()
    plt.show()

def draw_lcp_histo(in_dict: dict):
    """
        draws histogram with in_dict values
        in_dict: full dictionary with sa address and corresponding lcp value
    :param in_dict:
    :return:
    """

    full = deque(sorted(in_dict.values()))
    '''
    for key in in_dict.keys():
        for i in range(in_dict[key]):
            full.append(key)
    '''
    bin_width = freedman_bin_width(data=full)


    ############
    print("width of bars: ", str(bin_width))
    ############


    x_min = full[0]
    x_max = full[-1]

    ############
    print("Greatest lcp = ", str(x_max))
    print("Lowest lcp = ", str(x_min))
    ############

    n_bins = ceil( (x_max - x_min) / bin_width)

    ###########
    print("Number of bins: ", str(n_bins))
    ###########

    fig, ax = plt.subplots()
    x_ax = ax.xaxis

    # draw histogram
    n, bins, patches = plt.hist(x=full, bins=n_bins,density=1)


    plt.axis([0, x_max+1, 0, 1])
    ax.set_xlabel('LCP')
    ax.set_ylabel('Frequency')
    ax.set_title('Histogram of LCP')

    plt.grid(True)
    plt.tight_layout()
    plt.show()

def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False 
    otherwise.

    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.

    Returns:
    --------
        mask : A numobservations-length boolean array.

    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor. 
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def draw_bar_list(myl: list):

    """
        draws bar graph with in_dict values
        in_dict: sa_lcp_dict
    :param in_dict:
    :return:
    """
    count = {}
    for _ in range(len(myl)):
        curr = myl[_]
        count[curr] = count[curr] + 1 if curr in count else 1

    keys, values = sort_dict_to_numpy(count)

    values = values[~is_outlier(values)]

    #find smallest and greatest counts
    new_v = sorted(values)


    #############
    print(keys)
    #############

    # set the number of groups (along x-axis) equal to the number of keys in the count dictionary
    n_groups = len(keys)

    fig, ax = plt.subplots()
    x_ax = ax.xaxis

    index = np.arange(n_groups)

    # adjust bar width as necessary
    bar_width = 0.7

    opacity = 0.4

    error_config = {'ecolor':'0.3'}

    rects = ax.bar(x=index, height=values, width=bar_width, alpha=opacity, color='b')

    ax.set_xlabel('Size of Clusters')
    ax.set_ylabel('Frequency')
    #ax.set_xticks(keys)
    #ax.set_xticklabels(keys)
    ax.set_title('Bar Graph of k-mer clusters')
    x_ax.get_ticklines(minor=True)
    x_max = keys[-1]

    ###############
    print("x_max = ", str(x_max))
    ###############

    # set axes: must include min and max of both x and y
    #plt.axis([0, x_max+1, new_v[0], new_v[-1] + 1 ])

    # turn on grids
    plt.grid(True)

    # prevent clipping of y-label
    fig.tight_layout()
    plt.show()

    return


def draw_kmer_bar(in_dict: dict,bot:int, top:int):

    """
        draws bar graph with in_dict values
        in_dict: sa_lcp_dict
    :param in_dict:
    :return:
    """

    count = count_kmer(in_dict=in_dict,bot=bot, top=top)

    keys, values = sort_dict_to_numpy(count)

    #find smallest and greatest counts
    new_v = sorted(values)


    #############
    print(keys)
    #############

    # set the number of groups (along x-axis) equal to the number of keys in the count dictionary
    n_groups = len(keys)

    fig, ax = plt.subplots()
    x_ax = ax.xaxis

    index = np.arange(n_groups)

    # adjust bar width as necessary
    bar_width = 0.7

    opacity = 0.4

    error_config = {'ecolor':'0.3'}

    rects = ax.bar(x=index, height=values, width=bar_width, alpha=opacity, color='b')

    ax.set_xlabel('Length of k-mer')
    ax.set_ylabel('Frequency')
    #ax.set_xticks(keys)
    #ax.set_xticklabels(keys)
    ax.set_title('Bar Graph of K-mer')
    x_ax.get_ticklines(minor=True)
    x_max = keys[-1]

    ###############
    print("x_max = ", str(x_max))
    ###############

    # set axes: must include min and max of both x and y
    #plt.axis([0, x_max+1, new_v[0], new_v[-1] + 1 ])

    # turn on grids
    plt.grid(True)

    # prevent clipping of y-label
    fig.tight_layout()
    plt.show()

    return


def draw_motif(in_file:str, sequences=''):

    # ++++++++++++++++++
    # read custom file
    # ++++++++++++++++++
    sequences = read_kmer_custom_file(in_file)

    # ++++++++++++++++++
    # motif creation
    # ++++++++++++++++++

    # convert sequences to Bio.Seq
    instances = []
    for seq in sequences:
        instances.append(Seq(seq))

    assert instances
    # create motif
    m = motifs.create(instances)
    print(m.counts)
    m.weblogo("c22_21-mer.png")
    pwm = m.counts.normalize(pseudocounts=0.5)

    return


def test():

    draw_motif(in_file='../output/0403_c22_21-mer')
    return


if __name__ == '__main__':
    test()
