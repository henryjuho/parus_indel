from __future__ import print_function
import sys
import itertools
sys.path.append('..')
from anavar_analyses.sel_vs_neu_anavar import sfs2counts
import numpy as np


def correct_sfs(sfs_i, sfs_d, e_i=0.0110086484429, e_d=0.0166354937984, n=20):

    """
    uses the model estimates of polarisation error to correct the sfs
    :param sfs_i: list
    :param sfs_d: list
    :param e_i: float
    :param e_d: float
    :param n: int
    :return: list
    """

    # convert to counts
    sfs_i = sfs2counts(sfs_i, n)
    sfs_d = sfs2counts(sfs_d, n)
    freq_keys = [round(y/float(n), 3) for y in range(1, n)]
    corrected_i = []
    corrected_d = []

    # loop per frequency (forwards through sfs for ins and backwards for del)
    for i in range(0, len(sfs_d)):

        n_i, n_d = [sfs_i[i], sfs_d[-(i+1)]]  # get the frequency count pairs dels start at n-1 ins at 1

        # error matrix
        a = np.array([[1-e_i, e_d],
                      [e_i, 1-e_d]])

        # observed counts
        b = np.array([n_i, n_d])

        # corrected data
        i_true, d_true = np.linalg.solve(a, b)

        freq_i = freq_keys[i]
        freq_d = freq_keys[-(i+1)]

        correct_i = [freq_i for f in range(0, int(round(i_true)))]
        correct_d = [freq_d for f in range(0, int(round(d_true)))]

        corrected_i += correct_i
        corrected_d += correct_d

    return corrected_i, corrected_d


def add_error_to_sfs(sfs_i, sfs_d, e_i, e_d):

    """
    adds defined error on sfs
    :param sfs_i: list
    :param sfs_d: list
    :param e_i: float
    :param e_d: float
    :return: list
    """

    error_i = []
    error_d = []

    for i in range(0, len(sfs_d)):

        n_i = sfs_i[i] - (sfs_i[i] * e_i) + (sfs_d[-(i+1)] * e_d)
        n_d = sfs_d[i] - (sfs_d[i] * e_d) + (sfs_i[-(i+1)] * e_i)

        error_i.append(int(round(n_i)))
        error_d.append(int(round(n_d)))

    return error_i, error_d


def debug():

    """
    checks function for sfs correction using test data set
    :return: None
    """

    # 'true' data
    sfs_i_true = [200, 50, 30, 20, 5]
    sfs_d_true = [400, 80, 20, 10, 2]
    n = 6
    e_i = 0.05
    e_d = 0.01

    # error data - ie observed
    sfs_obs_i, sfs_obs_d = add_error_to_sfs(sfs_i_true, sfs_d_true, e_i, e_d)

    # estimate true

    sfs_i_freqs = [[round((i+1)/6.0, 3) for x in range(0, sfs_obs_i[i])] for i in range(0, len(sfs_obs_i))]
    sfs_d_freqs = [[round((i+1)/6.0, 3) for x in range(0, sfs_obs_d[i])] for i in range(0, len(sfs_obs_d))]
    #
    # sfs_i_one = 80139 - (80139*e_i) + (5005*e_d)
    # sfs_i_end = 5443 - (5443*e_i) + (113865*e_d)
    #
    new_ins, new_del = correct_sfs(itertools.chain(*sfs_i_freqs), itertools.chain(*sfs_d_freqs), e_i, e_d, n)

    print('# INSERTIONS #')
    print('True: ', sfs_i_true)
    print('Obs : ', sfs_obs_i)
    print('Est : ', sfs2counts(new_ins, 6))
    print('# DELETIONS #')
    print('True: ', sfs_d_true)
    print('Obs : ', sfs_obs_d)
    print('Est : ', sfs2counts(new_del, 6))

    assert sfs2counts(new_del, 6) == sfs_d_true
    assert sfs2counts(new_ins, 6) == sfs_i_true
