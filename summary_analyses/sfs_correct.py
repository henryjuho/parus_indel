from __future__ import print_function
import sys
import itertools
sys.path.append('..')
from anavar_analyses.sel_vs_neu_anavar import sfs2counts


def correct_sfs(sfs_i, sfs_d, e_i=0.0110086484429, e_d=0.0166354937984):

    """
    uses the model estimates of polarisation error to correct the sfs
    :param sfs_i: list
    :param sfs_d: list
    :param e_i: float
    :param e_d: float
    :return: list
    """

    # convert to counts
    sfs_i = sfs2counts(sfs_i, 20)
    sfs_d = sfs2counts(sfs_d, 20)
    freq_keys = [y/20.0 for y in range(1, 20)]
    corrected_i = []
    corrected_d = []

    for i in range(0, len(sfs_d)):

        n_i = sfs_i[i] - (sfs_i[i] * e_i) + (sfs_d[-(i+1)] * e_d)
        n_d = sfs_d[i] - (sfs_d[i] * e_d) + (sfs_i[-(i+1)] * e_i)

        freq = freq_keys[i]

        correct_i = [freq for f in range(0, int(round(n_i)))]
        correct_d = [freq for f in range(0, int(round(n_d)))]

        corrected_i += correct_i
        corrected_d += correct_d

    return corrected_i, corrected_d


def debug():

    """
    checks function for sfs correction using NC data
    :return:
    """

    sfs_i = [80139, 34617, 21586, 15479, 11971, 9639, 8173, 6929, 6160, 5589, 5048, 4617, 4401, 4260, 4097, 3977, 4313, 4514, 5443]
    sfs_d = [113865, 46173, 28183, 19922, 14912, 11503, 9731, 8355, 7215, 6485, 5970, 5426, 5120, 4947, 4753, 4497, 4480, 4677, 5005]

    n = 20
    e_i = 0.05
    e_d = 0.01

    sfs_i_freqs = [[(i+1)/20.0 for x in range(0, sfs_i[i])] for i in range(0, len(sfs_i))]
    sfs_d_freqs = [[(i+1)/20.0 for x in range(0, sfs_d[i])] for i in range(0, len(sfs_d))]

    sfs_i_one = 80139 - (80139*e_i) + (5005*e_d)
    sfs_i_end = 5443 - (5443*e_i) + (113865*e_d)

    new_ins, new_del = correct_sfs(itertools.chain(*sfs_i_freqs), itertools.chain(*sfs_d_freqs), e_i, e_d)

    assert round(sfs_i_one) == sfs2counts(new_ins, n)[0]
    assert round(sfs_i_end) == sfs2counts(new_ins, n)[-1]

    print('original reconstructed ins     : ', sfs2counts(itertools.chain(*sfs_i_freqs), n))
    print('corrected sfs insertions       : ', sfs2counts(new_ins, n))
    print('corrected singletons insertions: ', sfs_i_one)
    print('corrected n-1 freq insertions  : ', sfs_i_end)
