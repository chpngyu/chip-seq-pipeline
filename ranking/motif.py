# -*- coding: utf-8 -*-
"""
Created on Fri Feb 23 16:26:38 2018.

@author: Chun-Ping Yu

functions for motif (Position Weight matrix; PWM) management
"""

import numpy as np
from math import log, log2
from scipy.spatial.distance import correlation


def information_content(_pwm):
    """
    Calculate information content for a given PWM, and the correlation of small size of the PWM is used.

    Parameters
    ----------
    _pwm : matrix (a list of a list)
        position weight matrix of a motif.

    Returns
    -------
    __IC : float
        information content.
    """
    __IC = 0
    # set an epsilon to correct small sample size
    epsilon = (4.0 - 1.0) / (2.0 * _pwm.shape[0]) / log(2)
    for position_prob in _pwm:
        H = 0.0
        for prob in position_prob:
            if prob > 0:
                H -= prob * log(prob, 2)
        __IC += 2.0 - (H + epsilon)
    return __IC


def IC(_pwm):
    """Calculate IC value. it's same to the function information_content, but no small size correlation."""
    _ic = []
    for u in _pwm:
        __ic = 0
        for v in u:
            if v != 0.0:
                __ic += v * log2(v)
        _ic.append(2.0 + __ic)
    return _ic


def reverse_complement(_pwm):
    """Take reverse complement for PWM."""
    return np.flip(_pwm)


def normalize_pwm(_pwm):
    """Normalize row values into one."""
    _pred = np.zeros(_pwm.shape)
    _pwm[_pwm < 0] = 0
    for _idx, row in enumerate(_pwm):
        unit = 1.0 / np.sum(row)
        _pred[_idx] = row * unit
    return _pred


def norm(matrix, window_size):
    """
    Calculate the maximum norm value for a matrix in the length of window.

    Parameters
    ----------
    matrix : np.array
        a matrix to calculate its norm.
    window_size : int
        length for matrix norm.

    Returns
    -------
    max_query_idx : int
        the start site that has maximum norm value.
    max_subject_idx : int
        DESCRIPTION.
    max_score : float
        the maximum norm value.

    """
    query_w, subject_w = matrix.shape
    max_score = 0
    max_query_idx = 0
    max_subject_idx = 0
    for query_idx in range(0, query_w - window_size + 1):
        for subject_idx in range(0, subject_w - window_size + 1):
            _score = np.diag(matrix[query_idx:query_idx + window_size, subject_idx:subject_idx + window_size]).sum()
            if _score > max_score:
                max_query_idx = query_idx
                max_subject_idx = subject_idx
                max_score = _score
    return (max_query_idx, max_subject_idx), max_score


def align(query, subject, window_size):
    """Align one PWM (query) with the other PWM (subject)."""
    subject_w = subject.shape[0]
    query_w = query.shape[0]
    if window_size > min(subject_w, query_w):
        message = f'The window size ({window_size:d}) has to be less than the width of either subject ' \
                  f'({subject_w:d}) or query ({query_w:d}) '
        raise RuntimeError(message)
    position_score = np.matmul(query, subject.transpose())
    align_pos, max_score = norm(position_score, window_size)
    query_rc = reverse_complement(query)
    align_pos_rc, max_score_rc = norm(np.matmul(query_rc, subject.transpose()), window_size)
    new_query = np.empty(subject.shape)
    new_query.fill(0.25)
    if max_score_rc > max_score:
        shift = align_pos_rc[0] - align_pos_rc[1]
        start_new = max(0, -shift)
        end_new = min(subject_w, query_w - shift)
        start_q = max(0, shift)
        end_q = min(query_w, end_new + shift)
        new_query[start_new:end_new] = query_rc[start_q:end_q]
        align_score = max_score_rc
        _strand = 'reverse'
    else:
        shift = align_pos[0] - align_pos[1]
        start_new = max(0, -shift)
        end_new = min(subject_w, query_w - shift)
        start_q = max(0, shift)
        end_q = min(query_w, end_new + shift)
        new_query[start_new:end_new] = query[start_q: end_q]
        align_score = max_score
        _strand = 'forward'
    return new_query, align_score, _strand


def trim(_pwm, min_prob=0.3, isIC=True):
    """Trim one or more positions of PWM at the two ends, if each max probability (or IC) < min_prob."""
    if isIC:
        max_prob = IC(_pwm)
    else:
        max_prob = _pwm.max(axis=1)
    start_pos = 0
    for index, prob in enumerate(max_prob):
        if prob >= min_prob:
            start_pos = index
            break
    end_pos = 0
    for index in range(len(max_prob) - 1, -1, -1):
        if max_prob[index] >= min_prob:
            end_pos = index + 1
            break
    return _pwm[start_pos:end_pos]


def centralize(_pwm, take_rc=False, keep_size=None):
    """
    A given position weight matrix (PWM) is centralized to a new and wider PWM.

    (1) if the take_rc sets to True, the given pwm will be converted to reverse complement, and then
    (2) do nothing if keep_size is not set, or
    (3) if keep_size is set and smaller than the width of _pwm,
        only return a truncated pwm from its left-end (0-base) to the position of keep_size (not included), or
    (4) if keep_size >= the width of the _pwm, centralize the pwm in middle of a new and wider PWM.
        The new probabilities in the extended positions aer set to 1/4.

    Parameters
    ----------
    _pwm: ndarray (2D)
        the input PWM
    take_rc: bool
        whether to take reverse complement
    keep_size: int or None
        give a new width for a PWM to be created

    Returns
    -------
    ndarray (2D)
        a new PWM
    """
    if take_rc:
        _pwm = reverse_complement(_pwm)
    if keep_size:
        width, n_feature = _pwm.shape
        if width > keep_size:
            # drop out values at right-hand sides such that the width of pwm is equal to keep_size
            return _pwm[:keep_size]
        else:  # centralized the pwm
            new_pwm = np.empty((keep_size, n_feature))
            new_pwm.fill(0.25)
            new_start = int((keep_size - width) / 2)
            new_pwm[new_start:new_start + width] = _pwm
            return new_pwm
    return _pwm


def motif_correlation(motif1, motif2, has_rc=True):
    """
    Calculate pearson correlation value between two matrices. Return 0 if two are the same.

    precondition: m1.shape == m2.shape

    Parameters
    ----------
    motif1 : np array
        one matrix
    motif2 : np array
        the other matrix.
    has_rc : bool
        do or not do reverse complement
    Returns
    -------
    float
        correlation value, ranging from 0 (perfect correlated) to 2 (perfectly anti-correlated).

    """
    # correlation(x, x) == 0
    _x = correlation(motif1.flatten(), motif2.flatten())
    # taking reverse complement for m1 and m2
    if has_rc:
        m2_rc = motif2[::-1, ::-1]
        y = correlation(motif1.flatten(), m2_rc.flatten())
        if y < _x:
            return y, -1
    return _x, 1


def piece_correlation(motif1, motif2, has_rc=True):
    """
    Find the correlations between motif1 and motif2, where the comparison
    allowing motif2 can be is aligned with motif1 .

    Parameters
    ----------
    motif1: numpy.ndarray
        first motif used as a reference PWM.
    motif2: numpy.ndarray
        2nd motif used as a subject PWM.
    has_rc: boolean
        allow 2nd motif to be aligned with 1st motif in reverse complement

    Raises
    ------
    RuntimeError
        length of 2nd motif (m2) has to be shorter or equal to 1st motif (m1).

    Returns
    -------
     (numpy.ndarray, numpy.ndarray)
        1. correlation values for the alignments of motif2 with motif1 from 1st positions to last available position.
        2. which strand of motif2 is aligned: forward (+1) or reverse (-1)
    """
    width1, width2 = motif1.shape[0], motif2.shape[0]
    if width1 < width2:
        raise RuntimeError(
            'the length of 1st motif has to be longer or equal to the length of 2nd motif')
    res = [motif_correlation(motif1[start:start+width2,:], motif2, has_rc=has_rc) for start in range(0,width1-width2+1)]
    correlations = np.array([1-a for a, b in res])
    strands = np.array([b for a, b in res])
    return correlations, strands


def auto_correlation(motif1, motif2, fill=False):
    """
    Find the maximum correlation between motif1 and motif2.
    The motif2 is allows to be aligned with motif1 before doing the correlation.

    precondition: motif1.shape[0] >= motif2.shape[0]

    Parameters
    ----------
    motif1 : numpy.ndarray
        subject PWM with size of N*4.
    motif2 : numpy.ndarray
        target PWM with size of M*4, where N >= M.
    fill : bool
        whether to fill probabilities in the unaligned positions. default NO

    Returns
    -------
    (float, int, int)
        cos value. return 1 if motif1 == motif2
        strand (1 or -1). 1 (or -1) indicates that m2 is aligned by forward (or reverse) strand
        position (int). the position of motif1 that motif2 can be aligned to get the maximum correlation
    """
    _correlations, _aligned_strands = piece_correlation(motif1, motif2)
    _idx = np.argmax(_correlations)
    _aligned_strand = _aligned_strands[_idx]
    if not fill:
        return _correlations[_idx], _aligned_strand, _idx
    temp_motif2 = np.zeros(motif1.shape) + 0.25
    _end = _idx + motif2.shape[0]
    if _aligned_strand < 0:
        temp_motif2[_idx:_end,:] = reverse_complement(motif2)
    else:
        temp_motif2[_idx:_end, :] = motif2
    _correlation, _aligned_strand2 = motif_correlation(motif1, temp_motif2, has_rc=True)
    return 1-_correlation, _aligned_strand*_aligned_strand2, _idx


if __name__ == "__main__":
    print('example for manipulating the motifs')
    m1 = np.array([[1, 0, 0, 0], [0, 0, 0, 1], [0.9, 0, 0.1, 0], [0, 0.8, 0.2, 0]])
    m2 = np.abs(np.random.randn(10, 4))*0.1
    m2[3:3+m1.shape[0],:] += reverse_complement(m1)
    m2 = normalize_pwm(m2)
    print('m1:\n', m1.shape)
    print('m2:\n', m2.shape)

    print(auto_correlation(m2, m1))
    print(auto_correlation(m2, m1, fill=True))
