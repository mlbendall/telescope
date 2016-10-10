# -*- coding: utf-8 -*-
""" Helper functions
"""

import numpy as np

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


def phred(P):
    """ Calculate phred quality score for given probability/accuracy value

    Phred scores (Q) are logarithmically related to probabilities (P):

    Q = -10 * log10(1-P)

    Args:
        P (float): Probability between 0 and 1, inclusive

    Returns:
        int: Phred score

    Examples:
        >>> phred(0.9)
        10
        >>> phred(0.999999)
        60
        >>> phred(0)
        0
        >>> phred(1)
        255
    """
    return int(round(-10 * np.log10(1 - P))) if P < 1.0 else 255


def eprob(Q):
    """ Calculate probability/accuracy value for given phred quality score

    Probabilities (P) are logarithmically related to phred scores(Q):

    P = 1 - 10 ^ (-Q / 10)

    Args:
        Q (int): Phred score

    Returns:
        float: probability

    Examples:
        >>> eprob(10)
        0.9
        >>> eprob(60)
        0.999999
        >>> eprob(0)
        0.0
        >>> eprob(255)
        1.0
        >>> eprob(ord('@')-33)
        0.9992056717652757
    """
    return 1 - (10**(float(Q) / -10))


def format_minutes(seconds):
    mins = seconds / 60
    secs = seconds % 60
    return '%d minutes and %d secs' % (mins,secs)
