__author__ = 'bendall'

import numpy as np

def phred(f):
    """ Calculate phred quality score for given error probability
    :param f: Error probability (float)
    :return:  Phred score (int)
    """
    return int(round(-10 * np.log10(1 - f))) if f < 1.0 else 255

def eprob(q):
    """ Calculate error probability for given phred quality score
    :param q: Phred score (int)
    :return:  Error probability (int)
    """
    return 1 - (10**(float(q) / -10))

def format_minutes(seconds):
    mins = seconds / 60
    secs = seconds % 60
    return '%d minutes and %d secs' % (mins,secs)
