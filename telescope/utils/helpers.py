# -*- coding: utf-8 -*-
""" Helper functions
"""
from __future__ import division

from past.utils import old_div
import numpy as np
from itertools import zip_longest

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
    return 1 - (10**(old_div(float(Q), -10)))


def format_minutes(seconds):
    mins = old_div(seconds, 60)
    secs = seconds % 60
    return '%d minutes and %d secs' % (mins,secs)


def merge_blocks(ivs, dist=0):
    """ Merge blocks

    Args:
        ivs (list): List of intervals. Each interval is represented by a tuple of
            integers (start, end) where end > start.
        dist (int): Distance between intervals to be merged. Setting dist=1 will merge
            adjacent intervals

    Returns:
        list: Merged list of intervals

    Examples:
        >>> merge_interval_list([])
        []
        >>> merge_interval_list([(1,10)])
        [(1, 10)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)])
        [(1, 3), (4, 9), (10, 14)]
        >>> merge_interval_list([(4, 9), (10, 14), (1, 3)], dist=1)
        [(1, 14)]
    """
    if len(ivs)<= 1: return ivs
    ivs.sort(key=lambda x:x[0])
    ret = [ivs[0]]
    for iv in ivs[1:]:
        if iv[0] - ret[-1][1] > dist:
            ret.append(iv)
        else:
           ret[-1] = (ret[-1][0], max(iv[1],ret[-1][1]))
    return ret

class GenomeRegion:
    """
    """
    def __init__(self, chrom=None, start=None, end=None, region=None):
        if region is not None:
            m = re.match('(\w+):(\d+)-(\d+)', region)
            self.chrom = m.group(1)
            self.start, self.end = int(m.group(2)), int(m.group(3))
        else:
            self.chrom = chrom
            if start is None or end is None:
                self.start = self.end = None
            else:
                self.start, self.end = int(start), int(end)
                if self.end < self.start:
                    self.start, self.end = self.end, self.start

    def contains(self, chrom, pos):
        if self.chrom is None: return True
        if self.start is None: return chrom == self.chrom
        return self.chrom == chrom and self.start <= int(pos) <= self.end

    def __str__(self):
        if self.chrom is None: return 'genome'
        if self.start is None: return '%s' % self.chrom
        return '%s:%d-%d' % (self.chrom, self.start, self.end)


def region_iter(refs, lengths, winsize=1e7, overlap=0):
    winsize, overlap = map(int, (winsize, overlap))
    for ref,reflen in zip(refs,lengths):
        for i in range(0, reflen, winsize):
            regmin = max(0, i-overlap)
            regmax = min(i+winsize+overlap, reflen)
            yield (ref, regmin, regmax)


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

def str2int(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s
