# -*- coding: utf-8 -*-
""" Parse SAM/BAM alignment files

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import os
from collections import Counter
import logging as lg

import pysam

from telescope.utils.calignment import AlignedPair
# from ._alignment import AlignedPair
from . import model
from . import BIG_INT

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

CODES = [
    ('SU', 'single_unmapped'),
    ('SM', 'single_mapped'),
    ('PU', 'pair_unmapped'),
    ('PM', 'pair_mapped'),
    ('PX', 'pair_mixed'),
    ('PX*', 'pair_mixed_unmapped'),
]

CODE_INT = {t[0]:i for i,t in enumerate(CODES)}

def readkey(aln):
    ''' Key that is unique for alignment '''
    return (aln.query_name, aln.is_read1,
            aln.reference_id, aln.reference_start,
            aln.next_reference_id, aln.next_reference_start,
            abs(aln.template_length))


def matekey(aln):
    return (aln.query_name, not aln.is_read1,
            aln.next_reference_id, aln.next_reference_start,
            aln.reference_id, aln.reference_start,
            abs(aln.template_length))

def mate_before(aln):
    """ Check if mate is before (to the left) of aln

    Alignment A is before alignment B if A appears before B in a sorted BAM
    alignment. If A and B are on different chromosomes, the reference ID is
    compared.

    Args:
        aln (:obj:`pysam.AlignedSegment`): An aligned segment

    Returns:
        bool: True if alignment's mate is before, False otherwise

    """
    if aln.next_reference_id == aln.reference_id:
        return aln.next_reference_start < aln.reference_start
    return aln.next_reference_id < aln.reference_id

def mate_after(aln):
    """ Check if mate is after (to the right) of aln

    Alignment A is after alignment B if A appears after B in a sorted BAM
    alignment. If A and B are on different chromosomes, the reference ID is
    compared.

    Args:
        aln (:obj:`pysam.AlignedSegment`): An aligned segment

    Returns:
        bool: True if alignment's mate is after, False otherwise

    """
    if aln.next_reference_id == aln.reference_id:
        return aln.next_reference_start > aln.reference_start
    return aln.next_reference_id > aln.reference_id

def mate_same(aln):
    """ Check if mate has same start position

    Args:
        aln (:obj:`pysam.AlignedSegment`): An aligned segment

    Returns:
        bool: True if alignment's mate has same start position, False otherwise

    """
    return aln.next_reference_id == aln.reference_id and \
           aln.next_reference_start == aln.reference_start

def mate_in_region(aln, regtup):
    """ Check if mate is found within region

    Return True if mate is found within region or region is None

    Args:
        aln (:obj:`pysam.AlignedSegment`): An aligned segment
        regtup (:tuple: (chrom, start, end)): Region
    Returns:
        bool: True if mate is within region

    """
    if regtup is None: return True
    return aln.next_reference_id == regtup[0] and \
           regtup[1] < aln.next_reference_start < regtup[2]


""" Sequential read"""
def fetch_bundle(samfile, **kwargs):
    """ Iterate over alignment over reads with same ID """
    samiter = samfile.fetch(**kwargs)
    bundle = [ next(samiter) ]
    for aln in samiter:
        if aln.query_name == bundle[0].query_name:
            bundle.append(aln)
        else:
            yield bundle
            bundle = [aln]
    yield bundle


def pair_bundle(alniter):
    readcache = {}
    for aln in alniter:
        if not aln.is_paired:
            yield AlignedPair(aln)
        else:
            mate = readcache.pop(matekey(aln), None)
            if mate is not None:  # Mate found in cache
                if aln.is_read1:
                    yield AlignedPair(aln, mate)
                else:
                    yield AlignedPair(mate, aln)
            else:  # Mate not found in cache
                readcache[readkey(aln)] = aln

    # Yield the remaining reads in the cache as unpaired
    for aln in readcache.values():
        yield AlignedPair(aln)


def fetch_fragments_seq(samfile, **kwargs):
    b_iter = fetch_bundle(samfile, **kwargs)
    for alns in b_iter:
        if not alns[0].is_paired:
            _code = CODE_INT['SU'] if alns[0].is_unmapped else CODE_INT['SM']
            yield (_code, [AlignedPair(a) for a in alns])
        else:
            if alns[0].is_proper_pair:
                yield (CODE_INT['PM'], list(pair_bundle(alns)))
            else:
                if len(alns) == 2 and all(a.is_unmapped for a in alns):
                    yield (CODE_INT['PU'], [AlignedPair(alns[0], alns[1]), ])
                else:
                    yield (CODE_INT['PX'], [AlignedPair(a) for a in alns])

""" Parallel Read """

def fetch_pairs_sorted(alniter, regtup=None):
    readcache = {}
    for aln in alniter:
        if not aln.is_paired:
            _code = CODE_INT['SU'] if aln.is_unmapped else CODE_INT['SM']
            yield (_code, AlignedPair(aln))
        else:
            if aln.is_proper_pair:
                mate = readcache.pop(matekey(aln), None)
                if mate:
                    if aln.is_read1:
                        yield (CODE_INT['PM'], AlignedPair(aln, mate))
                    else:
                        yield (CODE_INT['PM'], AlignedPair(mate, aln))
                else:
                    readcache[readkey(aln)] = aln
            else:
                assert not (aln.is_unmapped and aln.mate_is_unmapped), 'Found unmapped pair in sorted BAM file'
                _code = CODE_INT['PX*'] if aln.is_unmapped else CODE_INT['PX']
                yield (_code, AlignedPair(aln))

    for aln in readcache.values():
        #TODO: deal with these somehow
        yield ('cached', AlignedPair(aln))

def fetch_region(samfile, annotation, opts, region):
    lg.info('processing {}:{}-{}'.format(*region))

    _nfkey = opts['no_feature_key']
    _omode, _othresh = opts['overlap_mode'], opts['overlap_threshold']
    _tempdir = opts['tempdir']

    assign = model.Assigner(annotation, _nfkey, _omode, _othresh).assign_func()

    _minAS, _maxAS = BIG_INT, -BIG_INT
    _unaligned = 0

    mfile = os.path.join(_tempdir, 'tmp_map.{}.{}.{}.txt'.format(*region))

    fh = open(mfile, 'w')
    with pysam.AlignmentFile(samfile) as sf:
        samiter = sf.fetch(*region, multiple_iterators=True)
        regtup = (sf.get_tid(region[0]), region[1], region[2])
        for ci, aln in fetch_pairs_sorted(samiter, regtup):
            if aln.is_unmapped:
                assert CODES[ci][0] == 'PX*'
                _unaligned += 1
                continue

            m = (ci, aln.query_id, assign(aln), aln.alnscore, aln.alnlen)
            _minAS = min(_minAS, m[3])
            _maxAS = max(_maxAS, m[3])
            print('\t'.join(map(str, m)), file=fh)

    fh.close()
    return mfile, (_minAS, _maxAS), _unaligned
