# -*- coding: utf-8 -*-
""" Parse SAM/BAM alignment files

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *
from collections import Counter

import pysam

from telescope.utils.calignment import AlignedPair
# from ._alignment import AlignedPair

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


def readkey(aln):
    ''' Key that is unique for alignment '''
    return (aln.query_name, aln.is_read1,
            aln.reference_id, aln.reference_start,
            aln.next_reference_id, aln.next_reference_start,
            abs(aln.template_length))


def matekey(aln):
    ''' Expected key of mate '''
    # if aln.reference_id == aln.next_reference_id and aln.reference_start == aln.next_reference_start:
    #     tlen = aln.template_length
    # else:
    #     tlen = -aln.template_length
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

    Alignment A is before alignment B if A appears before B in a sorted BAM
    alignment. If A and B are on different chromosomes, the reference ID is
    compared.

    Args:
        aln (:obj:`pysam.AlignedSegment`): An aligned segment

    Returns:
        bool: True if alignment's mate is before, False otherwise

    """
    if aln.next_reference_id == aln.reference_id:
        return aln.next_reference_start > aln.reference_start
    return aln.next_reference_id > aln.reference_id

def mate_same(aln):
    """ Check if mate has same start position

    Alignment A is before alignment B if A appears before B in a sorted BAM
    alignment. If A and B are on different chromosomes, the reference ID is
    compared.

    Args:
        aln (:obj:`pysam.AlignedSegment`): An aligned segment

    Returns:
        bool: True if alignment's mate is before, False otherwise

    """
    return aln.next_reference_id == aln.reference_id and \
           aln.next_reference_id == aln.reference_id

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

#
# def pair_alignments(alniter):
#     readcache = {}
#     for aln in alniter:
#         if not aln.is_paired:
#             yield AlignedPair(aln)
#         else:
#             mate = readcache.pop(matekey(aln), None)
#             if mate is not None:  # Mate found in cache
#                 if aln.is_read1:
#                     yield AlignedPair(aln, mate)
#                 else:
#                     yield AlignedPair(mate, aln)
#             else:  # Mate not found in cache
#                 readcache[readkey(aln)] = aln
#     # Yield the remaining reads in the cache as unpaired
#     for aln in readcache.values():
#         yield AlignedPair(aln)
#
#
# def organize_bundle(alns):
#     if not alns[0].is_paired:
#         assert all(not a.is_paired for a in alns), 'mismatch pair flag'
#         return [AlignedPair(a) for a in alns], None
#     else:
#         if len(alns) == 2 and alns[0].is_unmapped and alns[1].is_unmapped:
#             # Unmapped fragment
#             return [AlignedPair(alns[0], alns[1])], None
#         elif alns[0].is_proper_pair:
#             return list(pair_alignments(alns)), None
#         else:
#             return [AlignedPair(aln) for aln in alns], None
#             # ret1, ret2 = list(), list()
#             # for aln in alns:
#             #     if aln.is_read1:
#             #         ret1 += [AlignedPair(aln)]
#             #     else:
#             #         assert aln.is_read2
#             #         ret2 += [AlignedPair(aln)]
#             # return ret1, ret2
#
#
#
#
# def fetch_fragments(samfile, **kwargs):
#     biter = fetch_bundle(samfile, **kwargs)
#     for bundle in biter:
#         f1, f2 = organize_bundle(bundle)
#         yield f1
#         if f2 is not None:
#             print('here')
#             yield f2

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
            if alns[0].is_unmapped:
                assert len(alns) == 1
                yield 'SU', [AlignedPair(alns[0]), ]
            else:
                yield 'SM', [AlignedPair(a) for a in alns]
        else:
            if len(alns) == 2 and alns[0].is_unmapped and alns[1].is_unmapped:
                yield 'PU', [AlignedPair(alns[0], alns[1]), ]
            elif alns[0].is_proper_pair:
                yield 'PM', list(pair_bundle(alns))
            else:
                yield 'PX', [AlignedPair(a) for a in alns]


""" Parallel Read """
#
# def readkeyP(aln):
#     ''' Key that is unique for alignment '''
#     return (aln.query_name, aln.is_read1,
#             aln.reference_id, aln.reference_start,
#             aln.next_reference_id, aln.next_reference_start,
#             abs(aln.template_length))
#
#
# def matekeyP(aln):
#     ''' Expected key of mate '''
#     # if aln.reference_id == aln.next_reference_id and aln.reference_start == aln.next_reference_start:
#     #     tlen = aln.template_length
#     # else:
#     #     tlen = -aln.template_length
#     return (aln.query_name, not aln.is_read1,
#             aln.next_reference_id, aln.next_reference_start,
#             aln.reference_id, aln.reference_start,
#             abs(aln.template_length))

# READSTOPRINT = ["HWI-ST423:3:4:12572:24668:",]#  #  'HWI-ST423:3:4:8117:35965:']
# READSTOPRINT = ['HWI-ST423:3:4:16700:8733:',]
# READSTOPRINT = ['HWI-ST423:3:4:11856:57492:',]
READSTOPRINT = []
PRINTCACHE = True



def fetch_pairs_sorted(alniter, regtup=None):
    readcache = {}
    outcache = {}
    for aln in alniter:
        if not aln.is_paired:
            if aln.is_unmapped:
                yield ('SU', AlignedPair(aln))
            else:
                yield ('SM', AlignedPair(aln))
        else:
            if aln.is_proper_pair:
                mate = readcache.pop(matekey(aln), None)
                if mate:
                    if aln.is_read1:
                        yield ("PM", AlignedPair(aln, mate))
                    else:
                        yield ("PM", AlignedPair(mate, aln))
                else:
                    readcache[readkey(aln)] = aln
            else:
                if aln.is_unmapped:
                    yield ('PU', AlignedPair(aln))
                else:
                    yield ('PX', AlignedPair(aln))

    for aln in readcache.values():
        yield ('cached', AlignedPair(aln))

    #
    # print(len(readcache))


    #             if aln.mate_is_unmapped and aln.is_secondary:
    #                 ''' Mate will not be present since it did not align.
    #                     Unmapped mates are only present for the primary
    #                     alignment.'''
    #                 yield ('PX', (aln, None))
    #             else:
    #
    #
    #
    #
    #         if mate:
    #             ''' Mate was found in cache'''
    #             if aln.is_read1:
    #                 yield ('PM', (aln, mate))
    #             else:
    #                 yield ('P', (mate, aln))
    #         else:
    #             ''' Mate was not found '''
    #             if mate_before(aln):
    #                 ''' Mate is before. Do not cache current alignment since
    #                     mate will not be visited with the current iterator.'''
    #                 yield ('MB', (aln, None))
    #             elif mate_in_region(aln, regtup):
    #                 ''' Mate is within region. Cache current alignment since
    #                     mate might still be visited with current iterator.'''
    #                 readcache[readkeyP(aln)] = aln
    #             else:
    #                 ''' Mate is outside region. Cache current alignment in
    #                     outcache and handle later.'''
    #                 outcache[readkeyP(aln)] = aln
    #
    # for aln in readcache.values():
    #     # assert aln.get_tag('YT') == 'UP', '{}'.format(aln)
    #     yield ('readcache', (aln, None))
    # for aln in outcache.values():
    #     yield ('outcache', (aln, None))

"""

    print('{} reads in outcache'.format(len(outcache)))
    print('cache contains:')

    yt_tags = Counter()
    for key,aln in readcache.items():
        if aln.has_tag('YT'):
            yt_tags[aln.get_tag('YT')] += 1
        else:
            yt_tags['no_YT'] += 1
    for t,c in yt_tags.most_common():
        print('{}: {}'.format(t,c))

    is_unmapped = sum(a.is_unmapped for a in readcache.values())
    print('{} unmapped'.format(is_unmapped))
    mate_unmapped = sum(a.mate_is_unmapped for a in readcache.values())
    print('{} mate unmapped'.format(mate_unmapped))
    # vis = sum(readkeyP(a) in visited for a in readcache.values())
    # print('{} visited mate'.format(vis))
    bef = sum(mate_before(a) for a in readcache.values())
    print('{} mate before'.format(bef))
    if PRINTCACHE:
        print('**outcache**')
        x = list(outcache.values())
        for _ in x[:10]:
            print(_)
            print(readkeyP(_))
            print(matekeyP(_))
            print()

    print('outcache contains:')
    yt_tags = Counter()
    for key, aln in outcache.items():
        if aln.has_tag('YT'):
            yt_tags[aln.get_tag('YT')] += 1
        else:
            yt_tags['no_YT'] += 1
    for t, c in yt_tags.most_common():
         print('{}: {}'.format(t, c))
"""


from . import model
import os

def fetch_region(samfile, annotation, opts, region):
    print('processing {}:{}-{}'.format(*region))

    _nfkey = opts['no_feature_key']
    _omode, _othresh = opts['overlap_mode'], opts['overlap_threshold']
    _tempdir = opts['tempdir']

    assign = model.Assigner(annotation, _nfkey, _omode, _othresh).assign_func()

    alninfo = Counter()
    alninfo['minAS'] = 999999999
    alninfo['maxAS'] = -999999999

    mfile = os.path.join(_tempdir, 'tmp_map.{}.{}.{}.txt'.format(*region))
    mappings = []

    fh = open(mfile, 'w')
    with pysam.AlignmentFile(samfile) as sf:
        samiter = sf.fetch(*region)
        regtup = (sf.get_tid(region[0]), region[1], region[2])
        for code, aln in fetch_pairs_sorted(samiter, regtup):
            alninfo[code] += 1
            if code == 'PU' or code == 'SU':
                continue

            m = (code, aln.query_id, assign(aln), aln.alnscore, aln.alnlen)
            alninfo['minAS'] = min(alninfo['minAS'], m[3])
            alninfo['maxAS'] = max(alninfo['maxAS'], m[3])
            print('\t'.join(map(str, m)), file=fh)
            mappings.append(m)

    fh.close()
    return alninfo, mfile, mappings
