# -*- coding: utf-8 -*-
""" Parse SAM/BAM alignment files

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

from telescope.utils.calignment import AlignedPair
# from ._alignment import AlignedPair


def readkey(aln):
    ''' Key that is unique for alignment '''
    return (aln.query_name, aln.is_read1,
            aln.reference_id, aln.reference_start,
            aln.next_reference_id, aln.next_reference_start)


def matekey(aln):
    ''' Expected key of mate '''
    return (aln.query_name, not aln.is_read1,
            aln.next_reference_id, aln.next_reference_start,
            aln.reference_id, aln.reference_start)


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


def pair_alignments(alniter):
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


def organize_bundle(alns):
    if not alns[0].is_paired:
        assert all(not a.is_paired for a in alns), 'mismatch pair flag'
        return [AlignedPair(a) for a in alns], None
    else:
        if len(alns) == 2 and alns[0].is_unmapped and alns[1].is_unmapped:
            # Unmapped fragment
            return [AlignedPair(alns[0], alns[1])], None
        elif alns[0].is_proper_pair:
            return list(pair_alignments(alns)), None
        else:
            ret1, ret2 = list(), list()
            for aln in alns:
                if aln.is_read1:
                    ret1 += [AlignedPair(aln)]
                else:
                    assert aln.is_read2
                    ret2 += [AlignedPair(aln)]
            return ret1, ret2


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


def fetch_fragments(samfile, **kwargs):
    biter = fetch_bundle(samfile, **kwargs)
    for bundle in biter:
        f1, f2 = organize_bundle(bundle)
        yield f1
        if f2 is not None:
            yield f2
