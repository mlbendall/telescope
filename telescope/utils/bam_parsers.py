# -*- coding: utf-8 -*-
""" Parse SAM/BAM alignment files

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import logging
import pysam
from collections import namedtuple, Counter

from .helpers import GenomeRegion
from .helpers import merge_blocks

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

AlignmentData = namedtuple('AlignmentData', ('readid','featid','ascore','alen'))

def readkey(aln):
    ''' Key for read '''
    return (aln.query_name, aln.is_read1,
            aln.reference_id, aln.reference_start,
            aln.next_reference_id, aln.next_reference_start)


def matekey(aln):
    ''' Key for mate '''
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


def fetch_pairs(samfile, reference=None, start=None, end=None, region=None):
    _region = GenomeRegion(reference, start, end, region)
    readcache = {}
    outregion = []
    samiter = samfile.fetch(_region.chrom, _region.start, _region.end)
    for aln in samiter:
        if not aln.is_paired:
            yield (aln, )
        else:
            mate = readcache.pop(matekey(aln), None)
            if mate is not None:                               # Mate found in cache
                yield (aln, mate) if aln.is_read1 else (mate, aln)
            else:                                              # Mate not found in cache
                if _region.contains(aln.next_reference_name, aln.next_reference_start):
                    if mate_before(aln):                       # Mate was expected before this
                        yield (aln, )
                    else:                                      # Mate is still to come
                        readcache[readkey(aln)] = aln
                else:                                          # Mate is outside region
                    outregion.append(aln)
    # Resolve reads whose mate is outside the region
    for aln in outregion:
        mate = None
        for m in samfile.fetch(_region.chrom, aln.next_reference_start, aln.next_reference_start+1):
            if readkey(aln)==matekey(m) and matekey(aln)==readkey(m):
                mate = m
                break
        if mate is None:
            yield (aln, )
        else:
            yield (aln, mate) if aln.is_read1 else (mate, aln)
    # Yield the remaining reads in the cache as unpaired
    for aln in readcache.values():
        yield (aln, )


def iter_pairs(samfile):
    readcache = {}
    samiter = samfile.fetch(until_eof=True)
    for aln in samiter:
        if not aln.is_paired:
            yield (aln,)
        else:
            mate = readcache.pop(matekey(aln), None)
            if mate is not None:                               # Mate found in cache
                yield (aln, mate) if aln.is_read1 else (mate, aln)
            else:                                              # Mate not found in cache
                readcache[readkey(aln)] = aln
    # Yield the remaining reads in the cache as unpaired
    for aln in readcache.values():
        yield (aln, )

""" New functions """

def refblocks(aligned_pair):
    if len(aligned_pair) == 1:
        blocks = aligned_pair[0].get_blocks()
    else:
        blocks = aligned_pair[0].get_blocks() + aligned_pair[1].get_blocks()
    return merge_blocks(blocks, 1)

def alnlen(aligned_pair):
    return sum(b[1] - b[0] for b in refblocks(aligned_pair))

def fragid(aligned_pair):
    if len(aligned_pair) == 1:
        if aligned_pair[0].is_read1:
            return '{}/1'.format(aligned_pair[0].query_name)
        else:
            return '{}/2'.format(aligned_pair[0].query_name)
    else:
        return aligned_pair[0].query_name

def alnscore(aligned_pair):
    if len(aligned_pair) == 1:
        return aligned_pair[0].get_tag('AS')
    else:
        return aligned_pair[0].get_tag('AS') + aligned_pair[1].get_tag('AS')


def pair_alignments(alniter):
    readcache = {}
    for aln in alniter:
        if not aln.is_paired:
            yield (aln,)
        else:
            mate = readcache.pop(matekey(aln), None)
            if mate is not None:  # Mate found in cache
                yield (aln, mate) if aln.is_read1 else (mate, aln)
            else:  # Mate not found in cache
                readcache[readkey(aln)] = aln
    # Yield the remaining reads in the cache as unpaired
    for aln in readcache.values():
        yield (aln,)

def organize_bundle(alns):
    if not alns[0].is_paired:
        assert all(not a.is_paired for a in alns), 'mismatch pair flag'
        return [(a,) for a in alns], None
    else:
        if len(alns) == 2 and alns[0].is_unmapped and alns[1].is_unmapped:
            # Unmapped fragment
            return [(alns[0], alns[1])], None
        elif alns[0].is_proper_pair:
            return list(pair_alignments(alns)), None
        else:
            ret1, ret2 = list(), list()
            for aln in alns:
                if aln.is_read1:
                    ret1 += [(aln,)]
                else:
                    assert aln.is_read2
                    ret2 += [(aln,)]
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


def load_unsorted(sam_fn, annot, opts, alninfo):
    assigner = Assigner(annot, opts)
    _nfkey = opts.no_feature_key
    with pysam.AlignmentFile(sam_fn) as sf:
        for pairs in fetch_fragments(sf, until_eof=True):
            if pairs[0][0].is_unmapped:
                alninfo['unmap_{}'.format(len(pairs[0]))] += 1
                continue
            _ambig = len(pairs) > 1
            alninfo['map_{}'.format(len(pairs[0]))] += 1
            overlap_feats = list(map(assigner.assign_pair, pairs))
            has_overlap = any(f != _nfkey for f in overlap_feats)
            if not has_overlap:
                alninfo['nofeat_{}'.format('A' if _ambig else 'U')] += 1
                continue

            alninfo['feat_{}'.format('A' if _ambig else 'U')] += 1
            d = {}
            for pair, feat in zip(pairs, overlap_feats):
                _alnscore = alnscore(pair)
                alninfo['minAS'] = min(alninfo['minAS'], _alnscore)
                alninfo['maxAS'] = max(alninfo['maxAS'], _alnscore)
                yield fragid(pair), feat, _alnscore, alnlen(pair)

class Assigner:
    def __init__(self, annotation, opts):
        self.annotation = annotation
        self.opts = opts

    def _assign_pair_threshold(self, pair):
        assert not pair[0].is_unmapped, 'ERROR: only mapped reads are allowed'
        blocks = refblocks(pair)
        f = self.annotation.intersect_blocks(pair[0].reference_name, blocks)
        if not f:
            return self.opts.no_feature_key
        # Calculate the percentage of fragment mapped
        alnlen = sum(b[1]-b[0] for b in blocks)
        fname, overlap = f.most_common()[0]
        if overlap > alnlen * self.opts.overlap_threshold:
            return fname
        else:
            return self.opts.no_feature_key

    def assign_pair(self, pair):
        if self.opts.overlap_mode == 'threshold':
            return self._assign_pair_threshold(pair)
        else:
            assert False
