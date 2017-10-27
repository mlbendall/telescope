# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

import pysam
import HTSeq

from .helpers import merge_blocks

def readkey(aln):
    ''' This is a key to identify read '''
    return (aln.query_name, aln.is_read1,
            aln.reference_id, aln.reference_start,
            aln.next_reference_id, aln.next_reference_start)

def matekey(aln):
    ''' The expected key of the mate for this read '''
    return (aln.query_name, not aln.is_read1,
            aln.next_reference_id, aln.next_reference_start,
            aln.reference_id, aln.reference_start)

def pair_alignments(alns):
    readcache = {}
    for aln in alns:
        if aln.is_proper_pair:
            mate = readcache.pop(matekey(aln), None)
            if mate is not None:
                yield (aln, mate) if aln.is_read1 else (mate, aln)
            else:
                readcache[readkey(aln)] = aln
        else:
            yield (aln, None)
    assert len(readcache) == 0


def alnblocks(aligned_pair):
    if aligned_pair[1] is None:
        blocks = aligned_pair[0].get_blocks()
    else:
        blocks = aligned_pair[0].get_blocks() + aligned_pair[1].get_blocks()
    return merge_blocks(blocks)

def alnscore(aligned_pair):
    if aligned_pair[1] is None:
        return aligned_pair[0].get_tag('AS')
    else:
        return aligned_pair[0].get_tag('AS') + aligned_pair[1].get_tag('AS')

def fragid(aligned_pair):
    if aligned_pair[1] is None:
        if aligned_pair[0].is_read1:
            return '{}/1'.format(aligned_pair[0].query_name)
        else:
            return '{}/2'.format(aligned_pair[0].query_name)
    else:
        return aligned_pair[0].query_name

def alnlen(aligned_pair):
    return sum(b[1] - b[0] for b in alnblocks(aligned_pair))



"""
def pair_alignments(alns):
    readcache = {}
    for aln in alns:
        if not aln.is_paired:
            yield (aln, None)
        else:
            alnmode = aln.get_tag("YT")
            if alnmode == 'UP':
                yield (aln, None)
            elif aln.reference_id != aln.next_reference_id:
                yield (aln, None)
            else:
                mate = readcache.pop(matekey(aln), None)
                if mate is not None:
                    yield (aln, mate) if aln.is_read1 else (mate, aln)
                else:
                    readcache[readkey(aln)] = aln
    if len(readcache) != 0:
        print('here')
"""


# def _process_bundle(alns):
#     pairs = list(pair_alignments(alns))
#     if
#     for aln in alns:
#         if aln.is_pro
#         pairs = list(pair_alignments(b))
#
#
# def fetch_fragment_alignments(samfile, **kwargs):
#
#     samiter = samfile.fetch(**kwargs)
#     bundle = [ next(samiter) ]
#     for aln in samiter:
#         if aln.query_name == bundle[0].query_name:
#             bundle.append(aln)
#         else:
#
#             bundle = [ aln ]
#     yield pair_alignments(bundle)

def fetch_bundle_pairs(samfile, **kwargs):
    """ Iterate over alignment over reads with same ID """
    samiter = samfile.fetch(**kwargs)
    bundle = [ next(samiter) ]
    for aln in samiter:
        if aln.query_name == bundle[0].query_name:
            bundle.append(aln)
        else:
            yield pair_alignments(bundle)
            bundle = [ aln ]
    yield pair_alignments(bundle)


def fetch_bundle(samfile, **kwargs):
    """ Iterate over alignment over reads with same ID """
    samiter = samfile.fetch(**kwargs)
    bundle = [ next(samiter) ]
    for aln in samiter:
        if aln.query_name == bundle[0].query_name:
            bundle.append(aln)
        else:
            yield bundle
            bundle = [ aln ]
    yield bundle



def get_overlap_function(annotation, overlap_mode, overlap_threshold,
                            no_feature_key):

    def _assign_pair_threshold(pair):
        assert not pair[0].is_unmapped, 'ERROR: only mapped reads are allowed'
        blocks = get_blocks_pair(pair)
        f = annotation.intersect_blocks(pair[0].reference_name, blocks)
        if not f:
            return no_feature_key
        # Calculate the percentage of fragment mapped
        alnlen = sum(b[1]-b[0] for b in blocks)
        fname, overlap = f.most_common()[0]
        if overlap > alnlen * overlap_threshold:
            return fname
        else:
            return no_feature_key

    if overlap_mode == 'threshold':
        return _assign_pair_threshold


"""

class Assigner(object):
    def __init__(self, annotation, overlap_mode, overlap_threshold,
                 no_feature_key, not_aligned_key="__not_aligned"):
        self.annotation = annotation
        self.overlap_mode = overlap_mode
        self.overlap_threshold = overlap_threshold
        self.no_feature_key = no_feature_key
        self.not_aligned_key = not_aligned_key
        self._skipfeat = set([no_feature_key, not_aligned_key])

    def _assign_pair(self, pair):
        if pair[0].is_unmapped:
            return self.not_aligned_key
        else:
            blocks = get_blocks_pair(pair)
            f = self.annotation.intersect_blocks(pair[0].reference_name, blocks)
            if not f:
                return self.no_feature_key

            alnlen = sum(b[1] - b[0] for b in blocks)
            # Use a threshold to determine overlapping features. Fragment is
            # counted for this features if the proportion of positions contained
            # within the feature exceeds threshold. If multiple features meet this
            # criteria, feature with largest proportion is selected.
            if self.overlap_mode == 'threshold':
                fname, overlap = f.most_common()[0]
                if overlap > alnlen * self.overlap_threshold:
                    return fname

            # Strict intersection of all overlapping features, all positions in
            # alignment must be contained within feature. Ambiguous if this criteria
            # applies to more than one feature.
            if self.overlap_mode == 'intersection-strict':
                pass

            # Union of all overlapping features. Ambiguous if alignment overlaps
            # more than one feature.
            if self.overlap_mode == 'union':
                if len(feats) > 1:
                    logging.debug("Overlap is ambiguous: %s" % str(feats))
                return feats[0]

    def assign_bundle(self, bundle):
        pairs = list(bundle)
        # Number of fragments: len(pairs),
        # Number of alignments: sum(1 if p[1] is None else 2 for p in pairs)
        # Number unmapped:
        features = list(map(self._assign_pair, pairs))
        if all(f in self._skipfeat for f in features):
            return None
        alnscores = list(map(alnscore, pairs))
        fragids = list(map(fragid, pairs))
        return zip(features, alnscores, fragids)
"""

def bundle_stats(pairs):
    return(
    len(pairs), # Number of fragments
    sum(1 if p[1] is None else 2 for p in pairs), # Number of alignments
    sum(p[0].is_unmapped for p in pairs), # Number unmapped
    )

"""

                pairs = list(bundle)
                ''' Determine if any alignments overlap annotation '''
                feats = []
                for pair in pairs:
                    r = assigner._assign_pair(pair)
                    print(r)

                    alninfo['num_aln'] += 2 if pair[1] else 1
                    if pair[0].is_unmapped:
                        alninfo['num_unmap'] += 1
                        assert pair[1] is None
                        feats.append(Counter())
                    else:
                        feats.append(self.annotation.intersect_alignment(pair))

                ''' Finished if no overlaps found '''
                if sum(len(f) for f in feats) == 0:
                    alninfo['num_nonherv'] += len(pairs)
                    continue

                ''' Get the alignment scores '''
                scores = [alignment.alnscore(p) for p in pairs]
                if pairs[0][1] is None:
                    pass

        logging.info('{} fragments'.format(alninfo['num_frag']))
        logging.info('{} alignments'.format(alninfo['num_aln']))
        logging.info('{} unmapped'.format(alninfo['num_unmap']))
        logging.info('{} mapped'.format(alninfo['num_map']))

    def resolve_overlap(self, feats, alen):
        # No overlapping features found
        if len(feats) == 0:
            return self.nofeat_str, alen

        # Use a threshold to determine overlapping features. Fragment is
        # counted for this features if the proportion of positions contained
        # within the feature exceeds threshold. If multiple features meet this
        # criteria, feature with largest proportion is selected.
        if self.overlap_mode == 'threshold':
            if float(feats[0][1]) / alen >= self.threshold:
                return feats[0]
            else:
                logging.debug('Overlap fails to meet threshold: (len=%d, %s)'
                              % (alen, feats))
                return self.nofeat_str, alen
            #
            # ok_feats = [(f, ol) for f, ol in feats
            #             if float(ol) / alen >= self.threshold]
            # if not ok_feats:
            #     logging.debug('Discarded (len=%d, %s)' % (alen, feats))
            #     return (self.nofeat_str, alen)
            # else:
            #     # Check if value is repeated
            #     if sum(ol==feats[0][1] for f, ol in feats) > 1:
            #         logging.warning('Multiple with max OL (len=%d, %s)' % (alen, feats))
            #     return ok_feats[0]

        # Strict intersection of all overlapping features, all positions in
        # alignment must be contained within feature. Ambiguous if this criteria
        # applies to more than one feature.
        if self.overlap_mode == 'intersection-strict':
            if feats[0][1] >= alen:
                return feats[0]
            else:
                logging.debug('Overlap fails strict intersection: (len=%d, %s)'
                              % (alen, feats))
                return self.nofeat_str, alen

        # Union of all overlapping features. Ambiguous if alignment overlaps
        # more than one feature.
        if self.overlap_mode == 'union':
            if len(feats) > 1:
                logging.debug("Overlap is ambiguous: %s" % str(feats))
            return feats[0]
"""



"""
def load_unsorted(samfile_path, annot, opts):

    # assigner = Assigner(annot, opts)
    ret = []
    with pysam.AlignmentFile(samfile_path) as sf:
        for bundle in fetch_bundle_pairs(sf):
            alninfo['num_frag'] += 1
            frags = list(bundle)
            for pair in frags:
                alninfo['num_aln'] += len(pair)
                for a in pair:
                    if a.is_unmapped:
                        alninfo['num_unmap'] += 1
                    else:
                        alninfo['num_map'] += 1
    return alninfo
"""




"""
        for pair in iter_pairs(sf):
            for a in pair:
                if a.is_unmapped:
                    alninfo['num_unmap'] += 1
                else:
                    alninfo['num_map'] += 1
            for frag in assigner.assign_alignment(pair):
                alninfo['minAS'] = min(frag[0].ascore, alninfo['minAS'])
                alninfo['maxAS'] = max(frag[0].ascore, alninfo['maxAS'])
                if len(frag[1]) == 2:
                    alninfo['aligned_pairs'] += 1
                elif len(frag[1]) == 1:
                    alninfo['aligned_singletons'] += 1
                ret.append(frag)
            # ret.extend(assigner.assign_alignment(pair))
    return ret, alninfo
"""



"""
class Assigner(object):
    def __init__(self, annotation, opts):
        self.annotation = annotation
        self.nofeat_str = opts.no_feature_key
        self.overlap_mode = opts.overlap_mode
        self.threshold = opts.overlap_threshold

    def resolve_overlap(self, feats, alen):
        # No overlapping features found
        if len(feats) == 0:
            return self.nofeat_str, alen

        # Use a threshold to determine overlapping features. Fragment is
        # counted for this features if the proportion of positions contained
        # within the feature exceeds threshold. If multiple features meet this
        # criteria, feature with largest proportion is selected.
        if self.overlap_mode == 'threshold':
            if float(feats[0][1]) / alen >= self.threshold:
                return feats[0]
            else:
                logging.debug('Overlap fails to meet threshold: (len=%d, %s)'
                              % (alen, feats))
                return self.nofeat_str, alen
            #
            # ok_feats = [(f, ol) for f, ol in feats
            #             if float(ol) / alen >= self.threshold]
            # if not ok_feats:
            #     logging.debug('Discarded (len=%d, %s)' % (alen, feats))
            #     return (self.nofeat_str, alen)
            # else:
            #     # Check if value is repeated
            #     if sum(ol==feats[0][1] for f, ol in feats) > 1:
            #         logging.warning('Multiple with max OL (len=%d, %s)' % (alen, feats))
            #     return ok_feats[0]

        # Strict intersection of all overlapping features, all positions in
        # alignment must be contained within feature. Ambiguous if this criteria
        # applies to more than one feature.
        if self.overlap_mode == 'intersection-strict':
            if feats[0][1] >= alen:
                return feats[0]
            else:
                logging.debug('Overlap fails strict intersection: (len=%d, %s)'
                              % (alen, feats))
                return self.nofeat_str, alen

        # Union of all overlapping features. Ambiguous if alignment overlaps
        # more than one feature.
        if self.overlap_mode == 'union':
            if len(feats) > 1:
                logging.debug("Overlap is ambiguous: %s" % str(feats))
            return feats[0]

    def _assign_single(self, aln):
        rkey = (readkey(aln),)
        # Get read name and add mate indicator if paired
        query_name = aln.query_name
        if aln.is_paired:
            query_name += '\x97%d' % (1 if aln.is_read1 else 2)
        ref_id = aln.reference_id

        # Get the alignment score
        aln_score = aln.get_tag('AS')

        # List of aligned gapless blocks
        blocks = merge_interval_list(aln.get_blocks(), 1)
        aln_len = sum(bend - bstart for bstart, bend in blocks)

        # Map to features
        feats = map_to_features(ref_id, blocks, self.annotation)
        if not feats:
            return AlignmentData(query_name, self.nofeat_str, aln_score,
                                 aln_len), rkey
        else:
            best_feat, best_olen = self.resolve_overlap(feats, aln_len)
            return AlignmentData(query_name, best_feat, aln_score,
                                 aln_len), rkey

    def _assign_pair(self, aln1, aln2):
        rkey = (readkey(aln1), readkey(aln2))
        query_name = aln1.query_name
        ref_id = aln1.reference_id

        # Calculate total alignment score
        aln_score = aln1.get_tag('AS') + aln2.get_tag('AS')

        # List of aligned gapless blocks
        blocks = merge_interval_list(aln1.get_blocks() + aln2.get_blocks(), 1)
        aln_len = sum(bend - bstart + 1 for bstart, bend in blocks)

        # Map to features
        feats = map_to_features(ref_id, blocks, self.annotation)
        if not feats:
            return AlignmentData(query_name, self.nofeat_str, aln_score,
                                 aln_len), rkey
        else:
            best_feat, best_olen = self.resolve_overlap(feats, aln_len)
            return AlignmentData(query_name, best_feat, aln_score,
                                 aln_len), rkey

    def assign_alignment(self, alns):
        if len(alns) == 1:
            return [self._assign_single(alns[0]), ]
        a1, a2 = alns
        if a1.is_unmapped and a2.is_unmapped:
            return []
        if a1.is_unmapped or a2.is_unmapped:
            # One of the mates is unmapped
            if a1.is_unmapped:
                return [self._assign_single(a2), ]
            else:
                return [self._assign_single(a1), ]
        # Both mates are mapped
        if a1.is_proper_pair:
            return [self._assign_pair(a1, a2), ]
        else:
            return [self._assign_single(a1), self._assign_single(a2)]

"""



"""

        # alninfo = alignment.load_unsorted(self.opts.samfile, self.annotation, self.opts)

        alninfo = {'num_map': 0,
                   'num_unmap': 0,
                   'minAS': float('inf'),
                   'maxAS': float('-inf'),
                   'aligned_pairs': 0,
                   'aligned_singletons': 0,
                   }
        assigner = alignment.Assigner(annot, opts)
        ret = []
        with pysam.AlignmentFile(self.opts.samfile) as sf:
            for bundle in alignment.fetch_bundle_pairs(sf):
                aln_info['num_frag'] += 1
                pairs = list(bundle)
                if pairs[0][0].is_unmapped:
                    assert pairs[0][1] is None


                if pairs[0][1] is None:
                    if pairs[0][0].is_unmapped:
                        aln_info['num_unmap'] += 1
                if pairs[0][1] is None:
                    assert all(t[1] is None for t in pairs)
                    aln_info['num_aln'] += len(pairs)
                else:
                    assert all(t[1] is not None for t in pairs)
                    aln_info['num_aln'] += (len(pairs) * 2)


                # aligned_locs = {}

                for r1, r2 in pairs:
                    nbundles[r1.query_name] += 1
                    if r2 is not None:
                        nbundles[r1.query_name] += 1
                    f  = self.annotation.intersect(r1, r2)
                    if f:
                        killer += 1
                        print(f)
                    if killer > 10:
                        sys.exit()


import pysam
samfile = pysam.AlignmentFile('/Users/bendall/Development/telescope_test/data/aligned.100K.sam')
alnmodes = defaultdict(list)
for bundle in fetch_bundle(samfile):
    pp = bundle[0].is_proper_pair
    alnmode = bundle[0].get_tag("YT")
    assert all(alnmode==a.get_tag("YT") for a in bundle)
    assert all(pp == a.is_proper_pair for a in bundle)
    for aln in bundle:
        alnmodes[alnmode].append(aln)

sum(a.is_proper_pair for a in alnmodes['CP'])
[0].is_proper_pair


    for aln in bundle:



        if aln.get_tag('YT') == 'DP':
            alnmodes[]


        if aln.is_proper_pair:
            alnmodes.add(aln.get_tag("YT"))

        alnmodes.add(aln.get_tag("YT"))
    if len(alnmodes) != 1:
        print('found')
        break
"""

class Assigner:
    def __init__(self, annotation, opts):
        self.annotation = annotation
        self.nofeat_str = opts.no_feature_key
        self.overlap_mode = opts.overlap_mode
        self.threshold = opts.overlap_threshold

    def resolve_overlap(self, feats, alen):
        # No overlapping features found
        if len(feats) == 0:
            return self.nofeat_str, alen

        # Use a threshold to determine overlapping features. Fragment is
        # counted for this features if the proportion of positions contained
        # within the feature exceeds threshold. If multiple features meet this
        # criteria, feature with largest proportion is selected.
        if self.overlap_mode == 'threshold':
            if float(feats[0][1]) / alen >= self.threshold:
                return feats[0]
            else:
                logging.debug('Overlap fails to meet threshold: (len=%d, %s)'
                              % (alen, feats))
                return self.nofeat_str, alen
            #
            # ok_feats = [(f, ol) for f, ol in feats
            #             if float(ol) / alen >= self.threshold]
            # if not ok_feats:
            #     logging.debug('Discarded (len=%d, %s)' % (alen, feats))
            #     return (self.nofeat_str, alen)
            # else:
            #     # Check if value is repeated
            #     if sum(ol==feats[0][1] for f, ol in feats) > 1:
            #         logging.warning('Multiple with max OL (len=%d, %s)' % (alen, feats))
            #     return ok_feats[0]

        # Strict intersection of all overlapping features, all positions in
        # alignment must be contained within feature. Ambiguous if this criteria
        # applies to more than one feature.
        if self.overlap_mode == 'intersection-strict':
            if feats[0][1] >= alen:
                return feats[0]
            else:
                logging.debug('Overlap fails strict intersection: (len=%d, %s)'
                              % (alen, feats))
                return self.nofeat_str, alen

        # Union of all overlapping features. Ambiguous if alignment overlaps
        # more than one feature.
        if self.overlap_mode == 'union':
            if len(feats) > 1:
                logging.debug("Overlap is ambiguous: %s" % str(feats))
            return feats[0]

    def _assign_single(self, aln):
        rkey = (readkey(aln),)
        # Get read name and add mate indicator if paired
        query_name = aln.query_name
        if aln.is_paired:
            query_name += '\x97%d' % (1 if aln.is_read1 else 2)
        ref_id = aln.reference_id

        # Get the alignment score
        aln_score = aln.get_tag('AS')

        # List of aligned gapless blocks
        blocks = merge_interval_list(aln.get_blocks(), 1)
        aln_len = sum(bend - bstart for bstart, bend in blocks)

        # Map to features
        feats = map_to_features(ref_id, blocks, self.annotation)
        if not feats:
            return AlignmentData(query_name, self.nofeat_str, aln_score,
                                 aln_len), rkey
        else:
            best_feat, best_olen = self.resolve_overlap(feats, aln_len)
            return AlignmentData(query_name, best_feat, aln_score,
                                 aln_len), rkey

    def _assign_pair(self, aln1, aln2):
        rkey = (readkey(aln1), readkey(aln2))
        query_name = aln1.query_name
        ref_id = aln1.reference_id

        # Calculate total alignment score
        aln_score = aln1.get_tag('AS') + aln2.get_tag('AS')

        # List of aligned gapless blocks
        blocks = merge_interval_list(aln1.get_blocks() + aln2.get_blocks(), 1)
        aln_len = sum(bend - bstart + 1 for bstart, bend in blocks)

        # Map to features
        feats = map_to_features(ref_id, blocks, self.annotation)
        if not feats:
            return AlignmentData(query_name, self.nofeat_str, aln_score,
                                 aln_len), rkey
        else:
            best_feat, best_olen = self.resolve_overlap(feats, aln_len)
            return AlignmentData(query_name, best_feat, aln_score,
                                 aln_len), rkey

    def assign_alignment(self, alns):
        if len(alns) == 1:
            return [self._assign_single(alns[0]), ]
        a1, a2 = alns
        if a1.is_unmapped and a2.is_unmapped:
            return []
        if a1.is_unmapped or a2.is_unmapped:
            # One of the mates is unmapped
            if a1.is_unmapped:
                return [self._assign_single(a2), ]
            else:
                return [self._assign_single(a1), ]
        # Both mates are mapped
        if a1.is_proper_pair:
            return [self._assign_pair(a1, a2), ]
        else:
            return [self._assign_single(a1), self._assign_single(a2)]

def load_unsorted(samfile_path, annot, opts):
    alninfo = {'num_map': 0,
               'num_unmap': 0,
               'minAS': float('inf'),
               'maxAS': float('-inf'),
               'aligned_pairs': 0,
               'aligned_singletons': 0,
               }
    assigner = Assigner(annot, opts)
    ret = []
    with pysam.AlignmentFile(samfile_path) as sf:
        for pair in iter_pairs(sf):
            for a in pair:
                if a.is_unmapped:
                    alninfo['num_unmap'] += 1
                else:
                    alninfo['num_map'] += 1
            for frag in assigner.assign_alignment(pair):
                alninfo['minAS'] = min(frag[0].ascore, alninfo['minAS'])
                alninfo['maxAS'] = max(frag[0].ascore, alninfo['maxAS'])
                if len(frag[1]) == 2:
                    alninfo['aligned_pairs'] += 1
                elif len(frag[1]) == 1:
                    alninfo['aligned_singletons'] += 1
                ret.append(frag)
            # ret.extend(assigner.assign_alignment(pair))
    return ret, alninfo
