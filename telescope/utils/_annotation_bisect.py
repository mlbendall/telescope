# -*- coding: utf-8 -*-

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

from __future__ import division
from builtins import range
from past.utils import old_div
from builtins import object
__author__ = 'bendall'

import re
from collections import defaultdict, namedtuple, Counter, OrderedDict
from bisect import bisect_left,bisect_right

GTFRow = namedtuple('GTFRow', ['chrom','source','feature','start','end','score','strand','frame','attribute'])
# BEDRow = namedtuple('BEDRow', ['chrom','start','end','name','score','strand','strand','frame','attribute'])

# class _AnnotationLinear(object):
#     def __init__(self, gtf_file, attribute_name):
#         self.key = attribute_name
#         self.loci = OrderedDict()
#         self.intervals = {}
#
#         fh = open(gtf_file,'rU') if isinstance(gtf_file,str) else gtf_file
#         features = (GTFRow(*l.strip('\n').split('\t')) for l in fh if not l.startswith('#'))
#         for f in features:
#             attr = dict(re.findall('(\w+)\s+"(.+?)";', f.attribute))
#             ''' Add to locus list '''
#             if attr[self.key] not in self.loci:
#                 self.loci[attr[self.key]] = list()
#             self.loci[attr[self.key]].append(f)
#
#     def _addlocus(self, chrom, spos, epos, name):
#         if chrom not in self.intervals:
#             self.intervals[chrom] = [
#                 [float('-inf'), spos, set(), ],
#                 [spos,          epos, set([name])],
#                 [epos, float('inf'), set()],
#             ]
#
#     def _lookup(self, chrom, pos):
#
#
#
#

class _AnnotationBisect(object):
    def __init__(self, gtffile,  min_overlap=None, attr_name="locus"):
        self.key = attr_name

        # Instance variables
        self._locus = []                      # List of locus names
        self._locus_lookup = defaultdict(list) # {}
        self._intervals = defaultdict(list)   # Dictionary containing lists of intervals for each reference
        self._intS = {}                       # Dictionary containing lists of interval start positions for each reference
        self._intE = {}                       # Dictionary containing lists of interval end positions for each reference

        # GTF filehandle
        fh = open(gtffile,'rU') if isinstance(gtffile,str) else gtffile
        features = (GTFRow(*l.strip('\n').split('\t')) for l in fh if not l.startswith('#'))
        for i,f in enumerate(features):
            attr = dict(re.findall('(\w+)\s+"(.+?)";', f.attribute))
            _locus_name = attr[self.key] if self.key in attr else 'TELE%04d' % i
            if _locus_name not in self._locus:
                self._locus.append(_locus_name)
            # else:
            #     assert False, "Non-unique locus name found: %s" % _locus_name
            #self._locus.append( attr[attr_name] if attr_name in attr else 'PSRE%04d' % i )
            self._locus_lookup[_locus_name].append( (f.chrom, int(f.start), int(f.end)) )
            # self._intervals[f.chrom].append((int(f.start), int(f.end), i))
            self._intervals[f.chrom].append((int(f.start), int(f.end), _locus_name))

        # Sort intervals by start position
        for chrom in list(self._intervals.keys()):
            self._intervals[chrom].sort(key=lambda x:x[0])
            self._intS[chrom] = [s for s,e,i in self._intervals[chrom]]
            self._intE[chrom] = [e for s,e,i in self._intervals[chrom]]

    def lookup(self, chrom, pos):
        ''' Return the feature for a given reference and position '''
        if chrom not in self._intervals:
            return None

        sidx = bisect_right(self._intS[chrom], pos)   # Return index of leftmost interval where start > pos
        # If the end position is inclusive (as in GTF) use bisect_left
        eidx = bisect_left(self._intE[chrom], pos)   # Return index of leftmost interval where end >= pos

        # If sidx == eidx, the position is between intervals at (sidx-1) and (sidx)
        # If eidx < sidx, the position is within eidx
        feats = [self._intervals[chrom][i] for i in range(eidx,sidx)]
        if len(feats) == 0:
            return None
        else:
            possible = set(f[2] for f in feats)
            assert len(possible) == 1, '%s' % feats
            return possible.pop()

    def lookup_interval(self, chrom, spos, epos):
        ''' Resolve the feature that overlaps or contains the given interval
            NOTE: Only tests the start and end positions. This means that it does not handle
                  cases where a feature lies completely within the interval. This is OK when the
                  fragment length is expected to be smaller than the feature length.

                  Fragments where ends map to different features are resolved by
                  assigning the larger of the two overlaps.
        '''
        featL = self.lookup(chrom, spos)
        featR = self.lookup(chrom, epos)
        if featL is None and featR is None:     # Neither start nor end is within a feature
            return None
        else:
            if featL is None or featR is None:    # One (and only one) end is within a feature
                return featL if featL is not None else featR
            elif featL == featR:                  # Both ends within the same feature
                return featL
            else:                                 # Ends in different features
                locL = self._locus_lookup[featL][-1]   # Assume last fragment
                locR = self._locus_lookup[featR][0]    # Assume first fragment
                overlapL = locL[2] - spos
                overlapR = epos - locR[1]
                if overlapL >= overlapR:
                    return featL
                else:
                    return featR

    def feature_length(self):
        _ret = {}
        for chr,ilist in self._intervals.items():
            for spos,epos,locus_idx in ilist:
                _ret[self._locus[locus_idx]] = epos-spos
        return _ret

    def intersect_blocks(self, ref, blocks):
        raise NotImplementedError()
        _result = Counter()
        for b_start, b_end in blocks:
            pass
        return _result

def overlap_length(a,b):
    return max(0, min(a.end,b.end) - max(a.begin,b.begin))


def merge_intervals(a, b, d=None):
    return Interval(min(a.begin,b.begin), max(a.end,b.end), d)

class _AnnotationIntervalTree(object):

    def __init__(self, gtf_file, attribute_name):
        self.loci = OrderedDict()
        self.key = attribute_name
        self.itree = defaultdict(IntervalTree)

        # GTF filehandle
        fh = open(gtf_file,'rU') if isinstance(gtf_file,str) else gtf_file
        features = (GTFRow(*l.strip('\n').split('\t')) for l in fh if not l.startswith('#'))
        for f in features:
            attr = dict(re.findall('(\w+)\s+"(.+?)";', f.attribute))
            ''' Add to locus list '''
            if attr[self.key] not in self.loci:
                self.loci[attr[self.key]] = list()
            self.loci[attr[self.key]].append(f)
            ''' Add to interval tree '''
            new_iv = Interval(int(f.start), int(f.end)+1, attr)
            # Merge overlapping intervals from same locus
            if True:
                overlap = self.itree[f.chrom][new_iv]
                if len(overlap) > 0:
                    mergeable = [iv for iv in overlap if iv.data[self.key]==attr[self.key]]
                    if mergeable:
                        assert len(mergeable) == 1, "Error"
                        new_iv = merge_intervals(mergeable[0], new_iv, {self.key: attr[self.key]})
                        self.itree[f.chrom].remove(mergeable[0])
            self.itree[f.chrom].add(new_iv)

    def feature_length(self):
        """ Get feature lengths

        Returns:
            (dict of str: int): Feature names to feature lengths

        """
        ret = Counter()
        for chrom in list(self.itree.keys()):
            for iv in list(self.itree[chrom].items()):
                ret[iv.data[self.key]] += iv.length()
        return ret

    def intersect_blocks(self, ref, blocks):
        _result = Counter()
        for b_start, b_end in blocks:
            query = Interval(b_start, (b_end + 1))
            for iv in self.itree[ref][query]:
                _result[iv.data[self.key]] += overlap_length(iv, query)
        return _result


class _AnnotationHTSeq(object):
    def __init__(self, gtf_file, attribute_name):
        self.loci = OrderedDict()
        self.features =  HTSeq.GenomicArrayOfSets( "auto", stranded=False )
        for f in HTSeq.GFF_Reader(gtf_file, end_included = True ):
            if f.type == 'exon':
                self.features[f.iv] += f.attr[attribute_name]
                if f.attr[attribute_name] not in self.loci:
                    self.loci[f.attr[attribute_name]] = list()
                self.loci[f.attr[attribute_name]].append(f)

    def feature_length(self):
        """ Get feature lengths

        Returns:
            (dict of str: int): Feature names to feature lengths

        """
        ret = Counter()
        for ref_iv, step_set in self.features.steps():
            for feat_id in step_set:
                ret[feat_id] += ref_iv.length
        return ret

    def intersect_blocks(self, ref, blocks):
        _result = Counter()
        for b_start, b_end in blocks:
            iv = HTSeq.GenomicInterval(ref, b_start, b_end)
            for ref_iv, step_set in self.features[iv].steps():
                for feat_id in step_set:
                    _result[feat_id] += ref_iv.length
        return _result


ANNOTATION_CLASS='intervaltree'

# Import Annotation class based on ANNOTATION_CLASS
if ANNOTATION_CLASS == 'bisect':
    Annotation = _AnnotationBisect
elif ANNOTATION_CLASS == 'intervaltree':
    from intervaltree import Interval, IntervalTree
    Annotation = _AnnotationIntervalTree
elif ANNOTATION_CLASS == 'htseq':
    import HTSeq
    Annotation = _AnnotationHTSeq

else:
    raise ImportError('Could not import Annotation "%s"' % ANNOTATION_CLASS)
