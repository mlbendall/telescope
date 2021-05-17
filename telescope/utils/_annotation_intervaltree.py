# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import absolute_import
from builtins import object

import re
from collections import defaultdict, namedtuple, Counter, OrderedDict
import logging as lg
import pickle


from intervaltree import Interval, IntervalTree


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


GTFRow = namedtuple('GTFRow', ['chrom','source','feature','start','end','score','strand','frame','attribute'])

def overlap_length(a,b):
    return max(0, min(a.end,b.end) - max(a.begin,b.begin))

def merge_intervals(a, b, d=None):
    return Interval(min(a.begin,b.begin), max(a.end,b.end), d)

class _AnnotationIntervalTree(object):

    def __init__(self, gtf_file, attribute_name, feature_type='exon'):
        lg.debug('Using intervaltree for annotation.')
        self.loci = OrderedDict()
        self.key = attribute_name
        self.itree = defaultdict(IntervalTree)

        # GTF filehandle
        fh = open(gtf_file,'rU') if isinstance(gtf_file,str) else gtf_file
        for rownum, l in enumerate(fh):
            if l.startswith('#'): continue
            f = GTFRow(*l.strip('\n').split('\t'))
            if f.feature != feature_type: continue
            attr = dict(re.findall('(\w+)\s+"(.+?)";', f.attribute))
            if self.key not in attr:
                lg.warning('Skipping row %d: missing attribute "%s"' % (rownum, self.key))
                continue

            ''' Add to locus list '''
            if attr[self.key] not in self.loci:
                self.loci[attr[self.key]] = list()
            self.loci[attr[self.key]].append(f)
            ''' Add to interval tree '''
            new_iv = Interval(int(f.start), int(f.end)+1, attr)
            # Merge overlapping intervals from same locus
            if True:
                overlap = self.itree[f.chrom].overlap(new_iv)
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

    def subregion(self, ref, start_pos=None, end_pos=None):
        _subannot = type(self).__new__(type(self))
        _subannot.key = self.key
        _subannot.itree = defaultdict(IntervalTree)

        if ref in self.itree:
            _subtree = self.itree[ref].copy()
            if start_pos is not None:
                _subtree.chop(_subtree.begin(), start_pos)
            if end_pos is not None:
                _subtree.chop(end_pos, _subtree.end() + 1)
            _subannot.itree[ref] = _subtree
        return _subannot

    def intersect_blocks(self, ref, blocks):
        _result = Counter()
        for b_start, b_end in blocks:
            query = Interval(b_start, (b_end + 1))
            for iv in self.itree[ref].overlap(query):
                _result[iv.data[self.key]] += overlap_length(iv, query)
        return _result

    def save(self, filename):
        with open(filename, 'wb') as outh:
            pickle.dump({
                'key': self.key,
                'loci': self.loci,
                'itree': self.itree,
            }, outh)

    @classmethod
    def load(cls, filename):
        with open(filename, 'rb') as fh:
            loader = pickle.load(fh)
        obj = cls.__new__(cls)
        obj.key = loader['key']
        obj.loci = loader['loci']
        obj.itree = loader['itree']

        return obj
