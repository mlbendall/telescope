# -*- coding: utf-8 -*-
""" Telescope TAG

"""

import os
import sys

import pysam

import utils
from utils.alignment_parsers import TelescopeRead
from utils.colors import c2str, DARK2_PALETTE, GREENS

from utils.annotation_parsers import Annotation

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


class TagOpts:
    option_fields = ['verbose', 'gtffile', 'samfile', 'outfile', 'min_overlap', 'version', ]

    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            # if v is None: continue
            setattr(self, k, v)

    def __str__(self):
        _ret = "TagOpts:\n"
        _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
        return _ret


def set_color_tags(telescope_read):
    if telescope_read.is_unique:
        telescope_read.alignments[0].set_tag('YC', c2str(DARK2_PALETTE['vermilion']))
    else:
        for a in telescope_read.alignments:
            if a.AS == telescope_read.bestAS:
                a.set_tag('YC', c2str(DARK2_PALETTE['teal']))
            else:
                pct = float(a.AS) / telescope_read.bestAS
                if pct > 0.95:
                    a.set_tag('YC', c2str(GREENS[0]))
                elif pct > 0.90:
                    a.set_tag('YC', c2str(GREENS[1]))
                elif pct > 0.85:
                    a.set_tag('YC', c2str(GREENS[2]))
                else:
                    a.set_tag('YC', c2str(GREENS[3]))


def set_optional_tags(telescope_read):
    ''' Sets optional tags for read
            XC:i  Alignment count - total number of alignments for this read
            XT:Z  Alignment transcript - name of the transcript containing this alignment
            ZC:i  Best count - number of "best" alignments sharing the same maximum score
            ZS:i  Best score - maximum alignment score for this read
            ZT:Z  Best transcript - names(s) of transcript(s) containing best alignment(s)
    '''
    num_best  = sum(a.AS==telescope_read.bestAS for a in telescope_read.alignments)
    tags = [
             ('XC', len(telescope_read.alignments)),
             ('ZC', num_best),
             ('ZS', telescope_read.bestAS),
           ]
    if telescope_read.features:
        bestfeats = set([f for a,f in zip(telescope_read.alignments, telescope_read.features) if a.AS == telescope_read.bestAS])
        tags.append(('ZT',','.join(sorted(bestfeats))))
    for i,a in enumerate(telescope_read.alignments):
        if telescope_read.features:
            a.set_tags(tags + [('XT', telescope_read.features[i])])
        else:
            a.set_tags(tags)

def run_telescope_tag(args):
    opts = TagOpts(**vars(args))

    if opts.verbose:
        print >>sys.stderr, opts

    samfile = pysam.AlignmentFile(opts.samfile)
    refnames = dict(enumerate(samfile.references))

    has_features = opts.gtffile is not None
    if has_features:
        flookup = Annotation(opts.gtffile, min_overlap=opts.min_overlap)

    outfile = pysam.AlignmentFile(opts.outfile, 'wh', header=samfile.header)

    for rname,segments in utils.iterread(samfile):
        r = TelescopeRead(rname,segments)
        if not r.is_unmapped:
            if has_features:
                r.assign_feats(refnames, flookup)
                # r.assign_feats(refnames, annotation)
                #for aln,feat in zip(r.alignments,r.features):
                #    c,s,e = aln.coordinates()
                #    foundfeat = annotation.lookup_interval(refnames[c], s, e)
                #    if foundfeat is None:
                #        if feat != TelescopeRead.nofeature:
                #            print >>sys.stderr, '%s %s %s' % (feat, foundfeats, aln.coordinates())
                #    else:
                #        print >>sys.stderr, '%s %s' % (feat, foundfeat)
            set_color_tags(r)
            set_optional_tags(r)

        for a in r.alignments:
            a.write_samfile(outfile)

    samfile.close()
    outfile.close()
