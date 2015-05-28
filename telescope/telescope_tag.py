__author__ = 'bendall'

import os, sys

import pysam

import utils
from utils.alignment_parsers import TelescopeRead
from utils.annotation_parsers import AnnotationLookup


class TagOpts:
    option_fields = ['verbose', 'gtffile', 'samfile', 'outfile', ]

    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            # if v is None: continue
            setattr(self, k, v)

    def __str__(self):
        _ret = "TagOpts:\n"
        _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
        return _ret

def run_telescope_tag(args):
    opts = TagOpts(**vars(args))

    if opts.verbose:
        print >>sys.stderr, opts

    samfile = pysam.AlignmentFile(opts.samfile)
    refnames = dict(enumerate(samfile.references))

    has_features = opts.gtffile is not None
    if has_features:
        flookup = AnnotationLookup(opts.gtffile)

    outfile = pysam.AlignmentFile(opts.outfile, 'wh', header=samfile.header)

    for rname,segments in utils.iterread(samfile):
        r = TelescopeRead(rname,segments)
        if has_features: r.assign_feats(refnames, flookup)
        r.set_color_tags()
        r.set_optional_tags()
        for a in r.alignments:
            a.write_samfile(outfile)

    samfile.close()
    outfile.close()
