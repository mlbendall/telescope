# -*- coding: utf-8 -*-
""" Telescope id

"""

import sys
import os
from time import time

import pysam

import utils
from utils.alignment_parsers import TelescopeRead
from utils.model import TelescopeModel
from utils.colors import c2str, DARK2_PALETTE, GREENS
from utils import format_minutes

from utils.annotation_parsers import Annotation

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


class IDOpts:
    option_fields = ['ali_format','samfile','gtffile',
                     'verbose', 'outdir', 'exp_tag', 'out_matrix', 'updated_sam',
                     'checkpoint', 'checkpoint_interval',
                     'min_prob',
                     'conf_prob', 'reassign_mode',
                     'piPrior', 'thetaPrior', 'min_overlap',
                     'emEpsilon','maxIter',
                     'version',
                     ]

    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            if v is None: continue
            setattr(self, k, v)

    def generate_filename(self,suffix):
        basename = '%s-%s' % (self.exp_tag, suffix)
        return os.path.join(self.outdir,basename)

    def __str__(self):
        _ret = "IDOpts:\n"
        _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
        return _ret


def load_alignment(samfile, flookup, opts=None):
    _verbose = opts.verbose if opts is not None else True

    counts = {'unmapped':0, 'nofeat':0, 'ambig':0, 'unique':0}
    mapped   = {}

    # Lookup reference name from reference ID
    refnames = dict(enumerate(samfile.references))

    for rname,segments in utils.iterread(samfile):
        r = TelescopeRead(rname,segments)
        if r.is_unmapped:
            counts['unmapped'] += 1
        else:
            r.assign_feats(refnames, flookup)
            if r.aligns_to_feat():
                r.assign_best()
                mapped[rname] = r
                if r.is_unique:
                    counts['unique'] += 1
                else:
                    counts['ambig'] += 1
            else:
                counts['nofeat'] += 1

        if _verbose and sum(counts.values()) % 500000 == 0:
            print >>sys.stderr, "...Processed %d fragments" % sum(counts.values())

    return mapped, counts


def update_alignment(tm, mapped, newsam, opts):
    # Read x Transcript matrix = 1 if output alignment 0 otherwise
    # output_mat = tm.x_hat.apply_func(lambda x: 1 if x >= opts.min_prob else 0)
    output_mat = tm.reassign_to_best(opts.reassign_mode)
    for rownum in xrange(output_mat.shape[0]):
        # Continue if there are no alignments to output
        if output_mat[rownum,].maxr()[0,0] == 0: continue
        _rname = tm.readnames[rownum]
        _read  = mapped[_rname]
        try:
            # Sorted list of (transcript_index, prob)
            tidx_probs = sorted(((_, tm.x_hat[rownum,_]) for _ in output_mat[rownum,].nonzero()[1]), key=lambda x:x[1], reverse=True)
            best_transcript = tm.txnames[tidx_probs[0][0]]
            for colnum,prob in tidx_probs:
                transcript_name = tm.txnames[colnum]
                primary, alternates = _read.aligned_to_transcript(transcript_name)

                # Print primary alignment
                primary.set_mapq(utils.phred(prob))
                primary.set_tag('XP',int(round(prob*100)))
                primary.set_tag('XT', transcript_name)
                primary.set_tag('ZT', best_transcript)

                if transcript_name == best_transcript:            # best transcript
                    primary.set_secondary(False)
                    if len(tidx_probs)==1:                    # only one transcript has prob > min_prob
                        if prob >= opts.conf_prob:                     # high confidence
                            primary.set_tag('YC', c2str(DARK2_PALETTE['vermilion']))
                        else:                                     # low confidence
                            primary.set_tag('YC', c2str(DARK2_PALETTE['yellow']))
                    else:                                     # multiple transcripts have prob > min_prob
                        primary.set_tag('YC', c2str(DARK2_PALETTE['teal']))
                        assert prob < .9, "If there are multiple nonzero transcripts, qual must be < .9"
                else:
                    primary.set_tag('YC', c2str(GREENS[2]))    # Not best genome
                    primary.set_secondary(True)

                primary.write_samfile(newsam)

                # Print alternate alignments
                for altaln in alternates:
                    altaln.set_mapq(0)
                    altaln.set_tag('XP',0)
                    altaln.set_tag('XT', transcript_name)
                    altaln.set_tag('ZT', best_transcript)
                    altaln.set_tag('YC', c2str((248,248,248)))
                    altaln.set_secondary(True)
                    altaln.write_samfile(newsam)
        except IOError as e:
            print >>sys.stderr, e
            print >>sys.stderr, "Unable to write %s" % _rname


def make_report(tm, aln_counts, txlens, opts, sortby='final_count'):
    # Body of report has values for each transcript
    report_fmt = [('transcript','%s'),('transcript_length','%d'),
                  ('final_count','%d'), ('final_conf','%d'), ('final_prop','%.6g'),
                  ('unique_count','%d'), ('init_aligned','%d'),
                  ('init_best','%d'), ('init_best_random','%d'), ('init_best_avg','%.6g'),
                  ('init_prop','%.6g'),
                 ]
    columns = {}
    columns['transcript']  = tm.txnames
    columns['transcript_length'] = [txlens[tx] for tx in tm.txnames]
    columns['final_count'] = tm.reassign_to_best(opts.reassign_mode).sumc().A1
    columns['final_conf']  =  tm.reassign_to_best('conf', thresh=opts.conf_prob).sumc().A1
    columns['final_prop']  =  tm.pi

    # Number of unambiguous alignments
    columns['unique_count'] = tm.reassign_to_best('unique').sumc().A1

    # Initial number of alignments
    columns['init_aligned'] = tm.x_init.ceil().sumc().A1
    # Initial number of best alignments
    columns['init_best'] = tm.reassign_to_best('exclude', initial=True).sumc().A1
    columns['init_best_random'] = tm.reassign_to_best('choose', initial=True).sumc().A1
    columns['init_best_avg'] = tm.reassign_to_best('average', initial=True).sumc().A1
    # Initial proportions
    columns['init_prop'] = tm.pi_0

    R,T = tm.shape
    colheader = [h[0] for h in report_fmt]
    _fmt = [[f % columns[n][j] for n,f in report_fmt] for j in range(T)]
    _fmt.sort(key=lambda x: float(x[colheader.index(sortby)]), reverse=True)

    # Commented section of the report contains some overall run metrics
    runinfo = aln_counts
    runinfo['transcripts']       = len(tm.txnames)
    runinfo['telescope_version'] = opts.version
    comment = ['## RunInfo'] + ['%s:%s' % (k,v) for k,v in runinfo.iteritems()]
    # comment = ['## RunInfo'] + ['%s:%d' % (k,v) for k,v in aln_counts.iteritems()] + ['transcripts:%d' % len(tm.txnames)]

    return [comment, colheader] + _fmt


def run_telescope_id(args):
    opts = IDOpts(**vars(args))
    if opts.verbose:
        print >>sys.stderr, opts

    """ Load annotation """
    if opts.verbose:
        print >>sys.stderr, "Loading annotation file (%s):" % opts.gtffile

    flookup = Annotation(opts.gtffile, min_overlap=opts.min_overlap)

    """ Load alignment """
    if opts.verbose:
        print >>sys.stderr, "Loading alignment file (%s):" % opts.samfile
        substart = time()

    samfile = pysam.AlignmentFile(opts.samfile)
    mapped, aln_counts = load_alignment(samfile, flookup, opts)

    if opts.verbose:
        print >>sys.stderr, "Time to load alignment:".ljust(40) + format_minutes(time() - substart)
        print >>sys.stderr, "Alignment counts:"
        print >>sys.stderr, "Processed %d fragments" % sum(aln_counts.values())
        print >>sys.stderr, "\t%d fragments were unmapped" % aln_counts['unmapped']
        print >>sys.stderr, "\t%d fragments mapped to one or more positions on reference genome" % (aln_counts['unique'] + aln_counts['ambig'] + aln_counts['nofeat'])
        print >>sys.stderr, "\t\t%d fragments mapped to reference but did not map to any transcripts" % aln_counts['nofeat']
        print >>sys.stderr, "\t\t%d fragments have at least one alignment to a transcript" % (aln_counts['unique'] + aln_counts['ambig'])
        print >>sys.stderr, "\t\t\t%d fragments align uniquely to a single transcript" % aln_counts['unique']
        print >>sys.stderr, "\t\t\t%d fragments align ambiguously to multiple transcripts" % aln_counts['ambig']

    """ Create data structure """
    if opts.verbose:
        print >>sys.stderr, "Creating data structure... " ,
        substart = time()

    ridx = {}
    gidx = {}
    d = []
    for rname,r in mapped.iteritems():
        i = ridx.setdefault(rname,len(ridx))
        for gname,aln in r.feat_aln_map.iteritems():
            j = gidx.setdefault(gname,len(gidx))
            d.append((i, j, aln.AS + aln.query_length))

    tm = TelescopeModel(ridx, gidx, data=d)

    if opts.verbose:
        print >>sys.stderr, "done."
        print >>sys.stderr, "Time to create data structure:".ljust(40) + format_minutes(time() - substart)

    # Save some memory if you are not creating an updated SAM:
    if not opts.updated_sam:
        if opts.verbose:
            print >>sys.stderr, "Clearing reads from memory."
        mapped = None
        samfile.close()

    """ Initial checkpoint """
    if opts.checkpoint:
        if opts.verbose:
            print >>sys.stderr, "Checkpointing... " ,
            substart = time()

        with open(opts.generate_filename('checkpoint.init.p'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write checkpoint:".ljust(40) + format_minutes(time() - substart)

    """ Initial output matrix """
    if opts.out_matrix:
        if opts.verbose:
            print >>sys.stderr, "Writing initial model matrix...",
            substart = time()

        with open(opts.generate_filename('matrix_init.p'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write initial matrix:".ljust(40) + format_minutes(time() - substart)

    """ Reassignment """
    if opts.verbose:
        print >>sys.stderr, "Reassiging reads:"
        print >>sys.stderr, "(Reads,Transcripts)=%dx%d" % (tm.shape)
        print >>sys.stderr, "Delta Change:"
        substart = time()

    tm.matrix_em(opts)

    if opts.verbose:
        print >>sys.stderr, "Time for EM iteration:".ljust(40) + format_minutes(time() - substart)

    """ Checkpoint 2 """
    if opts.checkpoint:
        if opts.verbose:
            print >>sys.stderr, "Checkpointing... ",
            substart = time()

        with open(opts.generate_filename('checkpoint.final.p'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write checkpoint:".ljust(40) + format_minutes(time() - substart)

    """ Final output matrix """
    if opts.out_matrix:
        if opts.verbose:
            print >>sys.stderr, "Writing final model matrix...",
            substart = time()

        with open(opts.generate_filename('matrix_final.p'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write final matrix:".ljust(40) + format_minutes(time() - substart)

    """ Generate report """
    if opts.verbose:
        print >>sys.stderr, "Generating report... " ,
        substart = time()

    # report = tm.make_report(opts.conf_prob)
    report = make_report(tm, aln_counts, flookup.feature_length(), opts)
    # comment = ['%s:%d' % (k,v) for k,v in aln_counts.iteritems()] + ['transcripts:%d' % len(tm.txnames)]
    with open(opts.generate_filename('telescope_report.tsv'),'w') as outh:
        # print >>outh, '## RunInfo\t%s' % '\t'.join(comment)
        for row in report:
            print >>outh, '\t'.join(f for f in row)

    if opts.verbose:
        print >>sys.stderr, "done."
        print >>sys.stderr, "Time to generate report:".ljust(40) + format_minutes(time() - substart)

    """ Update alignment """
    if opts.updated_sam:
        if opts.verbose:
            print >>sys.stderr, "Updating alignment...",
            substart = time()

        updated_samfile = pysam.AlignmentFile(opts.generate_filename('updated.sam'), 'wh', header=samfile.header)
        update_alignment(tm, mapped, updated_samfile, opts)
        updated_samfile.close()

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to update alignment:".ljust(40) + format_minutes(time() - substart)

    samfile.close()
    return
