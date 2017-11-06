# -*- coding: utf-8 -*-
""" Telescope id

"""
from __future__ import print_function
from __future__ import absolute_import

from builtins import range
from builtins import object
from builtins import super
import sys
import os
from time import time
import logging as lg
from collections import OrderedDict, defaultdict, Counter
from operator import itemgetter

import pysam
import numpy as np

from . import utils
# from .utils.alignment_parsers import TelescopeRead
# from .utils import alignment

USE_CYTHON = True
if USE_CYTHON:
    from telescope.cTelescope import load_unsorted
else:
    from .utils.bam_parsers import load_unsorted

from .utils.model import TelescopeModel
from .utils.colors import c2str, DARK2_PALETTE, GREENS
from .utils import format_minutes as fmtmins

from .utils.annotation_parsers import Annotation
from .utils.annotation_parsers import ANNOTATION_CLASS

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


class IDOptions(utils.SubcommandOptions):
    OPTS = """
    - Input Options:
        - samfile:
            positional: True
            help: Path to alignment file. Alignment file can be in SAM or BAM
                  format. File must be collated so that all alignments for a
                  read pair appear sequentially in the file.
        - gtffile:
            positional: True
            help: Path to annotation file (GTF format)
        - attribute:
            default: locus
            help: GTF attribute that defines a transposable element locus. GTF
                  features that share the same value for --attribute will be
                  considered as part of the same locus.
        - ncpu:
            default: 1
            type: int
            help: Number of cores to use. (Multiple cores not supported yet).
        - no_feature_key:
            default: __no_feature
            help: Used internally to represent alignments. Must be different
                  from all other feature names.
        - no_feature_key:
            default: __no_feature
            help: Used internally to represent alignments. Must be different
                  from all other feature names.
    - Reporting Options:
        - quiet:
            action: store_true
            help: Silence (most) output.
        - debug:
            action: store_true
            help: Print debug messages.
        - logfile:
            type: argparse.FileType('r')
            help: Log output to this file.
        - outdir:
            default: .
            help: Output directory.
        - exp_tag:
            default: telescope
            help: Experiment tag
        - out_matrix:
            action: store_true
            help: Output alignment matrix
        - updated_sam:
            action: store_true
            help: Generate an updated alignment file.
    - Run Modes:
        - reassign_mode:
            default: exclude
            choices:
                - exclude
                - choose
                - average
                - conf
                - unique
            help: >
                  Reassignment mode. After EM is complete, each fragment is
                  reassigned according to the expected value of its membership
                  weights. The reassignment method is the method for resolving
                  the "best" reassignment for fragments that have multiple
                  possible reassignments.
                  Available modes are: "exclude" - fragments with multiple best
                  assignments are excluded from the final counts; "choose" -
                  the best assignment is randomly chosen from among the set of
                  best assignments; "average" - the fragment is divided evenly
                  among the best assignments; "conf" - only assignments that
                  exceed a certain threshold (see --conf_prob) are accepted;
                  "unique" - only uniquely aligned reads are included.
                  NOTE: Results using all assignment modes are included in the
                  Telescope report by default. This argument determines what
                  mode will be used for the "final counts" column.
        - conf_prob:
            type: float
            default: 0.9
            help: Minimum probability for high confidence assignment.
        - overlap_mode:
            default: threshold
            choices:
                - threshold
                - intersection-strict
                - union
            help: Overlap mode. The method used to determine whether a fragment
                  overlaps feature.
        - overlap_threshold:
            type: float
            default: 0.2
            help: Fraction of fragment that must be contained within a feature
                  to be assigned to that locus. Ignored if --overlap_method is
                  not "threshold".
        - bootstrap:
            hide: True
            type: int
            help: Set to an integer > 0 to turn on bootstrapping. Number of
                  bootstrap replicates to perform.
        - bootstrap_ci:
            hide: True
            type: float
            default: 0.95
            help: Size of bootstrap confidence interval
    - Model Parameters:
        - pi_prior:
            type: int
            default: 0
            help: Prior on π. Equivalent to adding n unique reads.
        - theta_prior:
            type: int
            default: 0
            help: Prior on θ. Equivalent to adding n non-unique reads.
        - em_epsilon:
            type: float
            default: 1e-7
            help: EM Algorithm Epsilon cutoff
        - max_iter:
            type: int
            default: 100
            help: EM Algorithm maximum iterations
    """

    def __init__(self, args):
        super().__init__(args)
        if self.logfile is None:
            self.logfile = sys.stderr

    def generate_filename(self, suffix):
        basename = '%s-%s' % (self.exp_tag, suffix)
        return os.path.join(self.outdir, basename)


class IDOpts(object):
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
        for k,v in kwargs.items():
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

        if sum(counts.values()) % 500000 == 0:
            lg.info("...Processed %d fragments" % sum(counts.values()))

    return mapped, counts


def update_alignment(tm, mapped, newsam, opts):
    # Read x Transcript matrix = 1 if output alignment 0 otherwise
    # output_mat = tm.x_hat.apply_func(lambda x: 1 if x >= opts.min_prob else 0)
    output_mat = tm.reassign_to_best(opts.reassign_mode)
    for rownum in range(output_mat.shape[0]):
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
            print(e, file=sys.stderr)
            print("Unable to write %s" % _rname, file=sys.stderr)


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
    comment = ['## RunInfo'] + ['%s:%s' % (k,v) for k,v in runinfo.items()]
    # comment = ['## RunInfo'] + ['%s:%d' % (k,v) for k,v in aln_counts.iteritems()] + ['transcripts:%d' % len(tm.txnames)]

    return [comment, colheader] + _fmt

def configure_logging(opts):
    loglev = lg.INFO
    if opts.quiet:
        loglev = lg.WARNING
    if opts.debug:
        loglev = lg.DEBUG

    logfmt = '%(asctime)s %(levelname)-8s %(message)-50s'
    logfmt += ' (from %(funcName)s in %(filename)s:%(lineno)d)'
    lg.basicConfig(level=loglev,
                        format=logfmt,
                        datefmt='%m-%d %H:%M',
                        stream=opts.logfile)


import scipy
import pickle
from .utils.helpers import grouper



class Telescope(object):
    def __init__(self, opts):
        self.opts       = opts
        self.read_index = {}
        self.readnames  = []
        self.feat_index = {}
        self.featnames  = []
        self.shape      = None
        self.readkeys   = {}

        self.raw_scores = None
        self.annotation = None

        # Results
        self.final_counts = None
        self.final_conf   = None
        self.final_prop   = None
        self.bs_counts = None
        self.bs_props = None

        self.run_info = OrderedDict()
        self.run_info['version'] = self.opts.version

        # Load initial information about the alignment file
        with pysam.AlignmentFile(opts.samfile) as sf:
            self.references = sf.references
            self.reflens = sf.lengths
            if sf.has_index():
                self.is_sorted = True
                self.run_info['alignments'] = sf.mapped
                self.run_info['unaligned'] = sf.unmapped
            else:
                self.is_sorted = False

    def load_annotation(self):
        lg.info('Loading annotation')
        self.annotation = Annotation(self.opts.gtffile, self.opts.attribute)
        lg.info('Loaded {} features'.format(len(self.annotation.loci)))
        self.run_info['annotated_features'] = len(self.annotation.loci)

    def load_alignment(self):
        alninfo = {'unmap_1': 0, 'unmap_2': 0,
                   'map_1': 0, 'map_2': 0,
                   'nofeat_U': 0, 'nofeat_A': 0,
                   'feat_U': 0, 'feat_A': 0,
                   'minAS': float('inf'),
                   'maxAS': float('-inf'),
                   }

        mappings = {}

        # Load unsorted reads
        iter = load_unsorted(self.opts.samfile, self.annotation, self.opts, alninfo)

        for pinfo in iter:
            _key = (pinfo[0], pinfo[1])
            if _key not in mappings:
                mappings[_key] = (pinfo[2], pinfo[3])
            else:
                if mappings[_key][0] < pinfo[2]:
                    mappings[_key] = (pinfo[2], pinfo[3])

        input_pairs = alninfo['map_2'] + alninfo['unmap_2'] + int(float(alninfo['map_1'] + alninfo['unmap_1']) / 2)
        total_frags = sum(alninfo[k] for k in ['unmap_1','unmap_2','map_1','map_2'])
        lg.info('Processed {} fragments.'.format(total_frags))
        lg.info('Input pairs: {}.'.format(input_pairs))
        for k, v in alninfo.items():
            lg.info('{:30}{}'.format(k + ':', v))

        _keys = sorted(mappings.keys(), key=itemgetter(0))
        tmpfn = 'tmp0.{}.txt'.format(ANNOTATION_CLASS)
        lg.info('annotation class: {}'.format(ANNOTATION_CLASS))
        with open(tmpfn, 'w') as outh:
            for k in _keys:
                print('id: {}; feat: {}; score: {}; len: {}'.format(*k, *mappings[k]), file=outh)

    def xload_alignment(self):
        alignment.load_unsorted(self.opts.samfile, self.annotation, self.opts)
        self.run_info['fragments'] = 0
        self.run_info['unmapped'] = 0
        self.run_info['alignments'] = 0
        self.run_info['nofeat'] = 0
        self.run_info['ambig'] = 0
        self.run_info['unique'] = 0
        self.run_info['maxAS'] = float('-inf')
        self.run_info['minAS'] = float('inf')

        _nfkey = self.opts.no_feature_key

        overlap_func = alignment.get_overlap_function(
            self.annotation,
            self.opts.overlap_mode,
            self.opts.overlap_threshold,
            _nfkey
        )

        with pysam.AlignmentFile(self.opts.samfile) as sf:
            for bundle in alignment.fetch_bundle_pairs(sf):
                self.run_info['fragments'] += 1
                pairs = list(bundle)

                mpairs = [p for p in pairs if not p[0].is_unmapped]
                if not mpairs:
                    self.run_info['unmapped'] += 1
                    continue

                overlap_feats = list(map(overlap_func, mpairs))
                has_overlap = any(f != _nfkey for f in overlap_feats)
                if not has_overlap:
                    self.run_info['nofeat'] += 1
                    continue

                d = {}
                for p,f in zip(mpairs, overlap_feats):
                    k = (alignment.fragid(p), f)
                    d[k] = (max(d[k], alignment.alnscore(p)), alignment.alnlen(p))

                self.run_info['maxAS'] = max(self.run_info['maxAS'], *d.values())
                if len(d) == 1:
                    self.run_info['unique'] += 1
                else:
                    self.run_info['ambig'] += 1

                dim = (1e8, self.run_info['annotated_features'])
                _m1 = scipy.sparse.dok_matrix(dim, dtype=np.int16)
                _ridx = {}
                _fidx = {}

                for (read_id, feat_id), ascore in d.items():
                    i = _ridx.setdefault(read_id, len(_ridx))
                    j = _fidx.setdefault(feat_id, len(_fidx))


    def _load_alignment(self):
        if self.is_sorted:
            if self.opts.ncpu <= 1:
                lg.info('Reading sorted alignment sequentially...')
                no_region = (None, None, None)
                _raw_alignments, metrics = bparse.load_region(
                    self.opts.samfile, self.annotation, self.opts, no_region
                )
            else:
                lg.info('Reading sorted alignment in parallel...')
                regions = region_iter(self.references, self.reflens)
                pool = Pool(processes=self.opts.ncpu)
                result = pool.map_async(
                    partial(bparse.load_region,
                            self.opts.samfile, self.annotation, self.opts
                            ),
                    regions
                )
                # Combine map results:
                _raw_alignments = []
                metrics = {'minAS': float('inf'), 'maxAS': float('-inf'),
                           'aligned_pairs': 0, 'aligned_singletons': 0,
                           }
                for raln, met in result.get():
                    _raw_alignments.extend(raln)
                    metrics['minAS'] = min(metrics['minAS'], met['minAS'])
                    metrics['maxAS'] = max(metrics['maxAS'], met['maxAS'])
                    metrics['aligned_pairs'] += met['aligned_pairs']
                    metrics['aligned_singletons'] += met['aligned_singletons']
                pool.close()
        else:
            lg.info('Reading unsorted alignment sequentially...')
            _raw_alignments, metrics = bparse.load_unsorted(
                self.opts.samfile, self.annotation, self.opts
            )
            self.run_info['alignments'] = metrics['num_map']
            self.run_info['unaligned'] = metrics['num_unmap']

        self.run_info['aligned_frags'] = len(_raw_alignments)

        c = Counter(len(ra[1]) for ra in _raw_alignments)
        assert c[1] == metrics['aligned_singletons']
        assert c[2] == metrics['aligned_pairs']

        self.run_info['aligned_singletons'] = metrics['aligned_singletons']
        self.run_info['aligned_pairs'] = metrics['aligned_pairs']
        scaleAS = metrics['minAS']
        # Create sparse matrix where the value m1[i, j] is the score for read i
        # aligning to feature j.
        # The size of the matrix must be larger than the expected data
        dim = (len(_raw_alignments), len(self.annotation.features())+1)
        _m1 = scipy.sparse.dok_matrix(dim, dtype=np.int)
        _ridx = {}
        _fidx = {}
        _rkeys = defaultdict(dict)

        # minscore = min(alndata.ascore for alndata, rkey in _raw_alignments)
        # maxscore = max(alndata.ascore for alndata, rkey in _raw_alignments)
        # assert minscore == alndata['minAS']
        # assert maxscore == alndata['maxAS']

        # rescaledmin = min(alndata.ascore + alndata.alen for alndata, rkey in _raw_alignments)
        # lg.debug('Rescaled min: %d' % rescaledmin)
        # rescaledmin = min(alndata.ascore - scaleAS + 1 for alndata, rkey in _raw_alignments)
        # lg.debug('Rescaled min: %d'  % rescaledmin)

        for metrics, rkey in _raw_alignments:
            i = _ridx.setdefault(metrics.readid, len(_ridx))
            j = _fidx.setdefault(metrics.featid, len(_fidx))
            # The raw score is the alignment score scaled to the minimum alignment score plus the alignment length
            score = metrics.ascore - scaleAS + 1 + metrics.alen
            assert score > 0
            if _m1[i, j] < score:
                _m1[i, j] = score
                _rkeys[i][j] = set([rkey])
            elif _m1[i, j] == score:
                _rkeys[i][j].add(rkey)

        # Trim matrix down to size of data
        _m1 = _m1[:len(_ridx), :len(_fidx)]
        _raw_alignments = None
        self.run_info['fragments'] = len(_ridx)

        # Remove reads that only align to "__nofeature__"
        # Indexes for reads that align to a feature
        _withfeat = sorted(set(r for r,c in zip(*_m1.nonzero())
                               if c != _fidx[self.opts.no_feature_key]))
        _rev_ridx = sorted(_ridx, key=_ridx.get) # Lookup read name given index
        self.raw_scores = csr_matrix(_m1[_withfeat, :])
        for newrow, oldrow in enumerate(_withfeat):
            self.read_index[_rev_ridx[oldrow]] = newrow
            self.readkeys[newrow] = _rkeys[oldrow]

        self.readnames = sorted(self.read_index, key=self.read_index.get)
        self.feat_index = _fidx
        self.featnames = sorted(self.feat_index, key=self.feat_index.get)
        self.shape = (len(self.readnames), len(self.featnames))

        alncounts = self.raw_scores.countr()
        self.run_info['unique_TE_frags'] = sum(alncounts == 1)
        self.run_info['ambiguous_TE_frags'] = sum(alncounts > 1)
        self.run_info['annotated_features'] = len(self.annotation.features())
        self.run_info['mapped_features'] = len(self.featnames)

        # Print run information
        run_info_str = ["Run Info:"]
        for k, v in self.run_info.iteritems():
            run_info_str.append('\t%s%s' % (k.ljust(30), v))
        lg.info('\n'.join(run_info_str))

    def reassign(self, tl):
        rmethod, rprob = self.opts.reassign_mode, self.opts.conf_prob
        if rmethod == 'conf':
            self.final_counts = tl.reassign(rmethod, rprob).sum(0).A1
        else:
            self.final_counts = tl.reassign(rmethod).sum(0).A1

        self.final_conf = tl.reassign('conf', rprob).sum(0).A1
        self.final_prop = tl.pi[-1]

    def run_bootstrap(self):
        lg.info("Running %d bootstrap replicates" % self.opts.bootstrap)

        if self.opts.ncpu <= 1:
            bs_results = []
            for b in xrange(self.opts.bootstrap):
                bs_results.append(
                    bootstrap_replicate(self.raw_scores, self.opts, b + 1))
        else:
            pool = Pool(processes=self.opts.ncpu)
            map_result = pool.map_async(
                partial(bootstrap_replicate, self.raw_scores, self.opts),
                (b + 1 for b in xrange(self.opts.bootstrap))
            )
            bs_results = map_result.get()
            pool.close()
        self.bs_counts = np.vstack([_[0] for _ in bs_results])
        self.bs_props = np.vstack([_[1] for _ in bs_results])

    def output_report(self, tl, bootdata=None):
        _flens = self.annotation.feature_length()
        _cols = OrderedDict()
        _cols['transcript']        = self.featnames
        _cols['transcript_length'] = [_flens[f] for f in self.featnames]
        _cols['final_counts']      = self.final_counts
        _cols['final_conf']        = self.final_conf
        _cols['final_prop']        = self.final_prop

        _cols['unique_count']      = tl.reassign('unique').sum(0).A1
        _cols['init_aligned']      = tl.z[0].ceil().sum(0).A1
        _cols['init_best']         = tl.reassign('exclude', iteration=0).sum(0).A1
        _cols['init_best_random']  = tl.reassign('choose', iteration=0).sum(0).A1
        _cols['init_best_avg']     = tl.reassign('average', iteration=0).sum(0).A1
        _cols['init_prop']         = tl.pi[1]

        if self.bs_counts is not None:
            ci_level = self.opts.bootstrap_ci
            boot_cols = bootstrap_stats(self.bs_counts, ci_level)
            _cols['bootFC_mean'] = boot_cols[0]
            _cols['bootFC_median'] = boot_cols[1]
            _cols['bootFC_min'] = boot_cols[2]
            _cols['bootFC_max'] = boot_cols[3]
            _cols['bootFC_lowerCI'] = boot_cols[4]
            _cols['bootFC_upperCI'] = boot_cols[5]
            _cols['bootFC_pvalue'] = boot_cols[6]

        _table = [[_cols[k][i] for k in _cols] for i in xrange(self.shape[1])]
        _table.sort(key=lambda x: -x[4])
        _table.sort(key=lambda x: -x[2])

        # Create comment line
        _comment = ['## Run Info']
        _comment += ['%s:%s' % (k,v) for k,v in self.run_info.iteritems()]
        _report = [ _comment, _cols.keys() ]
        _report += [[str(_) for _ in l] for l in _table]
        return _report

    def output_matrix(self, tl, iteration=-1):
        rows = [ ['readid'] + self.featnames ]
        for i, rn in enumerate(self.readnames):
            rows.append([rn] + [str(v) for v
                                in tl.z[iteration].getrow(i).toarray()[0, :]]
                        )
        return rows

    def update_sam(self, tl, outsam_name):
        if self.opts.reassign_mode == 'conf':
            mat = tl.reassign(self.opts.reassign_mode, self.opts.conf_prob)
        else:
            mat = tl.reassign(self.opts.reassign_mode)
        best_feats = {i: self.featnames[j] for i, j in zip(*mat.nonzero())}

        readkeys_lookup = {}
        for ridx, d in self.readkeys.iteritems():
            for fidx, rkset in d.iteritems():
                for pair in rkset:
                    readkeys_lookup[pair[0]] = (ridx, fidx)
                    if len(pair) == 2:
                        readkeys_lookup[pair[1]] = (ridx, fidx)

        sf = pysam.AlignmentFile(self.opts.samfile)
        outsam = pysam.AlignmentFile(outsam_name, 'wb', template=sf)
        for r in sf:
            ridx, fidx = readkeys_lookup.get(bparse.readkey(r), (None, None))
            if ridx is None:
                r.set_tag('YC', c2str((248, 248, 248)))
                r.mapping_quality = 0
            else:
                prob = tl.z[-1][ridx, fidx]
                r.set_tag('XP', int(round(prob * 100)))
                r.set_tag('XT', self.featnames[fidx])
                r.set_tag('ZT', best_feats[ridx])
                if mat[ridx, fidx] > 0: # Best
                    r.set_tag('YC', c2str(DARK2_PALETTE['vermilion']))
                elif prob >= 0.2:
                    r.set_tag('YC', c2str(DARK2_PALETTE['yellow']))
                else:
                    r.set_tag('YC', c2str(GREENS[2]))
            _ = outsam.write(r)

    def __str__(self):
        _ret = '<Telescope samfile=%s, gtffile=%s>' % (self.opts.samfile, self.opts.gtffile)
        return _ret


class TelescopeLikelihood(object):
    """

    """
    def __init__(self, score_matrix, opts):
        """
        """
        self.epsilon = opts.em_epsilon
        self.max_iter = opts.max_iter

        # N fragments x K transcripts
        self.N, self.K = score_matrix.shape

        # pi[j] is the proportion of fragments that originate from
        # transcript j. Initial value assumes that all transcripts contribute
        # equal proportions of fragments
        self.pi = [ np.repeat(1./self.K, self.K) ]

        # theta[j] is the proportion of non-unique fragments that need to be
        # reassigned to transcript j. Initial value assumes that all transcripts
        # are reassigned an equal proportion of fragments
        self.theta = [ np.repeat(1./self.K, self.K) ]


        # Q[i,] is the set of mapping qualities for fragment i, where Q[i,j]
        # represents the conditional probability of observing fragment i assuming
        # it was generated from transcript j. We calculate this by scaling the
        # raw alignment score by the maximum alignment score observed.
        max_score = score_matrix.max()
        self.Q = score_matrix.multiply(100. / max_score ).expm1()

        # z[i,] is the partial assignment weights for fragment i, where z[i,j]
        # is the expected value for fragment i originating from transcript j. The
        # initial estimate is the normalized mapping qualities:
        # z_init[i,] = Q[i,] / sum(Q[i,])
        self.z = [ self.Q.normr() ]

        # Y[i] is the ambiguity indicator for fragment i, where Y[i]=1 if
        # fragment i is aligned to multiple transcripts and Y[i]=0 otherwise.
        # Store as N x 1 matrix
        self.Y = np.where(self.Q.countr() != 1, 1, 0)[:,None]

        # Log-likelihood score
        self.lnl = []

        # Prior values
        self.pi_prior = opts.pi_prior
        self.theta_prior = opts.theta_prior

        # Precalculated values
        self._weights = self.Q.maxr()             # Weight assigned to each fragment
        self._total_wt = self._weights.sum()      # Total weight
        self._ambig_wt = self._weights.multiply(self.Y).sum() # Weight of ambig frags
        self._unique_wt = self._weights.multiply(1-self.Y).sum()

        # Weighted prior values
        self._pi_prior_wt = self.pi_prior * self._weights.max()
        self._theta_prior_wt = self.theta_prior * self._weights.max()
        #
        self._pisum0 = self.Q.multiply(1-self.Y).sumc()

    def estep(self):
        """ Calculate the expected values of z
                E(z[i,j]) = ( pi[j] * theta[j]**Y[i] * Q[i,j] ) /
        """
        # assert len(self.z) == len(self.pi) == len(self.theta)
        _pi = self.pi[-1]
        _theta = self.theta[-1]
        _numerator = self.Q.multiply(csr_matrix(_pi * (_theta**self.Y)))
        self.z.append(_numerator.normr())

    def mstep(self):
        """ Calculate the maximum a posteriori (MAP) estimates for pi and theta

        """
        # assert (len(self.z)-1) == len(self.pi) == len(self.theta)
        # The expected values of z weighted by mapping score
        _weighted = self.z[-1].multiply(self._weights)

        # Estimate theta_hat
        _thetasum = _weighted.multiply(self.Y).sumc()
        _theta_denom = self._ambig_wt + self._theta_prior_wt * self.K
        _theta_hat = (_thetasum + self._theta_prior_wt) / _theta_denom

        # Estimate pi_hat
        _pisum = self._pisum0 + _thetasum
        _pi_denom = self._ambig_wt + self._unique_wt + self._pi_prior_wt * self.K
        _pi_hat = (_pisum + self._pi_prior_wt) / _pi_denom

        self.theta.append(_theta_hat.A1)
        self.pi.append(_pi_hat.A1)

    def calculate_lnl(self):
        _z, _p, _t = self.z[-1], self.pi[-1], self.theta[-1]
        cur = _z.multiply(self.Q.multiply(_p * _t**self.Y).log1p()).sum()
        self.lnl.append(cur)

    def em(self, use_likelihood=False, loglev=lg.WARNING):
        msg = 'Iteration %d, lnl=%g, diff=%g'
        converged = False
        while not converged:
            self.estep()
            self.mstep()
            self.calculate_lnl()

            iidx = len(self.lnl)
            if iidx <= 1:
                diff_lnl = float('inf')
            else:
                diff_lnl = abs(self.lnl[-1] - self.lnl[-2])
            diff_est = abs(self.pi[-1] - self.pi[-2]).sum()

            lg.log(loglev, msg % (iidx, self.lnl[-1], diff_est))
            if use_likelihood:
                converged = diff_lnl < self.epsilon or iidx >= self.max_iter
            else:
                converged = diff_est < self.epsilon or iidx >= self.max_iter

    def reassign(self, method, thresh=0.9, iteration=-1):
        """ Reassign fragments to expected transcripts

        Running EM finds the expected fragment assignment weights at the MAP
        estimates of pi and theta. This function reassigns all fragments based
        on these assignment weights. A simple heuristic is to assign each
        fragment to the transcript with the highest assignment weight.

        In practice, not all fragments have exactly one best hit. The "method"
        argument defines how we deal with fragments that are not fully resolved
        after EM:
                exclude - reads with > 1 best hits are excluded
                choose  - one of the best hits is randomly chosen
                average - read is evenly divided among best hits
                conf    - only confident reads are reassigned
                unique  - only uniquely aligned reads
        Args:
            method:
            thresh:
            iteration:

        Returns:
            matrix where m[i,j] == 1 iff read i is reassigned to transcript j

        """
        _z = self.z[iteration]
        if method == 'exclude':
            # Identify best hit(s), then exclude rows with >1 best hits
            v = _z.rowmax_identity()
            return v.multiply(np.where(v.sumr()>1, 0, 1))
        elif method == 'choose':
            # Identify best hit(s), then randomly choose reassignment
            v = _z.rowmax_identity()
            return v.choose_random()
        elif method == 'average':
            # Identify best hit(s), then divide by row sum
            v = _z.rowmax_identity()
            return v.normr()
        elif method == 'conf':
            # Zero out all values less than threshold
            # If thresh > 0.5 then at most
            v = _z.apply_func(lambda x: x if x >= thresh else 0)
            # Average each row so each sums to 1.
            return v.normr()
        elif method == 'unique':
            # Zero all rows that are ambiguous
            return _z.multiply(1 - self.Y).ceil()


def run_telescope_id(args):
    opts = IDOptions(args)
    configure_logging(opts)
    lg.info('\n{}\n'.format(opts))
    total_time = time()

    ts = Telescope(opts)

    # Load annotation file
    stime = time()
    ts.load_annotation()
    lg.info("Loaded annotation in %s" % fmtmins(time() - stime))

    # Load alignment file
    stime = time()
    ts.load_alignment()
    lg.info("Loaded alignment in %s" % fmtmins(time() - stime))
    #
    # lg.info("Processed {} fragments".format(ts.run_info['fragments']))
    # lg.info("\t{} fragments were unmapped".format(ts.run_info['unmapped']))
    # lg.info("\t{} fragments did not map to features".format(
    #     ts.run_info['nofeat']))
    # lg.info("\t{} fragments mapped uniquely".format(
    #     ts.run_info['unique']))
    # lg.info("\t{} fragments mapped ambiguously".format(
    #     ts.run_info['ambig']))
    #
    return


def _run_telescope_id(args):
    opts = IDOptions(args)
    configure_logging(opts)
    lg.info('\n{}\n'.format(opts))
    total_time = time()

    """ Load annotation """
    stime = time()
    flookup = Annotation(opts.gtffile, min_overlap=opts.overlap_threshold)
    lg.info("Loaded annotation in %s" % fmtmins(time() - stime))

    """ Load alignment """
    stime = time()
    samfile = pysam.AlignmentFile(opts.samfile)
    mapped, aln_counts = load_alignment(samfile, flookup, opts)
    lg.info("Loaded alignment in %s" % fmtmins(time() - stime))

    lg.info("Processed %d fragments" % sum(aln_counts.values()))
    lg.info("\t%d fragments were unmapped" % aln_counts['unmapped'])
    lg.info("\t%d fragments mapped to one or more positions on reference genome" % (
        aln_counts['unique'] + aln_counts['ambig'] + aln_counts['nofeat']))
    lg.info(
        "\t\t%d fragments mapped to reference but did not map to any transcripts" %
        aln_counts['nofeat'])
    lg.info("\t\t%d fragments have at least one alignment to a transcript" % (
    aln_counts['unique'] + aln_counts['ambig']))
    lg.info("\t\t\t%d fragments align uniquely to a single transcript" %
          aln_counts['unique'])
    lg.info("\t\t\t%d fragments align ambiguously to multiple transcripts" %
          aln_counts['ambig'])

    """ Create data structure """
    stime = time()
    ridx = {}
    gidx = {}
    d = []
    for rname,r in mapped.items():
        i = ridx.setdefault(rname,len(ridx))
        for gname,aln in r.feat_aln_map.items():
            j = gidx.setdefault(gname,len(gidx))
            d.append((i, j, aln.AS + aln.query_length))

    tm = TelescopeModel(ridx, gidx, data=d)
    lg.info("Created data structure in %s" % fmtmins(time() - stime))

    #
    # if opts.verbose:
    #     print("done.", file=sys.stderr)
    #     print("Time to create data structure:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # # Save some memory if you are not creating an updated SAM:
    # if not opts.updated_sam:
    #     if opts.verbose:
    #         print("Clearing reads from memory.", file=sys.stderr)
    #     mapped = None
    #     samfile.close()
    #
    # """ Initial checkpoint """
    # if opts.checkpoint:
    #     if opts.verbose:
    #         print("Checkpointing... ", end=' ', file=sys.stderr)
    #         substart = time()
    #
    #     with open(opts.generate_filename('checkpoint.init.p'),'wb') as outh:
    #         tm.dump(outh)
    #
    #     if opts.verbose:
    #         print("done.", file=sys.stderr)
    #         print("Time to write checkpoint:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # """ Initial output matrix """
    # if opts.out_matrix:
    #     if opts.verbose:
    #         print("Writing initial model matrix...", end=' ', file=sys.stderr)
    #         substart = time()
    #
    #     with open(opts.generate_filename('matrix_init.p'),'wb') as outh:
    #         tm.dump(outh)
    #
    #     if opts.verbose:
    #         print("done.", file=sys.stderr)
    #         print("Time to write initial matrix:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # """ Reassignment """
    # if opts.verbose:
    #     print("Reassiging reads:", file=sys.stderr)
    #     print("(Reads,Transcripts)=%dx%d" % (tm.shape), file=sys.stderr)
    #     print("Delta Change:", file=sys.stderr)
    #     substart = time()

    tm.matrix_em(opts)

    # if opts.verbose:
    #     print("Time for EM iteration:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # """ Checkpoint 2 """
    # if opts.checkpoint:
    #     if opts.verbose:
    #         print("Checkpointing... ", end=' ', file=sys.stderr)
    #         substart = time()
    #
    #     with open(opts.generate_filename('checkpoint.final.p'),'wb') as outh:
    #         tm.dump(outh)
    #
    #     if opts.verbose:
    #         print("done.", file=sys.stderr)
    #         print("Time to write checkpoint:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # """ Final output matrix """
    # if opts.out_matrix:
    #     if opts.verbose:
    #         print("Writing final model matrix...", end=' ', file=sys.stderr)
    #         substart = time()
    #
    #     with open(opts.generate_filename('matrix_final.p'), 'wb') as outh:
    #         tm.dump(outh)
    #
    #     if opts.verbose:
    #         print("done.", file=sys.stderr)
    #         print("Time to write final matrix:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # """ Generate report """
    # if opts.verbose:
    #     print("Generating report... ", end=' ', file=sys.stderr)
    #     substart = time()
    #
    # # report = tm.make_report(opts.conf_prob)
    # report = make_report(tm, aln_counts, flookup.feature_length(), opts)
    # # comment = ['%s:%d' % (k,v) for k,v in aln_counts.iteritems()] + ['transcripts:%d' % len(tm.txnames)]
    # with open(opts.generate_filename('telescope_report.tsv'),'w') as outh:
    #     # print >>outh, '## RunInfo\t%s' % '\t'.join(comment)
    #     for row in report:
    #         print('\t'.join(f for f in row), file=outh)
    #
    # lg.info()
    # if opts.verbose:
    #     print("done.", file=sys.stderr)
    #     print("Time to generate report:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # """ Update alignment """
    # if opts.updated_sam:
    #     if opts.verbose:
    #         print("Updating alignment...", end=' ', file=sys.stderr)
    #         substart = time()
    #
    #     updated_samfile = pysam.AlignmentFile(opts.generate_filename('updated.sam'), 'wh', header=samfile.header)
    #     update_alignment(tm, mapped, updated_samfile, opts)
    #     updated_samfile.close()
    #
    #     if opts.verbose:
    #         print("done.", file=sys.stderr)
    #         print("Time to update alignment:".ljust(40) + format_minutes(time() - substart), file=sys.stderr)
    #
    # samfile.close()
    return
