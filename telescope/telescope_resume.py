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
import gc

from . import utils
from .utils.helpers import format_minutes as fmtmins

from .utils.model import Telescope, TelescopeLikelihood

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"



class ResumeOptions(utils.SubcommandOptions):
    OPTS = """
    - Input Options:
        - samfile:
            positional: True
            help: Path to intermediate alignment file.
        - gtffile:
            positional: True
            help: Path to annotation file (GTF format)
        - attribute:
            default: locus
            help: GTF attribute that defines a transposable element locus. GTF
                  features that share the same value for --attribute will be
                  considered as part of the same locus.
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

    def outfile_path(self, suffix):
        basename = '%s-%s' % (self.exp_tag, suffix)
        return os.path.join(self.outdir, basename)

    def tempfile_path(self, suffix):
        basename = suffix.format(os.getpid())
        return os.path.join(self.outdir, basename)


def run(args):
    """

    Args:
        args:

    Returns:

    """
    opts = ResumeOptions(args)
    utils.configure_logging(opts)
    lg.info('\n{}\n'.format(opts))
    total_time = time()

    ''' Create Telescope object '''
    ts = Telescope(opts)

    ''' Load annotation '''
    lg.info('Loading annotation...')
    stime = time()
    ts.load_annotation()
    lg.info("Loaded annotation in {}".format(fmtmins(time() - stime)))
    lg.info('Loaded {} features.'.format(len(ts.annotation.loci)))
    ts.annotation = None
    lg.debug('GC collected {:d} items.'.format(gc.collect()))

    ''' Reload alignments '''
    lg.info('Loading alignments from "%s"...' % ts.tmp_bam)
    stime = time()
    mappings = ts.load_mappings(ts.tmp_bam)
    lg.info("Loaded alignment in {}".format(fmtmins(time() - stime)))

    lg.debug('calling:  Telescope._mapping_to_matrix')
    stime = time()
    ts._mapping_to_matrix(mappings)
    lg.debug('complete: Telescope._mapping_to_matrix (in {})'.format(fmtmins(time() - stime)))

    # ''' Create likelihood '''
    lg.debug('calling:  TelescopeLikelihood.__init__')
    ts_model = TelescopeLikelihood(ts.raw_scores, opts)
    lg.debug('complete:  TelescopeLikelihood.__init__ (in {})'.format(fmtmins(time() - stime)))
    #
    # ''' Run Expectation-Maximization '''
    lg.info('Running EM...')
    stime = time()
    ts_model.em(loglev=lg.INFO)
    lg.info("EM converged in %s" % fmtmins(time() - stime))
    #
    # # Output final report
    lg.info("Generating Report...")
    report_out = opts.outfile_path('telescope_report.tsv')
    if os.path.exists(report_out):
        report_out = opts.outfile_path('telescope_report.resume.tsv')

    ts.output_report(ts_model, report_out)

    if opts.updated_sam:
        lg.info("Creating updated SAM file...")
        sam_out = opts.outfile_path('updated.bam')
        if os.path.exists(sam_out):
            sam_out = opts.outfile_path('updated.resume.bam')
        ts.update_sam(ts_model, sam_out)

    return
