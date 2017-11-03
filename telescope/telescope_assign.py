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
import logging

from . import utils
from .utils import format_minutes as fmtmins

from . import Telescope
from . import TelescopeLikelihood

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


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
    opts = IDOptions(args)
    utils.configure_logging(opts)
    logging.info('\n{}\n'.format(opts))
    total_time = time()

    ''' Create Telescope object '''
    ts = Telescope(opts)

    ''' Load annotation '''
    logging.info('Loading annotation...')
    stime = time()
    ts.load_annotation()
    logging.info("Loaded annotation in {}".format(fmtmins(time() - stime)))
    logging.info('Loaded {} features.'.format(len(ts.annotation.loci)))

    ''' Load alignments '''
    logging.info('Loading alignments...')
    stime = time()
    ts.load_alignment()
    logging.info("Loaded alignment in {}".format(fmtmins(time() - stime)))

    ''' Print alignment summary '''
    # _ainfo = ts.run_info['alignment_info']
    # input_pairs = _ainfo['map_2'] + _ainfo['unmap_2']
    # input_single = _ainfo['map_1'] + _ainfo['unmap_1']
    # total_frags = input_pairs + input_single
    _rinfo = ts.run_info
    logging.info("Alignment Summary:")
    logging.info('\t{} total fragments.'.format(_rinfo['total_fragments']))
    logging.info('\t\t{} mapped as pairs.'.format(_rinfo['mapped_pairs']))
    logging.info('\t\t{} mapped single.'.format(_rinfo['mapped_single']))
    logging.info('\t\t{} failed to map.'.format(_rinfo['unmapped']))
    logging.info('--')
    logging.info('\t{} fragments mapped to reference; of these'.format(
        _rinfo['mapped_pairs'] + _rinfo['mapped_single']))
    logging.info('\t\t{} had one unique alignment.'.format(_rinfo['unique']))
    logging.info('\t\t{} had multiple alignments.'.format(_rinfo['ambig']))
    logging.info('--')
    logging.info('\t{} fragments overlapped annotation; of these'.format(
        _rinfo['overlap_unique'] + _rinfo['overlap_ambig']))
    logging.info('\t\t{} had one unique alignment.'.format(_rinfo['overlap_unique']))
    logging.info('\t\t{} had multiple alignments.'.format(_rinfo['overlap_ambig']))
    logging.info('\n')

    ''' Create likelihood '''
    ts_model = TelescopeLikelihood(ts.raw_scores, opts)

    ''' Run Expectation-Maximization '''
    logging.info('Running EM...')
    stime = time()
    ts_model.em(loglev=logging.INFO)
    logging.info("EM converged in %s" % fmtmins(time() - stime))

    # Output final report
    logging.info("Generating Report...")
    ts.output_report(ts_model)

    if opts.updated_sam:
        logging.info("Creating updated SAM file...")
        ts.update_sam(ts_model)

    return
