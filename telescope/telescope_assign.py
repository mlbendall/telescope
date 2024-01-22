# -*- coding: utf-8 -*-
""" Telescope id

"""
from __future__ import print_function
from __future__ import absolute_import
from builtins import super

import sys
import os
from time import time
import logging as lg
import gc
import tempfile
import atexit
import shutil

import numpy as np

from . import utils
from .utils.helpers import format_minutes as fmtmins
from .utils.model import Telescope, TelescopeLikelihood
from .utils.annotation import get_annotation_class

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


class AssignOptions(utils.SubcommandOptions):
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
            - no_feature_key:
                default: __no_feature
                help: Used internally to represent alignments. Must be different
                      from all other feature names.
            - ncpu:
                default: 1
                type: int
                help: Number of cores to use. (Multiple cores not supported yet).
            - tempdir:
                help: Path to temporary directory. Temporary files will be stored
                      here. Default uses python tempfile package to create the
                      temporary directory.
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
                      statistics report by default. This argument determines what
                      mode will be used for the outputted counts file.
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
            - annotation_class:
                default: intervaltree
                choices:
                    - intervaltree
                    - htseq
                help: Annotation class to use for finding overlaps. Both htseq and
                      intervaltree appear to yield identical results. Performance
                      differences are TBD.
            - stranded_mode:
                type: str
                default: None
                choices:
                    - None
                    - RF
                    - R
                    - FR
                    - F
                help: Options for considering feature strand when assigning reads. 
                      If None, for each feature in the annotation, returns counts for the positive strand and negative strand. 
                      If not None, specifies the orientation of paired end reads (RF - read 1 reverse strand, read 2 forward strand) and
                      single end reads (F - forward strand). 
        - Model Parameters:
            - pi_prior:
                type: int
                default: 0
                help: Prior on π. Equivalent to adding n unique reads.
            - theta_prior:
                type: int
                default: 200000
                help: >
                      Prior on θ. Equivalent to adding n non-unique reads. NOTE: It
                      is recommended to set this prior to a large value. This
                      increases the penalty for non-unique reads and improves
                      accuracy.
            - em_epsilon:
                type: float
                default: 1e-7
                help: EM Algorithm Epsilon cutoff
            - max_iter:
                type: int
                default: 100
                help: EM Algorithm maximum iterations
            - use_likelihood:
                action: store_true
                help: Use difference in log-likelihood as convergence criteria.
            - skip_em:
                action: store_true
                help: Exits after loading alignment and saving checkpoint file.
        """

    old_opts = """
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
        - out_matrix:
            action: store_true
            help: Output alignment matrix
        """

    def __init__(self, args):
        super().__init__(args)
        if self.logfile is None:
            self.logfile = sys.stderr

        if hasattr(self, 'tempdir') and self.tempdir is None:
            if hasattr(self, 'ncpu') and self.ncpu > 1:
                self.tempdir = tempfile.mkdtemp()
                atexit.register(shutil.rmtree, self.tempdir)

    def outfile_path(self, suffix):
        basename = '%s-%s' % (self.exp_tag, suffix)
        return os.path.join(self.outdir, basename)


def run(args):
    """

    Args:
        args:

    Returns:

    """
    opts = AssignOptions(args)
    utils.configure_logging(opts)
    lg.info('\n{}\n'.format(opts))
    total_time = time()

    ''' Create Telescope object '''
    ts = Telescope(opts)

    ''' Load annotation '''
    Annotation = get_annotation_class(opts.annotation_class)
    lg.info('Loading annotation...')
    stime = time()
    annot = Annotation(opts.gtffile, opts.attribute, opts.stranded_mode)
    lg.info("Loaded annotation in {}".format(fmtmins(time() - stime)))
    lg.info('Loaded {} features.'.format(len(annot.loci)))

    # annot.save(opts.outfile_path('test_annotation.p'))

    ''' Load alignments '''
    lg.info('Loading alignments...')
    stime = time()
    ts.load_alignment(annot)
    lg.info("Loaded alignment in {}".format(fmtmins(time() - stime)))

    ''' Print alignment summary '''
    ts.print_summary(lg.INFO)
    # if opts.ncpu > 1:
    #     sys.exit('not implemented yet')

    ''' Exit if no overlap '''
    if ts.run_info['overlap_unique'] + ts.run_info['overlap_ambig'] == 0:
        lg.info("No alignments overlapping annotation")
        lg.info("telescope assign complete (%s)" % fmtmins(time() - total_time))
        return

    ''' Free up memory used by annotation '''
    annot = None
    lg.debug('garbage: {:d}'.format(gc.collect()))

    ''' Save object checkpoint '''
    ts.save(opts.outfile_path('checkpoint'))
    if opts.skip_em:
        lg.info("Skipping EM...")
        lg.info("telescope assign complete (%s)" % fmtmins(time()-total_time))
        return

    ''' Seed RNG '''
    seed = ts.get_random_seed()
    lg.debug("Random seed: {}".format(seed))
    np.random.seed(seed)

    ''' Create likelihood '''
    ts_model = TelescopeLikelihood(ts.raw_scores, opts)

    ''' Run Expectation-Maximization '''
    lg.info('Running Expectation-Maximization...')
    stime = time()
    ts_model.em(use_likelihood=opts.use_likelihood, loglev=lg.INFO)
    lg.info("EM completed in %s" % fmtmins(time() - stime))

    # Output final report
    lg.info("Generating Report...")
    ts.output_report(ts_model, opts.outfile_path('run_stats.tsv'), opts.outfile_path('TE_counts.tsv'))

    if opts.updated_sam:
        lg.info("Creating updated SAM file...")
        ts.update_sam(ts_model, opts.outfile_path('updated.bam'))

    lg.info("telescope assign complete (%s)" % fmtmins(time() - total_time))
    return
