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

import numpy as np

from . import utils
from .utils.helpers import format_minutes as fmtmins

from .utils.model import Telescope, TelescopeLikelihood
from .telescope_assign import IDOptions

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

class ResumeOptions(IDOptions):
    OPTS = """
    - Input Options:
        - checkpoint:
            positional: True
            help: Path to checkpoint file.
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
        - use_likelihood:
            action: store_true
            help: Use difference in log-likelihood as convergence criteria.
    """

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
    lg.info('Loading Telescope object from file...')
    ts = Telescope.load(opts.checkpoint)
    ts.opts = opts

    ''' Print alignment summary '''
    ts.print_summary(lg.INFO)

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

    ''' Output final report '''
    lg.info("Generating Report...")
    report_out = opts.outfile_path('telescope_report.tsv')
    ts.output_report(ts_model, report_out)

    # if opts.updated_sam:
    #     lg.info("Creating updated SAM file...")
    #     ts.update_sam(ts_model, opts.outfile_path('updated.bam'))

    lg.info("telescope resume complete (%s)" % fmtmins(time() - total_time))

    return
