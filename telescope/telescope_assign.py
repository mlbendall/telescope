# -*- coding: utf-8 -*-
""" Telescope assign

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


class TelescopeAssignOptions(utils.OptionsBase):

    OPTS = pkgutil.get_data('telescope', 'cmdopts/telescope_assign.yaml')

    def __init__(self, args):
        super().__init__(args)


def run(args):
    """

    Args:
        args:

    Returns:

    """
    option_class = TelescopeAssignOptions
    opts = option_class(args)
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
    ts.output_report(ts_model,
                     opts.outfile_path('run_stats.tsv'),
                     opts.outfile_path('TE_counts.tsv'))

    if opts.updated_sam:
        lg.info("Creating updated SAM file...")
        ts.update_sam(ts_model, opts.outfile_path('updated.bam'))

    lg.info("telescope assign complete (%s)" % fmtmins(time() - total_time))
    return
