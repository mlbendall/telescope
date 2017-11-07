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
from .telescope_assign import IDOptions

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

class ResumeOptions(IDOptions):
    # This class is exactly the same as IDOptions for now
    pass


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

    lg.info("telescope resume complete (%s)" % fmtmins(time() - total_time))

    return
