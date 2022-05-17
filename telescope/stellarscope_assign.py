# -*- coding: utf-8 -*-
""" Stellarscope assign

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
import pkgutil

import numpy as np

from . import utils
from .utils.helpers import format_minutes as fmtmins
from .utils.model import TelescopeLikelihood
from .utils.model_stellarscope import Stellarscope
from .utils.annotation import get_annotation_class
from .utils.sparse_plus import csr_matrix_plus as csr_matrix

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"

def fit_telescope_model(ts: Stellarscope, pooling_mode: str) -> TelescopeLikelihood:

    if pooling_mode == 'individual':
        ''' Initialise the z matrix for all reads '''
        z = ts.raw_scores.copy()
        for barcode in ts.barcodes:
            if barcode in ts.barcode_read_indices:
                _rows = ts.barcode_read_indices[barcode]
                ''' Create likelihood object using only reads from the cell '''
                _cell_raw_scores = csr_matrix(ts.raw_scores[_rows,:].copy())
                ts_model = TelescopeLikelihood(_cell_raw_scores, ts.opts)
                ''' Run EM '''
                ts_model.em(use_likelihood=ts.opts.use_likelihood, loglev=lg.DEBUG)
                ''' Add estimated posterior probs to the final z matrix '''
                z[_rows, :] = ts_model.z
        ts_model = TelescopeLikelihood(ts.raw_scores, ts.opts)
        ts_model.z = z
    elif pooling_mode == 'pseudobulk':
        ''' Create likelihood '''
        ts_model = TelescopeLikelihood(ts.raw_scores, ts.opts)
        ''' Run Expectation-Maximization '''
        ts_model.em(use_likelihood=ts.opts.use_likelihood, loglev=lg.INFO)
    elif pooling_mode == 'celltype':
        z = ts.raw_scores.copy()
        for celltype, df in ts.barcode_celltypes.groupby('celltype'):
            celltype_barcodes = set(df['barcode']).intersection(ts.barcodes)
            if celltype_barcodes:
                _rows = np.unique(np.concatenate([ts.barcode_read_indices[bc] for bc in celltype_barcodes]))
                ''' Create likelihood object using only reads from the celltype '''
                _celltype_raw_scores = csr_matrix(ts.raw_scores[_rows, :].copy())
                ts_model = TelescopeLikelihood(_celltype_raw_scores, ts.opts)
                ''' Run EM '''
                ts_model.em(use_likelihood=ts.opts.use_likelihood, loglev=lg.DEBUG)
                ''' Add estimated posterior probs to the final z matrix '''
                z[_rows, :] = ts_model.z
        ts_model = TelescopeLikelihood(ts.raw_scores, ts.opts)
        ts_model.z = z
    else:
        raise ValueError('Argument "pooling_mode" should be one of (individual, pseudobulk, celltype)')

    return ts_model

class StellarscopeAssignOptions(utils.OptionsBase):
    """
    import command options
    """
    OPTS = pkgutil.get_data('telescope', 'cmdopts/stellarscope_assign.yaml')

    def __init__(self, args):
        super().__init__(args)

        if self.pooling_mode == 'celltype':
            if self.celltypefile is None:
                raise ValueError('Pooling mode of "celltype" was specified but no cell type file provided.')
        elif self.celltypefile is not None:
            lg.info(f'Argument "celltypefile" provided but pooling_mode={self.pooling_mode}, '
                    f'not using provided cell type file')

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
    option_class = StellarscopeAssignOptions
    opts = option_class(args)
    utils.configure_logging(opts)
    lg.info('\n{}\n'.format(opts))
    total_time = time()

    ''' Create Telescope object '''
    ts = Stellarscope(opts)

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

    ''' Free up memory used by annotation '''
    annot = None
    lg.debug('garbage: {:d}'.format(gc.collect()))

    ''' Save object checkpoint '''
    ts.save(opts.outfile_path('checkpoint'))
    if opts.skip_em:
        lg.info("Skipping EM...")
        lg.info("stellarscope assign complete (%s)" % fmtmins(time()-total_time))
        return

    ''' Seed RNG '''
    seed = ts.get_random_seed()
    lg.debug("Random seed: {}".format(seed))
    np.random.seed(seed)

    lg.info('Running Expectation-Maximization...')
    stime = time()
    ts_model = fit_telescope_model(ts, opts.pooling_mode)
    lg.info("EM completed in %s" % fmtmins(time() - stime))

    # Output final report
    lg.info("Generating Report...")
    ts.output_report(ts_model,
                     opts.outfile_path('run_stats.tsv'),
                     opts.outfile_path('TE_counts.mtx'),
                     opts.outfile_path('barcodes.tsv'),
                     opts.outfile_path('features.tsv'))

    if opts.updated_sam:
        lg.info("Creating updated SAM file...")
        ts.update_sam(ts_model, opts.outfile_path('updated.bam'))

    lg.info("stellarscope assign complete (%s)" % fmtmins(time() - total_time))
    return
