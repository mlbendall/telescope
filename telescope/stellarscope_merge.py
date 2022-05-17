# -*- coding: utf-8 -*-
""" Stellarscope merge

"""

from __future__ import print_function
from __future__ import absolute_import
from builtins import super

import os
import logging as lg
import pkgutil
from time import time
import scipy
from scipy import io
import pandas as pd

from . import utils
from .utils.helpers import format_minutes as fmtmins

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"

class StellarscopeMergeOptions(utils.OptionsBase):

    # import command options from the yaml file
    OPTS = pkgutil.get_data('telescope', 'cmdopts/stellarscope_merge.yaml')

    def __init__(self, args):
        super().__init__(args)


def run(args):
    """
    args: command line arguments for stellarscope merge
    """
    opts = StellarscopeMergeOptions(args)
    utils.configure_logging(opts)
    lg.debug('\n{}\n'.format(opts))
    total_time = time()

    # import the relevant data files
    lg.info('Loading gene counts...')
    gene_counts = scipy.sparse.csr_matrix(io.mmread(opts.gene_counts))
    gene_features = pd.read_csv(opts.gene_features, sep='\t', header=None)
    gene_barcodes = pd.read_csv(opts.gene_barcodes, sep='\t')

    lg.info('Loading transposable element counts...')
    TE_counts = scipy.sparse.csr_matrix(io.mmread(opts.TE_counts))
    TE_features = pd.read_csv(opts.TE_features, sep='\t', header=None)
    TE_barcodes = pd.read_csv(opts.TE_barcodes, sep='\t')

    if len(TE_barcodes) != len(gene_barcodes) or not TE_barcodes.iloc[:, 0].isin(gene_barcodes.iloc[:, 0]).all():
        lg.warning('Barcode sets do not match, only keeping barcodes present in both sets.')

    # align barcode sets to each other, only keeping barcodes present in both sets
    gene_bc_aligned, TE_bc_aligned = gene_barcodes.align(TE_barcodes, join='inner', axis=0)
    gene_count_rows, TE_count_rows = gene_bc_aligned.index.to_numpy(), TE_bc_aligned.index.to_numpy()

    # use common barcodes to combine count matrices
    lg.info('Combining count matrices...')
    merged_mtx = scipy.sparse.hstack(
        [gene_counts[gene_count_rows,:], TE_counts[TE_count_rows,:]]
    )

    # if the gene features data frame is bigger than the TE features, adjust the size of the TE features
    if TE_features.shape[1] < gene_features.shape[1]:
        TE_features = pd.concat([TE_features.iloc[:, 0]] * gene_features.shape[1],
                                axis=1, ignore_index=True)

    merged_features = gene_features.append(TE_features)
    merged_barcodes = pd.Series(gene_bc_aligned.iloc[:,0].values)


    # save files
    io.mmwrite(opts.outfile_path('merged_counts.mtx'), merged_mtx)
    merged_features.to_csv(opts.outfile_path('merged_features.tsv'), sep='\t', index=False, header=False)
    merged_barcodes.to_csv(opts.outfile_path('merged_barcodes.tsv'), sep='\t', index=False, header=False)

    lg.info("stellarscope merge complete (%s)" % fmtmins(time() - total_time))

    return

