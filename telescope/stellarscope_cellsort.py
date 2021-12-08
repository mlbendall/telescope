# -*- coding: utf-8 -*-
""" Stellarscope assign

"""

from __future__ import print_function
from __future__ import absolute_import
from builtins import super

import os
import logging as lg
import pkgutil
from time import time
import subprocess

from . import utils
from .utils.helpers import format_minutes as fmtmins


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"


class StellarscopeCellSortOptions(utils.OptionsBase):
    """

    """
    OPTS = pkgutil.get_data('telescope', 'cmdopts/stellarscope_cellsort.yaml')

    def __init__(self, args):
        super().__init__(args)

def run(args):
    """

    Args:
        args:

    Returns:

    """
    opts = StellarscopeCellSortOptions(args)
    utils.configure_logging(opts)
    lg.debug('\n{}\n'.format(opts))
    total_time = time()

    # Default arguments
    view_thread = 1
    sort_thread = opts.ncpu
    tempdir_arg = '' if opts.tempdir is None else '-T {}'.format(opts.tempdir)

    # Filter passing cell barcodes
    cmd1 = 'samtools view -@{ncpu:d} -u -F 4 -D CB:{bcfile:s} {inbam:s}'
    # Sort by cell barcode and read name
    cmd2 = 'samtools sort -@{ncpu:d} -n -t CB {tempdir_arg:s}'

    cmd = ' '.join([
        cmd1.format(ncpu=view_thread, bcfile=opts.barcodes, inbam=opts.infile),
        '|',
        cmd2.format(ncpu=sort_thread, tempdir_arg=tempdir_arg),
        '>',
        opts.outfile
    ])

    lg.info('Filtering and sorting by cell barcode')
    lg.debug('Command: {}'.format(cmd))
    output = subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell=True)
    lg.debug('Command output:\n{}\n'.format(output.decode()))
    lg.info("stellarscope cellsort complete (%s)" % fmtmins(time() - total_time))
    return

