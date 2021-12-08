#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Main functionality of Telescope

"""
from __future__ import absolute_import

import sys
import os
import argparse
import errno

from ._version import VERSION
from . import telescope_assign
from . import telescope_resume
from . import stellarscope_cellsort
from . import stellarscope_assign
from . import stellarscope_merge
from . import stellarscope_resume

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"

def generate_test_command(args, singlecell = False):
    try:
        _ = FileNotFoundError()
    except NameError:
        class FileNotFoundError(OSError):
            pass

    _base = os.path.dirname(os.path.abspath(__file__))
    _data_path = os.path.join(_base, 'data')
    _alnpath = os.path.join(_data_path, 'alignment.bam')
    _gtfpath = os.path.join(_data_path, 'annotation.gtf')
    if not os.path.exists(_alnpath):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), _alnpath
        )
    if not os.path.exists(_gtfpath):
        raise FileNotFoundError(
            errno.ENOENT, os.strerror(errno.ENOENT), _gtfpath
        )
    if singlecell == False:
        print('telescope assign %s %s' % (_alnpath, _gtfpath), file=sys.stdout)
    else:
        print('stellarscope assign %s %s' % (_alnpath, _gtfpath), file=sys.stdout)

TS_USAGE = ''' %(prog)s <command> [<args>]

The most commonly used commands are:
    assign    Reassign ambiguous fragments that map to repetitive elements
    resume    Resume previous run from checkpoint file
    test      Generate a command line for testing
'''

def telescope():

    parser = argparse.ArgumentParser(
        description='telescope: Locus-specific quantification of transposable element expression from RNA-seq data',
        usage=TS_USAGE
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser.add_argument('--version',
        action='version',
        version=VERSION,
        default=VERSION,
    )

    subparsers = parser.add_subparsers(help='sub-command help')

    ''' Parser for assign '''
    assign_parser = subparsers.add_parser('assign',
        description='''Reassign ambiguous fragments that map to repetitive
                       elements''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_assign.TelescopeAssignOptions.add_arguments(assign_parser)
    assign_parser.set_defaults(func=telescope_assign.run)

    ''' Parser for resume '''
    resume_parser = subparsers.add_parser('resume',
        description='''Resume a previous telescope run''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_resume.TelescopeResumeOptions.add_arguments(resume_parser)
    resume_parser.set_defaults(func=telescope_resume.run)

    test_parser = subparsers.add_parser('test',
        description='''Print a test command''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    test_parser.set_defaults(func=lambda args: generate_test_command(args, singlecell = False))

    args = parser.parse_args()
    args.func(args)


ST_USAGE = ''' %(prog)s <command> [<args>]

The most commonly used commands are:
    cellsort  Sort and filter BAM file according to cell barcode    
    assign    Reassign ambiguous fragments that map to repetitive elements
'''

def stellarscope():
    parser = argparse.ArgumentParser(
        description='''
            stellarscope: Locus-specific quantification of transposable element 
            expression in single-cell RNA-seq data''',
        usage=ST_USAGE
    )

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    parser.add_argument('--version',
        action='version',
        version=VERSION,
        default=VERSION,
    )

    subparsers = parser.add_subparsers(help='sub-command help')

    ''' Parser for assign '''
    assign_parser = subparsers.add_parser(
        'assign',
        description='''Reassign ambiguous fragments that map to repetitive elements (scRNA-seq)''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_assign.StellarscopeAssignOptions.add_arguments(assign_parser)
    assign_parser.set_defaults(func=stellarscope_assign.run)

    ''' Parser for resume '''
    resume_parser = subparsers.add_parser(
        'resume',
        description='''Resume a previous stellarscope run''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_resume.StellarscopeResumeOptions.add_arguments(resume_parser)
    resume_parser.set_defaults(func=stellarscope_resume.run)

    ''' Parser for cellsort '''
    cellsort_parser = subparsers.add_parser('cellsort',
        description='''Sort and filter BAM file according to cell barcode''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_cellsort.StellarscopeCellSortOptions.add_arguments(cellsort_parser)
    cellsort_parser.set_defaults(func=stellarscope_cellsort.run)

    ''' Parser for merge '''
    merge_parser = subparsers.add_parser(
        'merge',
        description='''Merge stellarscope-generated single-cell transposable element counts with a 
        single-cell gene count matrix''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_merge.StellarscopeMergeOptions.add_arguments(merge_parser)
    merge_parser.set_defaults(func=stellarscope_merge.run)

    ''' Parser for test '''
    test_parser = subparsers.add_parser('test',
                                        description='''Print a test command''',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        )
    test_parser.set_defaults(func=lambda args: generate_test_command(args, singlecell=True))

    args = parser.parse_args()
    args.func(args)


# if __name__ == '__main__':
#     main()
