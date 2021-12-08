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

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"

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
    telescope_resume.BulkResumeOptions.add_arguments(resume_parser)
    resume_parser.set_defaults(func=telescope_resume.run)

    test_parser = subparsers.add_parser('test',
        description='''Print a test command''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    test_parser.set_defaults(func=generate_test_command)

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

    ''' Parser for cellsort '''
    cellsort_parser = subparsers.add_parser('cellsort',
        description='''Sort and filter BAM file according to cell barcode''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_cellsort.CmdOpts.add_arguments(cellsort_parser)
    cellsort_parser.set_defaults(func=stellarscope_cellsort.run)

    ''' Parser for cellsort '''
    merge_parser = subparsers.add_parser(
        'merge',
        description='''Merge telescope-generated transposable element counts with a gene count matrix''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_merge.CmdOpts.add_arguments(merge_parser)
    merge_parser.set_defaults(func=stellarscope_merge.run)

    args = parser.parse_args()
    args.func(args)


# if __name__ == '__main__':
#     main()
