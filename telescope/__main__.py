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
# from . import telescope_resume

from . import stellarscope_cellsort


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"


def generate_test_command(args, seq_mode):
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
    print('telescope %s assign %s %s' % (seq_mode, _alnpath, _gtfpath), file=sys.stdout)

def main():
    if len(sys.argv) == 1:
        empty_parser = argparse.ArgumentParser(
            description='Tools for analysis of repetitive DNA elements',
            usage=USAGE,
        )
        empty_parser.print_help(sys.stderr)
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description='Tools for analysis of repetitive DNA elements',
    )
    parser.add_argument('--version',
        action='version',
        version=VERSION,
        default=VERSION,
    )

    subparsers = parser.add_subparsers(help='Sequencing modality help', dest = 'sc_or_bulk')

    ''' Parser for scRNA-seq '''
    sc_parser = subparsers.add_parser('sc',
        description='''Telescope for single-cell RNA-sequencing data sets''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ''' Parser for bulk RNA-seq '''
    bulk_parser = subparsers.add_parser('bulk',
        description='''Telescope for bulk RNA-sequencing data sets''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    sc_subparser = sc_parser.add_subparsers(help='scRNA-seq sub-command help', dest = 'subcommand')

    ''' Parser for scRNA-seq assign '''
    sc_assign_parser = sc_subparser.add_parser('assign',
        description='''Reassign ambiguous fragments that map to repetitive
                       elements (scRNA-seq)''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_assign.scIDOptions.add_arguments(sc_assign_parser)
    sc_assign_parser.set_defaults(func=lambda args: telescope_assign.run(args, sc = True))

    ''' Parser for scRNA-seq resume '''
    sc_resume_parser = sc_subparser.add_parser('resume',
        description='''Resume a previous telescope run (scRNA-seq)''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_resume.scResumeOptions.add_arguments(sc_resume_parser)
    sc_resume_parser.set_defaults(func=lambda args: telescope_resume.run(args, sc = True))

    sc_test_parser = sc_subparser.add_parser('test',
        description='''Print a test command''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sc_test_parser.set_defaults(func=lambda args: generate_test_command(args, 'sc'))

    bulk_subparser = bulk_parser.add_subparsers(help='Bulk RNA-seq sub-command help', dest='subcommand')

    ''' Parser for bulk RNA-seq assign '''
    bulk_assign_parser = bulk_subparser.add_parser('assign',
        description='''Reassign ambiguous fragments that map to repetitive elements (bulk RNA-seq)''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_assign.BulkIDOptions.add_arguments(bulk_assign_parser)
    bulk_assign_parser.set_defaults(func=lambda args: telescope_assign.run(args, sc = False))

    ''' Parser for bulk RNA-seq resume '''
    bulk_resume_parser = bulk_subparser.add_parser('resume',
        description='''Resume a previous telescope run (scRNA-seq)''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_resume.BulkResumeOptions.add_arguments(bulk_resume_parser)
    bulk_resume_parser.set_defaults(func=lambda args: telescope_resume.run(args, sc = False))

    bulk_test_parser = bulk_subparser.add_parser('test',
        description='''Print a test command''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    bulk_test_parser.set_defaults(func=lambda args: generate_test_command(args, 'bulk'))

    args = parser.parse_args()
    args.func(args)


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
    # resume_parser = subparsers.add_parser('resume',
    #     description='''Resume a previous telescope run''',
    #     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    # )
    # telescope_resume.BulkResumeOptions.add_arguments(resume_parser)
    # resume_parser.set_defaults(func=telescope_resume.run)

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

    ''' Parser for cellsort '''
    cellsort_parser = subparsers.add_parser('cellsort',
        description='''Sort and filter BAM file according to cell barcode''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    stellarscope_cellsort.CmdOpts.add_arguments(cellsort_parser)
    cellsort_parser.set_defaults(func=stellarscope_cellsort.run)




    args = parser.parse_args()
    args.func(args)


# if __name__ == '__main__':
#     main()
