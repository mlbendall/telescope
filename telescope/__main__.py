#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Main functionality of Telescope

Example:
    $ python telescope/__main__.py --help
    $ python telescope/__main__.py id --help
    $ python telescope/__main__.py tag --help
    $ python telescope/__main__.py load --help  

"""
from __future__ import absolute_import

import sys
import argparse

from ._version import VERSION
from . import telescope_assign
from . import telescope_resume

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

# Set the usage string
USAGE   = ''' %(prog)s <command> [<args>]

The most commonly used commands are:
   assign    Reassign ambiguous fragments that map to repetitive elements
'''

def main():
    parser = argparse.ArgumentParser(
        description='Tools for analysis of repetitive DNA elements',
        usage=USAGE,
    )
    parser.add_argument('--version',
        action='version',
        version=VERSION,
        default=VERSION,
    )

    subparsers = parser.add_subparsers(help='sub-command help')

    ''' Parser for ID '''
    # id_parser = subparsers.add_parser('id',
    #                                   description='Reassign reads',
    #                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    #                                  )
    # telescope_id.IDOptions.add_arguments(id_parser)
    # id_parser.set_defaults(func=telescope_id.run_telescope_id)

    ''' Parser for assign '''
    assign_parser = subparsers.add_parser('assign',
        description='''Reassign ambiguous fragments that map to repetitive
                       elements''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_assign.IDOptions.add_arguments(assign_parser)
    assign_parser.set_defaults(func=telescope_assign.run)

    ''' Parser for resume '''
    resume_parser = subparsers.add_parser('resume',
        description='''Resume a previous telescope run''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_resume.ResumeOptions.add_arguments(resume_parser)
    resume_parser.set_defaults(func=telescope_resume.run)


    # ''' Parser for TAG '''
    # tag_parser = subparsers.add_parser('tag',
    #                                    description='Tag BAM file',
    #                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    #                                   )
    # tag_parser.add_argument('--verbose', action='store_true',
    #                         help='Prints verbose text while running')
    # tag_parser.add_argument('--gtffile', help='Path to annotation file (GTF format)')
    # tag_parser.add_argument('samfile', nargs="?", default="-", help='Path to alignment file (default is STDIN)')
    # tag_parser.add_argument('outfile', nargs="?", default="-", help='Output file (default is STDOUT)')
    # tag_parser.add_argument('--min_overlap', type=float, default=0.1,
    #                        help='Minimum fraction of read that must overlap to be assigned to feature.')
    #
    # tag_parser.set_defaults(func=telescope_tag.run_telescope_tag)
    #
    # ''' Parser for LOAD '''
    # load_parser = subparsers.add_parser('load',
    #                                    description='Load checkpoint file',
    #                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    #                                   )
    # load_parser.add_argument('--verbose', action='store_true',
    #                          help='Prints verbose text while running')
    # load_parser.add_argument('--outparam',
    #                          help='Output this parameter value')
    # load_parser.add_argument('--prec', type=int, default=6,
    #                          help='Output precision')
    # load_parser.add_argument('--float', action='store_true',
    #                          help='Force output as floats')
    # load_parser.add_argument('--exp', action='store_true',
    #                          help='Force output as exponential')
    # load_parser.add_argument('checkpoint', help='Checkpoint file')
    # load_parser.add_argument('outfile', nargs="?", type=argparse.FileType('w'), default=sys.stdout,
    #                          help='Output file (default is STDOUT)')
    #
    # load_parser.set_defaults(func=telescope_load.run_telescope_load)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
