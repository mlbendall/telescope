#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Main functionality of Telescope

"""
from __future__ import absolute_import

import sys
import argparse

from ._version import VERSION
from . import telescope_assign
from . import telescope_resume


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


USAGE   = ''' %(prog)s <command> [<args>]

The most commonly used commands are:
   assign    Reassign ambiguous fragments that map to repetitive elements
   resume    Resume previous run from checkpoint file
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

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
