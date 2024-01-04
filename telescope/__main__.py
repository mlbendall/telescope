#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Main functionality of Telescope

"""
import sys
import os
import argparse
import errno

from telescope import __version__
from . import telescope_assign
from . import telescope_resume


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2019 Matthew L. Bendall"


USAGE   = ''' %(prog)s <command> [<args>]

The most commonly used commands are:
   assign    Reassign ambiguous fragments that map to repetitive elements
   resume    Resume previous run from checkpoint file
   test      Generate a command line for testing
'''

def generate_test_command(args):
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
    print('telescope assign %s %s' % (_alnpath, _gtfpath), file=sys.stdout)

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
        version=__version__,
        default=__version__,
    )

    subparser = parser.add_subparsers(help='Sub-command help', dest='subcommand')

    ''' Parser for bulk RNA-seq assign '''
    assign_parser = subparser.add_parser('assign',
        description='''Reassign ambiguous fragments that map to repetitive elements''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_assign.AssignOptions.add_arguments(assign_parser)
    assign_parser.set_defaults(func=lambda args: telescope_assign.run(args))

    ''' Parser for bulk RNA-seq resume '''
    resume_parser = subparser.add_parser('resume',
        description='''Resume a previous telescope run''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    telescope_resume.ResumeOptions.add_arguments(resume_parser)
    resume_parser.set_defaults(func=lambda args: telescope_resume.run(args))

    test_parser = subparser.add_parser('test',
        description='''Print a test command''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    test_parser.set_defaults(func=lambda args: generate_test_command(args))

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
