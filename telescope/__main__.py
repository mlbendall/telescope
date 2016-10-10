#! /usr/bin/env python
# -*- coding: utf-8 -*-
""" Main functionality of Telescope

Example:
    $ python telescope/__main__.py --help
    $ python telescope/__main__.py id --help
    $ python telescope/__main__.py tag --help
    $ python telescope/__main__.py load --help  

"""

import sys
import argparse

from _version import VERSION
import telescope_id # from telescope_id import run_telescope_id
import telescope_tag # from telescope_tag import run_telescope_tag
import telescope_load # from telescope_load import run_telescope_load

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

# Set the usage string
USAGE   = ''' %(prog)s <command> [<args>]

The most commonly used commands are:
   id       Record changes to the repository
   tag      Add tags to an alignment
'''

def main():
    parser = argparse.ArgumentParser(description='Tools for analysis of repetitive DNA elements',
                                     usage=USAGE,
                                     )
    parser.add_argument('--version', action='version', version=VERSION, default=VERSION)

    subparsers = parser.add_subparsers(help='sub-command help')

    ''' Parser for ID '''
    id_parser = subparsers.add_parser('id',
                                      description='Reassign reads',
                                      formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     )

    inputopts = id_parser.add_argument_group('input', 'Input options')
    inputopts.add_argument('--ali_format', default='sam', help='Alignment Format. Only SAM is supported.')
    inputopts.add_argument('samfile', help='Path to alignment file')
    inputopts.add_argument('gtffile', help='Path to annotation file (GTF format)')
    inputopts.add_argument('--no_feature_key', default='__nofeature__',
                            help='Feature name for unassigned reads. Must not match any other feature name')

    outputopts = id_parser.add_argument_group('output', 'Output options')
    outputopts.add_argument('--verbose', action='store_true',
                            help='Prints verbose text while running')
    outputopts.add_argument('--outdir', default=".",
                             help='Output Directory')
    outputopts.add_argument('--exp_tag', default="telescope",
                            help='Experiment tag')
    #outputopts.add_argument('--min_final_guess', type=float, default=0.01,
    #                        help='Minimum final guess for genome to appear in report. Genomes with one or more final hits will always be included.')
    outputopts.add_argument('--out_matrix', action='store_true',
                            help='Output alignment matrix')
    outputopts.add_argument('--updated_sam', action='store_true', dest='updated_sam',
                            help='Generate an updated alignment file')
    outputopts.add_argument('--checkpoint', action='store_true', dest='checkpoint',
                            help='Enable checkpointing feature')
    outputopts.add_argument('--checkpoint_interval', type=int, default=10,
                            help='Number of EM iterations between checkpoints')
    outputopts.add_argument('--min_prob', type=float, default=0.2,
                            help='Minimum probability to be included in updated alignment file')
    outputopts.add_argument('--conf_prob', type=float, default=0.9,
                            help='Minimum probability for high confidence assignment')

    outputopts.add_argument('--reassign_mode', default='exclude', choices=['exclude','choose','random'],
                            help='Method for reassigning stubborn reads')

    modelopts = id_parser.add_argument_group('model', 'Model parameters')
    modelopts.add_argument('--piPrior', type=int, default=0,
                           help='Pi Prior equivalent to adding n unique reads')
    modelopts.add_argument('--thetaPrior', type=int, default=0,
                           help='Theta Prior equivalent to adding n non-unique reads')
    #modelopts.add_argument('--score_cutoff', type=float, default=0.01,
    #                       help='Minimum final probability score for alignment')
    modelopts.add_argument('--min_overlap', type=float, default=0.1,
                           help='Minimum fraction of read that must overlap to be assigned to feature.')


    emopts = id_parser.add_argument_group('em', 'EM parameters')
    emopts.add_argument('--emEpsilon', type=float, default=1e-7,
                        help='EM Algorithm Epsilon cutoff')
    emopts.add_argument('--maxIter', type=int, default=100,
                        help='EM Algorithm maximum iterations')

    id_parser.set_defaults(func=telescope_id.run_telescope_id)

    ''' Parser for TAG '''
    tag_parser = subparsers.add_parser('tag',
                                       description='Tag BAM file',
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                      )
    tag_parser.add_argument('--verbose', action='store_true',
                            help='Prints verbose text while running')
    tag_parser.add_argument('--gtffile', help='Path to annotation file (GTF format)')
    tag_parser.add_argument('samfile', nargs="?", default="-", help='Path to alignment file (default is STDIN)')
    tag_parser.add_argument('outfile', nargs="?", default="-", help='Output file (default is STDOUT)')
    tag_parser.add_argument('--min_overlap', type=float, default=0.1,
                           help='Minimum fraction of read that must overlap to be assigned to feature.')

    tag_parser.set_defaults(func=telescope_tag.run_telescope_tag)

    ''' Parser for LOAD '''
    load_parser = subparsers.add_parser('load',
                                       description='Load checkpoint file',
                                       formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                      )
    load_parser.add_argument('--verbose', action='store_true',
                             help='Prints verbose text while running')
    load_parser.add_argument('--outparam',
                             help='Output this parameter value')
    load_parser.add_argument('--prec', type=int, default=6,
                             help='Output precision')
    load_parser.add_argument('--float', action='store_true',
                             help='Force output as floats')
    load_parser.add_argument('--exp', action='store_true',
                             help='Force output as exponential')
    load_parser.add_argument('checkpoint', help='Checkpoint file')
    load_parser.add_argument('outfile', nargs="?", type=argparse.FileType('w'), default=sys.stdout,
                             help='Output file (default is STDOUT)')
    
    load_parser.set_defaults(func=telescope_load.run_telescope_load)

    args = parser.parse_args()
    args.func(args)

if __name__ == '__main__':
    main()
