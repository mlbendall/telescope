# -*- coding: utf-8 -*-
from __future__ import absolute_import

import yaml
import logging
from collections import OrderedDict
# These are needed for eval statements:
import sys
import argparse


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


class SubcommandOptions(object):
    OPTS = """
    - Input Options:
        - infile:
            positional: True
            help: Input file.
    - Output Options:
        - outfile:
            positional: True
            help: Output file.
    """

    def __init__(self, args):
        self.opt_names, self.opt_groups = self._parse_yaml_opts(self.OPTS)
        for k, v in vars(args).items():
            if k in self.opt_names:
                setattr(self, k, v)
            else:
                # print("extra option: {} = {}".format(k,v))
                setattr(self, k, v)

    @classmethod
    def add_arguments(cls, parser):
        opt_names, opt_groups = cls._parse_yaml_opts(cls.OPTS)
        for group_name, args in opt_groups.items():
            argparse_grp = parser.add_argument_group(group_name, '')
            for arg_name, arg_d in args.items():
                _d = dict(arg_d)
                if _d.pop('hide', False):
                    continue
                if _d.pop('positional', False):
                    _arg_name = arg_name
                elif len(arg_name) == 1:
                    _arg_name = '-{}'.format(arg_name)
                else:
                    _arg_name = '--{}'.format(arg_name)

                if 'type' in _d:
                    _d['type'] = eval(_d['type'])

                argparse_grp.add_argument(_arg_name, **_d)

    @staticmethod
    def _parse_yaml_opts(opts_yaml):
        _opt_names = []
        _opt_groups = OrderedDict()
        for grp in yaml.load(opts_yaml):
            grp_name, args = list(grp.items())[0]
            _opt_groups[grp_name] = OrderedDict()
            for arg in args:
                arg_name, d = list(arg.items())[0]
                _opt_groups[grp_name][arg_name] = d
                _opt_names.append(arg_name)
        return _opt_names, _opt_groups

    def __str__(self):
        ret = []
        if hasattr(self, 'version'):
            ret.append('{:34}{}'.format('Version:', self.version))
        for group_name, args in self.opt_groups.items():
            ret.append('{}'.format(group_name))
            for arg_name in args.keys():
                v = getattr(self, arg_name, "Not set")
                v  = getattr(v, 'name') if hasattr(v, 'name') else v
                ret.append('    {:30}{}'.format(arg_name + ':', v))
        return '\n'.join(ret)


def configure_logging(opts):
    """ Configure logging options

    Args:
        opts: SubcommandOptions object. Important attributes are "quiet",
              "debug", and "logfile"
    Returns:  None
    """
    loglev = logging.INFO
    if opts.quiet:
        loglev = logging.WARNING
    if opts.debug:
        loglev = logging.DEBUG

    logfmt = '%(asctime)s %(levelname)-8s %(message)-60s'
    logfmt += ' (from %(funcName)s in %(filename)s:%(lineno)d)'
    logging.basicConfig(level=loglev,
                        format=logfmt,
                        datefmt='%Y-%m-%d %H:%M:%S',
                        stream=opts.logfile)
    return

# A very large integer
BIG_INT = 2**32 - 1