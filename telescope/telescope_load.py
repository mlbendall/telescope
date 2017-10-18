# -*- coding: utf-8 -*-
""" Telescope load

"""
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import object
import sys

from .utils.model import TelescopeModel

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


class LoadOpts(object):
    option_fields = ['verbose', 'checkpoint', 'outparam', 'prec', 'float', 'exp', 'outfile', ]

    def __init__(self, **kwargs):
        for k,v in kwargs.items():
            # if v is None: continue
            setattr(self, k, v)

    def __str__(self):
        _ret = "LoadOpts:\n"
        _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
        return _ret


def run_telescope_load(args):
    opts = LoadOpts(**vars(args))

    if opts.verbose:
        print(opts, file=sys.stderr)

    # Formats floats to requested precision
    if opts.float:
        _prec_str = '%%.%df' % opts.prec
    elif opts.exp:
        _prec_str = '%%.%de' % opts.prec
    else:
        _prec_str = '%%.%dg' % opts.prec

    with open(opts.checkpoint,'r') as fh:
        new_tm = TelescopeModel.load(fh)

    if opts.verbose:
        print("Checkpoint %s loaded successfully" % opts.checkpoint, file=sys.stderr)

    if opts.outparam is None:
        print("No output was requested", file=sys.stderr)
        return

    if opts.outparam == 'pi':
        if new_tm.pi is None:
            print("ERROR: pi is not set in this model", file=sys.stderr)
            return
        print('\t'.join(['genome','pi']), file=opts.outfile)
        for tup in zip(new_tm.txnames, new_tm.pi):
            print(tup[0] + '\t' + _prec_str % tup[1], file=opts.outfile)
    elif opts.outparam == 'pi_0':
        if new_tm.pi_0 is None:
            print("ERROR: pi_0 is not set in this model", file=sys.stderr)
            return
        print('\t'.join(['genome','pi_0']), file=opts.outfile)
        for tup in zip(new_tm.txnames, new_tm.pi_0):
            print(tup[0] + '\t' + _prec_str % tup[1], file=opts.outfile)

    elif opts.outparam == 'theta':
        if new_tm.theta is None:
            print("ERROR: theta is not set in this model", file=sys.stderr)
            return
        print('\t'.join(['genome','theta']), file=opts.outfile)
        for tup in zip(new_tm.txnames, new_tm.theta):
            print(tup[0] + '\t' + _prec_str % tup[1], file=opts.outfile)

    elif opts.outparam == 'x_hat':
        if new_tm.x_hat is None:
            print("ERROR: x_hat is not set in this model", file=sys.stderr)
            return
        print('\t'.join(['readname'] + new_tm.txnames), file=opts.outfile)
        for i,r in enumerate(new_tm.readnames):
            fmtvals = [_prec_str % v for v in new_tm.x_hat.getrow(i).toarray()[0,:]]
            print('\t'.join([r] + fmtvals), file=opts.outfile)
    elif opts.outparam == 'x_init':
        if new_tm.x_init is None:
            print("ERROR: x_init is not set in this model", file=sys.stderr)
            return
        print('\t'.join(['readname'] + new_tm.txnames), file=opts.outfile)
        for i,r in enumerate(new_tm.readnames):
            fmtvals = [_prec_str % v for v in new_tm.x_init.getrow(i).toarray()[0,:]]
            print('\t'.join([r] + fmtvals), file=opts.outfile)
    elif opts.outparam == 'Q':
        if new_tm.Q is None:
            print("ERROR: Q is not set in this model", file=sys.stderr)
            return
        print('\t'.join(['readname'] + new_tm.txnames), file=opts.outfile)
        for i,r in enumerate(new_tm.readnames):
            fmtvals = [_prec_str % v for v in new_tm.Q.getrow(i).toarray()[0,:]]
            print('\t'.join([r] + fmtvals), file=opts.outfile)
    else:
        print("Unknown output parameter: %s" % opts.outparam, file=sys.stderr)

    return
