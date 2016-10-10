# -*- coding: utf-8 -*-
""" Telescope load

"""

import sys

from utils.model import TelescopeModel

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


class LoadOpts:
    option_fields = ['verbose', 'checkpoint', 'outparam', 'prec', 'float', 'exp', 'outfile', ]

    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            # if v is None: continue
            setattr(self, k, v)

    def __str__(self):
        _ret = "LoadOpts:\n"
        _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
        return _ret


def run_telescope_load(args):
    opts = LoadOpts(**vars(args))

    if opts.verbose:
        print >>sys.stderr, opts

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
        print >>sys.stderr, "Checkpoint %s loaded successfully" % opts.checkpoint

    if opts.outparam is None:
        print >> sys.stderr, "No output was requested"
        return

    if opts.outparam == 'pi':
        if new_tm.pi is None:
            print >> sys.stderr, "ERROR: pi is not set in this model"
            return
        print >>opts.outfile, '\t'.join(['genome','pi'])
        for tup in zip(new_tm.txnames, new_tm.pi):
            print >>opts.outfile, tup[0] + '\t' + _prec_str % tup[1]
    elif opts.outparam == 'pi_0':
        if new_tm.pi_0 is None:
            print >> sys.stderr, "ERROR: pi_0 is not set in this model"
            return
        print >>opts.outfile, '\t'.join(['genome','pi_0'])
        for tup in zip(new_tm.txnames, new_tm.pi_0):
            print >>opts.outfile, tup[0] + '\t' + _prec_str % tup[1]

    elif opts.outparam == 'theta':
        if new_tm.theta is None:
            print >> sys.stderr, "ERROR: theta is not set in this model"
            return
        print >>opts.outfile, '\t'.join(['genome','theta'])
        for tup in zip(new_tm.txnames, new_tm.theta):
            print >>opts.outfile, tup[0] + '\t' + _prec_str % tup[1]

    elif opts.outparam == 'x_hat':
        if new_tm.x_hat is None:
            print >> sys.stderr, "ERROR: x_hat is not set in this model"
            return
        print >>opts.outfile, '\t'.join(['readname'] + new_tm.txnames)
        for i,r in enumerate(new_tm.readnames):
            fmtvals = [_prec_str % v for v in new_tm.x_hat.getrow(i).toarray()[0,:]]
            print >>opts.outfile, '\t'.join([r] + fmtvals)
    elif opts.outparam == 'x_init':
        if new_tm.x_init is None:
            print >> sys.stderr, "ERROR: x_init is not set in this model"
            return
        print >>opts.outfile, '\t'.join(['readname'] + new_tm.txnames)
        for i,r in enumerate(new_tm.readnames):
            fmtvals = [_prec_str % v for v in new_tm.x_init.getrow(i).toarray()[0,:]]
            print >>opts.outfile, '\t'.join([r] + fmtvals)
    elif opts.outparam == 'Q':
        if new_tm.Q is None:
            print >> sys.stderr, "ERROR: Q is not set in this model"
            return
        print >>opts.outfile, '\t'.join(['readname'] + new_tm.txnames)
        for i,r in enumerate(new_tm.readnames):
            fmtvals = [_prec_str % v for v in new_tm.Q.getrow(i).toarray()[0,:]]
            print >>opts.outfile, '\t'.join([r] + fmtvals)
    else:
        print >> sys.stderr, "Unknown output parameter: %s" % opts.outparam

    return
