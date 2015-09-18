__author__ = 'bendall'

import sys
from utils.model import TelescopeModel


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

    print >>sys.stderr, 'format string is "' + _prec_str + '"'
    with open(opts.checkpoint,'r') as fh:
        new_tm = TelescopeModel.load(fh)

    if opts.outparam is None:
        print >> sys.stderr, "No output was requested"
        return

    if opts.outparam == 'pi' or opts.outparam == 'theta':
        lol = [ new_tm.colnames ]
        colheader = ['genome']
        if new_tm.theta is not None:
            lol.append(new_tm.theta)
            colheader.append('theta')
        if new_tm.pi_0 is not None:
            lol.append(new_tm.pi_0)
            colheader.append('pi_0')
        if new_tm.pi is not None:
            lol.append(new_tm.pi)
            colheader.append('pi')

        print >>opts.outfile, '\t'.join(colheader)
        # Iterate over rows
        for tup in zip(*lol):
            l = [tup[0]] + [_prec_str % v for v in tup[1:]]
            print >>opts.outfile, '\t'.join(l)

    elif opts.outparam == 'x_init':
        print >>opts.outfile, '\t'.join(['readname'] + new_tm.colnames)
        for i,r in enumerate(new_tm.rownames):
            fmtvals = [_prec_str % v for v in new_tm.x_init.getrow(i).toarray()[0,:]]
            print >>opts.outfile, '\t'.join([r] + fmtvals)

    elif opts.outparam == 'x_hat':
        print >>opts.outfile, '\t'.join(['readname'] + new_tm.colnames)
        for i,r in enumerate(new_tm.rownames):
            fmtvals = [_prec_str % v for v in new_tm.x_hat.getrow(i).toarray()[0,:]]
            print >>opts.outfile, '\t'.join([r] + fmtvals)

    elif opts.outparam == 'Q':
        print >>opts.outfile, '\t'.join(['readname'] + new_tm.colnames)
        for i,r in enumerate(new_tm.rownames):
            fmtvals = [_prec_str % v for v in new_tm.Q.getrow(i).toarray()[0,:]]
            print >>opts.outfile, '\t'.join([r] + fmtvals)

    else:
        print >> sys.stderr, "Unknown output parameter: %s" % opts.outparam

    return
