__author__ = 'bendall'

import os, sys
from time import time

import pysam

import utils
from utils.alignment_parsers import TelescopeRead
from utils.annotation_parsers import AnnotationLookup
from utils.model import TelescopeModel

class IDOpts:
    option_fields = ['verbose','score_cutoff','out_matrix','no_updated_sam','emEpsilon','maxIter','piPrior','thetaPrior',
                     'exp_tag','outdir','ali_format','samfile','gtffile', ]

    def __init__(self, **kwargs):
        for k,v in kwargs.iteritems():
            if v is None: continue
            setattr(self, k, v)

    def generate_filename(self,suffix):
        basename = '%s-%s' % (self.exp_tag, suffix)
        return os.path.join(self.outdir,basename)

    def __str__(self):
        _ret = "IDOpts:\n"
        _ret += '\n'.join('  %s%s' % (f.ljust(30),getattr(self,f)) for f in self.option_fields)
        return _ret

def load_alignment(samfile, flookup, opts=None):
  _verbose = opts.verbose if opts is not None else True

  counts = {'unmapped':0, 'nofeat':0, 'mapped':0}
  mapped   = {}

  # Lookup reference name from reference ID
  refnames = dict(enumerate(samfile.references))

  for rname,segments in utils.iterread(samfile):
    r = TelescopeRead(rname,segments)
    if r.is_unmapped:
      counts['unmapped'] += 1
    else:
      r.assign_feats(refnames, flookup)
      if r.aligns_to_feat():
        r.assign_best()
        mapped[rname] = r
        counts['mapped'] += 1
      else:
        counts['nofeat'] += 1

    if _verbose and sum(counts.values()) % 100000 == 0:
      print >>sys.stderr, "...Processed %d fragments" % sum(counts.values())

  if _verbose:
    print >>sys.stderr, "Processed %d fragments" % sum(counts.values())
    print >>sys.stderr, "\t%d fragments were unmapped" % counts['unmapped']
    print >>sys.stderr, "\t%d fragments mapped to one or more positions on reference genome" % (counts['mapped'] + counts['nofeat'])
    print >>sys.stderr, "\t\t%d fragments mapped to reference but did not map to annotation" % counts['nofeat']
    print >>sys.stderr, "\t\t%d fragments have at least one alignment within annotation" % counts['mapped']

  return mapped

def run_telescope_id(args):
    opts = IDOpts(**vars(args))
    opts.verbose = True

    flookup = AnnotationLookup(opts.gtffile)
    samfile = pysam.AlignmentFile(opts.samfile)

    if opts.verbose:
        print >>sys.stderr, opts
        print >>sys.stderr, "Loading alignment file (%s)" % opts.samfile
        loadstart = time()

    mapped = load_alignment(samfile, flookup, opts)

    if opts.verbose:
        print >>sys.stderr, "Time to load alignment:".ljust(40) + "%d seconds" % (time() - loadstart)

    # Build the matrix
    ridx = {}
    gidx = {}
    d = []
    for rname,r in mapped.iteritems():
        i = ridx.setdefault(rname,len(ridx))
        for gname,aln in r.feat_aln_map.iteritems():
            j = gidx.setdefault(gname,len(gidx))
            d.append((i, j, aln.AS + aln.query_length))

    tm = TelescopeModel(d,ridx,gidx)

    if opts.verbose:
        print >>sys.stderr, "EM iteration..."
        print >>sys.stderr, "(Reads,Genomes)=%dx%d" % (tm.shape)
        print >>sys.stderr, "Delta Change:"
        emtime = time()

    tm.pi_0, tm.pi, tm.theta, tm.x_hat = utils.matrix_em(tm.Q, opts)

    if opts.verbose:
        print >>sys.stderr, "Time for EM iteration:".ljust(40) +  "%d seconds" % (time() - emtime)

    """ Output results """

    # Results report
    report = tm.make_report()
    with open(opts.generate_filename('telescope_report.tsv'),'w') as outh:
      for row in report:
        print >>outh, '\t'.join(str(f) for f in row)

    # Output matrix
    if opts.out_matrix:
        with open(opts.generate_filename('xmat_initial.txt'),'w') as outh:
          print >>outh, tm.x_init.pretty_tsv(tm.rownames, tm.colnames)
        with open(opts.generate_filename('xmat_final.txt'),'w') as outh:
          print >>outh, tm.x_hat.pretty_tsv(tm.rownames, tm.colnames)

    """
    #tm.x_hat.max_row_matrix()

    t_0 = utils.telescope_best_hit(tm.x_init)
    t_F = utils.telescope_best_hit(tm.x_hat)

    best_assignments = utils.get_best_assignment(tm.x_hat)

    R,G = tm.shape
    report_data = {'genome': [k for k,v in sorted(gidx.iteritems(),key=lambda x:x[1])],
                   'init_pi':tm.pi_0,    'final_pi':tm.pi,
                   'init_best':t_0[0],  'init_l1': t_0[1],  'init_l2':t_0[2],  'init_best_hit':[float(v)/R for v in t_0[0]] ,
                   'final_best':t_F[0], 'final_l1': t_F[1], 'final_l2':t_F[2], 'final_best_hit':[float(v)/R for v in t_0[0]] ,
                   'weighted_counts': tm.calculate_weighted_counts(),
                   # 'fractional_counts': tm.calculate_fractional_counts(),
                   'unique_counts': tm.calculate_unique_counts(),
    }
    with open(opts.generate_filename('rs_report.tsv'),'w') as outh:
      for row in utils.telescope_report(report_data,R,G):
        print >>outh, '\t'.join(str(f) for f in row)

    if opts.out_matrix:
        with open(opts.generate_filename('xmat_initial.txt'),'w') as outh:
          print >>outh, utils.pretty_xmat(tm, tm.x_init)
        with open(opts.generate_filename('xmat_final.txt'),'w') as outh:
          print >>outh, utils.pretty_xmat(tm, tm.x_hat)
    """

    if not opts.no_updated_sam:
        updated_samfile = pysam.AlignmentFile(opts.generate_filename('updated.sam'), 'wh', header=samfile.header)
        mapq = tm.x_hat.apply_func(utils.phred)
        with open(opts.generate_filename('mapq.txt'),'w') as outh:
            print >>outh, mapq.pretty_tsv(tm.rownames, tm.colnames)
        for rname,r in mapped.iteritems():
            nzidx = mapq[tm.ridx[rname],].nonzero()[1]
            nzval = [mapq[tm.ridx[rname],v] for v in nzidx]
            nonzero = sorted(zip(nzidx,nzval),key=lambda x:x[1],reverse=True)
            for gidx,mq in nonzero:
                primary,alternates = r.aligned_to_genome(tm.colnames[gidx])
                primary.set_mapq(mq)
                primary.set_tag('ZF',tm.colnames[gidx])
                primary.write_samfile(updated_samfile)
                for altaln in alternates:
                    altaln.set_mapq(0)
                    altaln.set_tag('ZF',tm.colnames[gidx])
                    altaln.write_samfile(updated_samfile)

        updated_samfile.close()

    samfile.close()
    return
