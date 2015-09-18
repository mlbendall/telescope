__author__ = 'bendall'

import os, sys
from time import time

import pysam

import utils
from utils.alignment_parsers import TelescopeRead
from utils.annotation_parsers import AnnotationLookup
from utils.model import TelescopeModel
from utils.colors import c2str, DARK2_PALETTE, GREENS
from utils import format_minutes

class IDOpts:
    option_fields = ['ali_format','samfile','gtffile',
                     'verbose', 'outdir', 'exp_tag', 'out_matrix', 'updated_sam', 'checkpoint',
                     'min_prob', 'conf_prob',
                     'piPrior', 'thetaPrior',
                     'emEpsilon','maxIter',
                     ]

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

        if _verbose and sum(counts.values()) % 500000 == 0:
            print >>sys.stderr, "...Processed %d fragments" % sum(counts.values())

    if _verbose:
        print >>sys.stderr, "Processed %d fragments" % sum(counts.values())
        print >>sys.stderr, "\t%d fragments were unmapped" % counts['unmapped']
        print >>sys.stderr, "\t%d fragments mapped to one or more positions on reference genome" % (counts['mapped'] + counts['nofeat'])
        print >>sys.stderr, "\t\t%d fragments mapped to reference but did not map to annotation" % counts['nofeat']
        print >>sys.stderr, "\t\t%d fragments have at least one alignment within annotation" % counts['mapped']

    return mapped

def update_alignment(tm, mapped, newsam, min_prob=0.1, conf_prob=0.9):
    # Read x Genome matrix = 1 if output alignment 0 otherwise
    output_mat = tm.x_hat.apply_func(lambda x: 1 if x >= min_prob else 0)
    for rownum in xrange(output_mat.shape[0]):
        _rname = tm.rownames[rownum]
        _read  = mapped[_rname]
        try:
            # Sorted list of (genome_index, prob)
            gidx_probs = sorted(((_, tm.x_hat[rownum,_]) for _ in output_mat[rownum,].nonzero()[1]), key=lambda x:x[1],reverse=True)
            best_genome = tm.colnames[gidx_probs[0][0]]
            for colnum,prob in gidx_probs:
                genome_name = tm.colnames[colnum]
                primary, alternates = _read.aligned_to_genome(genome_name)

                # Print primary alignment
                primary.set_mapq(utils.phred(prob))
                primary.set_tag('XP',int(round(prob*100)))
                primary.set_tag('XF', genome_name)
                primary.set_tag('ZF', best_genome)

                if genome_name == best_genome:            # best genome
                    primary.set_secondary(False)
                    if len(gidx_probs)==1:                    # only one genome has prob > min_prob
                        if prob >= conf_prob:                     # high confidence
                            primary.set_tag('YC', c2str(DARK2_PALETTE['vermilion']))
                        else:                                     # low confidence
                            primary.set_tag('YC', c2str(DARK2_PALETTE['yellow']))
                    else:                                     # multiple genomes have prob > min_prob
                        primary.set_tag('YC', c2str(DARK2_PALETTE['teal']))
                        assert prob < .9, "If there are multiple nonzero genomes, qual must be < .9"
                else:
                    primary.set_tag('YC', c2str(GREENS[2]))    # Not best genome
                    primary.set_secondary(True)

                primary.write_samfile(newsam)

                # Print alternate alignments
                for altaln in alternates:
                    altaln.set_mapq(0)
                    altaln.set_tag('XP',0)
                    altaln.set_tag('XF', genome_name)
                    altaln.set_tag('ZF', best_genome)
                    altaln.set_tag('YC', c2str((248,248,248)))
                    altaln.set_secondary(True)
                    altaln.write_samfile(newsam)
        except IndexError:
            print >>sys.stderr, "Unable to write %s" % _rname

def run_telescope_id(args):
    opts = IDOpts(**vars(args))
    if opts.verbose:
        print >>sys.stderr, opts

    """ Load alignment """
    if opts.verbose:
        print >>sys.stderr, "Loading alignment file (%s):" % opts.samfile
        substart = time()

    flookup = AnnotationLookup(opts.gtffile)
    samfile = pysam.AlignmentFile(opts.samfile)
    mapped = load_alignment(samfile, flookup, opts)

    if opts.verbose:
        print >>sys.stderr, "Time to load alignment:".ljust(40) + format_minutes(time() - substart)

    """ Create data structure """
    if opts.verbose:
        print >>sys.stderr, "Creating data structure... " ,
        substart = time()

    ridx = {}
    gidx = {}
    d = []
    for rname,r in mapped.iteritems():
        i = ridx.setdefault(rname,len(ridx))
        for gname,aln in r.feat_aln_map.iteritems():
            j = gidx.setdefault(gname,len(gidx))
            d.append((i, j, aln.AS + aln.query_length))

    tm = TelescopeModel(ridx, gidx, data=d)

    if opts.verbose:
        print >>sys.stderr, "done."
        print >>sys.stderr, "Time to create data structure:".ljust(40) + format_minutes(time() - substart)

    # Save some memory if you are not creating an updated SAM:
    if opts.no_updated_sam:
        if opts.verbose:
            print >>sys.stderr, "Clearing reads from memory."
        mapped = None
        samfile.close()

    """ Checkpoint 1 """
    if opts.checkpoint:
        if opts.verbose:
            print >>sys.stderr, "Checkpointing... " ,
            substart = time()

        with open(opts.generate_filename('checkpoint.pickle'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write checkpoint:".ljust(40) + format_minutes(time() - substart)

    """ Initial output matrix """
    if opts.out_matrix:
        if opts.verbose:
            print >>sys.stderr, "Writing initial model matrix...",
            substart = time()

        with open(opts.generate_filename('model-initial.pickle'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write initial matrix:".ljust(40) + format_minutes(time() - substart)

    """ Reassignment """
    if opts.verbose:
        print >>sys.stderr, "Reassiging reads:"
        print >>sys.stderr, "(Reads,Genomes)=%dx%d" % (tm.shape)
        print >>sys.stderr, "Delta Change:"
        substart = time()

    tm.pi_0, tm.pi, tm.theta, tm.x_hat = utils.matrix_em(tm.Q, opts)

    if opts.verbose:
        print >>sys.stderr, "Time for EM iteration:".ljust(40) + format_minutes(time() - substart)

    """ Checkpoint 2 """
    if opts.checkpoint:
        if opts.verbose:
            print >>sys.stderr, "Checkpointing... " ,
            substart = time()

        with open(opts.generate_filename('checkpoint.pickle'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write checkpoint:".ljust(40) + format_minutes(time() - substart)

    """ Final output matrix """
    if opts.out_matrix:
        if opts.verbose:
            print >>sys.stderr, "Writing final model matrix...",
            substart = time()

        with open(opts.generate_filename('model-final.pickle'),'w') as outh:
            tm.dump(outh)

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to write final matrix:".ljust(40) + format_minutes(time() - substart)

    """ Generate report """
    if opts.verbose:
        print >>sys.stderr, "Generating report... " ,
        substart = time()

    report = tm.make_report(opts.conf_prob)
    with open(opts.generate_filename('telescope_report.tsv'),'w') as outh:
        for row in report:
            print >>outh, '\t'.join(str(f) for f in row)

    if opts.verbose:
        print >>sys.stderr, "done."
        print >>sys.stderr, "Time to generate report:".ljust(40) + format_minutes(time() - substart)

    """ Update alignment """
    if opts.updated_sam:
        if opts.verbose:
            print >>sys.stderr, "Updating alignment...",
            substart = time()

        updated_samfile = pysam.AlignmentFile(opts.generate_filename('updated.sam'), 'wh', header=samfile.header)
        update_alignment(tm, mapped, updated_samfile, min_prob=opts.min_prob, conf_prob=opts.conf_prob)
        updated_samfile.close()

        if opts.verbose:
            print >>sys.stderr, "done."
            print >>sys.stderr, "Time to update alignment:".ljust(40) + format_minutes(time() - substart)

    samfile.close()
    return
