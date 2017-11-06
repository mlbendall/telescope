# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import absolute_import

import sys
import os
import logging as lg
from collections import OrderedDict, defaultdict, Counter
import gc

import numpy as np
import scipy
import pysam


from .annotation import Annotation
from .sparse_plus import csr_matrix_plus as csr_matrix
from .colors import c2str, D2PAL, GPAL


__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


def process_overlap_frag(pairs, overlap_feats):
    ''' Find the best alignment for each locus '''
    assert all(pairs[0].query_id == p.query_id for p in pairs)
    ''' Organize by feature'''
    byfeature = defaultdict(list)
    for pair, feat in zip(pairs, overlap_feats):
        byfeature[feat].append(pair)

    _maps = []
    for feat, falns in byfeature.items():
        falns.sort(key=lambda x: x.alnscore + x.alnlen,
                  reverse=True)
        _topaln = falns[0]
        # Add best to mappings
        _maps.append(
            (_topaln.query_id, feat, _topaln.alnscore, _topaln.alnlen)
        )
        # Set tags
        falns[0].set_tag('ZF', feat)
        for aln in falns[1:]:
            aln.set_tag('ZF', feat)
            aln.set_tag('ZT', 'SEC')

    _maps.sort(key=lambda x: x[2], reverse=True)
    _topfeat = ','.join(t[1] for t in _maps if t[2] == _maps[0][2])
    for p in pairs:
        p.set_tag('ZB', _topfeat)

    return _maps

from .alignment import fetch_fragments

class Telescope(object):
    """

    """
    def __init__(self, opts):

        self.opts = opts               # Command line options
        self.run_info = OrderedDict()  # Information about the run
        self.annotation = None         # Anntation object
        self.feature_length = None    # Lengths of features
        self.read_index = {}           # {"fragment name": row_index}
        self.feat_index = {}           # {"feature_name": column_index}
        self.shape = None              # Fragments x Features
        self.raw_scores = None         # Initial alignment scores

        # BAM with non overlapping fragments (or unmapped)
        self.other_bam = opts.outfile_path('other.bam')
        # BAM with overlapping fragments
        self.tmp_bam = opts.outfile_path('tmp_tele.bam')

        # Set the version
        self.run_info['version'] = self.opts.version

    def cleanup(self):
        if os.path.exists(self.tmp_bam):
            os.unlink(self.tmp_bam)

    def load_annotation(self):
        self.annotation = Annotation(self.opts.gtffile, self.opts.attribute)
        self.run_info['annotated_features'] = len(self.annotation.loci)
        self.feature_length = self.annotation.feature_length().copy()

    def load_alignment(self):
        alninfo = Counter()

        _nfkey = self.opts.no_feature_key
        _update_sam = self.opts.updated_sam

        _mappings = []
        assign = Assigner(self.annotation, self.opts).get_assign_function()

        """ Load unsorted reads """
        with pysam.AlignmentFile(self.opts.samfile) as sf:
            # Create output temporary files
            if _update_sam:
                bam_u = pysam.AlignmentFile(self.other_bam, 'w', template=sf)
                bam_t = pysam.AlignmentFile(self.tmp_bam, 'w', template=sf)

            for pairs in fetch_fragments(sf, until_eof=True):
                alninfo['fragments'] += 1
                if alninfo['fragments'] % 500000 == 0:
                    lg.info('...processed {:.1f}M fragments'.format(alninfo['fragments']/1e6))

                ''' Check whether fragment is mapped '''
                if pairs[0].is_unmapped:
                    alninfo['unmap_{}'.format(pairs[0].numreads)] += 1
                    if _update_sam: pairs[0].write(bam_u)
                    continue

                ''' Fragment is mapped '''
                alninfo['map_{}'.format(pairs[0].numreads)] += 1

                ''' Fragment is ambiguous if multiple mappings'''
                _ambig = len(pairs) > 1

                ''' Check whether fragment overlaps annotation '''
                overlap_feats = list(map(assign, pairs))
                has_overlap = any(f != _nfkey for f in overlap_feats)

                ''' Fragment has no overlap '''
                if not has_overlap:
                    alninfo['nofeat_{}'.format('A' if _ambig else 'U')] += 1
                    if _update_sam:
                        [p.write(bam_u) for p in pairs]
                    continue

                ''' Fragment overlaps with annotation '''
                alninfo['feat_{}'.format('A' if _ambig else 'U')] += 1

                ''' Find the best alignment for each locus '''
                _mappings += process_overlap_frag(pairs, overlap_feats)

                if _update_sam:
                    [p.write(bam_t) for p in pairs]

        ''' Loading complete '''
        self.run_info['total_fragments'] = alninfo['fragments']
        self.run_info['mapped_pairs'] = alninfo['map_2']
        self.run_info['mapped_single'] = alninfo['map_1']
        self.run_info['unmapped'] = alninfo['unmap_2'] + alninfo['unmap_1']
        self.run_info['unique'] = alninfo['nofeat_U'] + alninfo['feat_U']
        self.run_info['ambig'] = alninfo['nofeat_A'] + alninfo['feat_A']
        self.run_info['overlap_unique'] = alninfo['feat_U']
        self.run_info['overlap_ambig'] = alninfo['feat_A']

        if _update_sam:
            bam_u.close()
            bam_t.close()

        self._mapping_to_matrix(_mappings)

    def load_mappings(self, samfile_path):
        _mappings = []
        with pysam.AlignmentFile(samfile_path) as sf:
            for pairs in fetch_fragments(sf, until_eof=True):
                for pair in pairs:
                    if pair.r1.has_tag('ZT'):
                        continue
                    _mappings.append((
                        pair.query_id,
                        pair.r1.get_tag('ZF'),
                        pair.alnscore,
                        pair.alnlen
                    ))
        return _mappings

    # @profile
    def _mapping_to_matrix(self, mappings):
        ''' '''
        _maxAS = max(t[2] for t in mappings)
        _minAS = min(t[2] for t in mappings)

        # Rescale integer alignment score to be greater than zero
        rescale = {s: (s - _minAS + 1) for s in range(_minAS, _maxAS + 1)}

        # Construct dok matrix with mappings
        if 'annotated_features' in self.run_info:
            ncol = self.run_info['annotated_features']
        else:
            ncol = len(set(t[1] for t in mappings))
        dim = (len(mappings), ncol)
        _m1 = scipy.sparse.dok_matrix(dim, dtype=np.uint16)
        _ridx = self.read_index
        _fidx = self.feat_index
        for rid, fid, ascr, alen in mappings:
            i = _ridx.setdefault(rid, len(_ridx))
            j = _fidx.setdefault(fid, len(_fidx))
            _m1[i, j] = max(_m1[i, j], (rescale[ascr] + alen))

        # Trim matrix to size
        _m1 = _m1[:len(_ridx), :len(_fidx)]

        # Convert dok matrix to csr
        self.raw_scores = csr_matrix(_m1)
        self.shape = (len(_ridx), len(_fidx))

    def output_report(self, tl, filename):
        _rmethod, _rprob = self.opts.reassign_mode, self.opts.conf_prob
        _fnames = sorted(self.feat_index, key=self.feat_index.get)
        _flens = self.feature_length #self.annotation.feature_length()
        _final_type = '{:.2f}' if _rmethod in ['average', 'conf'] else '{:d}'
        _dtype = [
            ('transcript', '{:s}'),
            ('transcript_length', '{:d}'),
            ('final_count', _final_type),
            ('final_conf', '{:.2f}'),
            ('final_prop', '{:.3g}'),
            ('init_aligned', '{:d}'),
            ('unique_count', '{:d}'),
            ('init_best', '{:d}'),
            ('init_best_random', '{:d}'),
            ('init_best_avg', '{:.2f}'),
            ('init_prop', '{:.3g}'),
        ]

        _report0 = [
            _fnames,                                       # transcript
            [_flens[f] for f in _fnames],                  # tx_len
            tl.reassign(_rmethod, _rprob).sum(0).A1,       # final_count
            tl.reassign('conf', _rprob).sum(0).A1,         # final_conf
            tl.pi[-1],                                     # final_prop
            tl.reassign('all', iteration=0).sum(0).A1,     # init_aligned
            tl.reassign('unique').sum(0).A1,               # unique_count
            tl.reassign('exclude', iteration=0).sum(0).A1, # init_best
            tl.reassign('choose', iteration=0).sum(0).A1,  # init_best_random
            tl.reassign('average', iteration=0).sum(0).A1, # init_best_avg
            tl.pi[1]                                       # init_prop
        ]

        # Rotate the report
        _report = [[r0[i] for r0 in _report0] for i in range(len(_fnames))]

        # Sort the report
        _report.sort(key=lambda x: x[4], reverse=True)
        _report.sort(key=lambda x: x[2], reverse=True)

        _fmtstr = '\t'.join(t[1] for t in _dtype)

        # Run info line
        _comment = ["## RunInfo", ]
        _comment += ['{}:{}'.format(*tup) for tup in self.run_info.items()]

        with open(filename, 'w') as outh:
            print('\t'.join(_comment), file=outh)
            print('\t'.join(t[0] for t in _dtype), file=outh)
            for row in _report:
                print(_fmtstr.format(*row), file=outh)
        return

    def update_sam(self, tl, filename):
        _rmethod, _rprob = self.opts.reassign_mode, self.opts.conf_prob
        _fnames = sorted(self.feat_index, key=self.feat_index.get)

        mat = csr_matrix(tl.reassign(_rmethod, _rprob))
        # best_feats = {i: _fnames for i, j in zip(*mat.nonzero())}

        with pysam.AlignmentFile(self.tmp_bam) as sf:
            header = sf.header
            header['PG'].append({
                'PN': 'telescope', 'ID': 'telescope',
                'VN': self.run_info['version'],
                'CL': ' '.join(sys.argv),
            })
            outsam = pysam.AlignmentFile(filename, 'wb', header=header)
            for pairs in fetch_fragments(sf, until_eof=True):
                if len(pairs) == 0: continue
                ridx = self.read_index[pairs[0].query_id]
                for aln in pairs:
                    if aln.r1.has_tag('ZT'):
                        aln.set_tag('YC', c2str((248, 248, 248)))
                        aln.set_mapq(0)
                    else:
                        fidx = self.feat_index[aln.r1.get_tag('ZF')]
                        prob = tl.z[-1][ridx, fidx]
                        aln.set_tag('XP', int(round(prob*100)))
                        if mat[ridx, fidx] > 0:
                            aln.unset_flag(pysam.FSECONDARY)
                            aln.set_tag('YC',c2str(D2PAL['vermilion']))
                        else:
                            aln.set_flag(pysam.FSECONDARY)
                            if prob >= 0.2:
                                aln.set_tag('YC', c2str(D2PAL['yellow']))
                            else:
                                aln.set_tag('YC', c2str(GPAL[2]))
                    aln.write(outsam)
            outsam.close()

    def __str__(self):
        _ret = '<Telescope samfile=%s, gtffile=%s>' % (self.opts.samfile, self.opts.gtffile)
        return _ret


class TelescopeLikelihood(object):
    """

    """
    def __init__(self, score_matrix, opts):
        """
        """
        # Raw scores
        self.raw_scores = score_matrix
        self.max_score = self.raw_scores.max()

        # N fragments x K transcripts
        self.N, self.K = self.raw_scores.shape

        # Q[i,] is the set of mapping qualities for fragment i, where Q[i,j]
        # represents the evidence for fragment i being generated by fragment j.
        # In this case the evidence is represented by an alignment score, which
        # is greater when there are more matches and is penalized for
        # mismatches
        # Scale the raw alignment score by the maximum alignment score
        # and multiply by a scale factor.
        self.scale_factor = 100.
        self.Q = self.raw_scores.scale().multiply(self.scale_factor).expm1()
        # self.Q = self.raw_scores.multiply(self.scale_factor).expm1()

        # z[i,] is the partial assignment weights for fragment i, where z[i,j]
        # is the expected value for fragment i originating from transcript j. The
        # initial estimate is the normalized mapping qualities:
        # z_init[i,] = Q[i,] / sum(Q[i,])
        self.z = [ self.Q.norm(1) ]

        self.epsilon = opts.em_epsilon
        self.max_iter = opts.max_iter

        # pi[j] is the proportion of fragments that originate from
        # transcript j. Initial value assumes that all transcripts contribute
        # equal proportions of fragments
        self.pi = [ np.repeat(1./self.K, self.K) ]

        # theta[j] is the proportion of non-unique fragments that need to be
        # reassigned to transcript j. Initial value assumes that all transcripts
        # are reassigned an equal proportion of fragments
        self.theta = [ np.repeat(1./self.K, self.K) ]

        # Y[i] is the ambiguity indicator for fragment i, where Y[i]=1 if
        # fragment i is aligned to multiple transcripts and Y[i]=0 otherwise.
        # Store as N x 1 matrix
        self.Y = (self.Q.count(1) > 1).astype(np.int)

        # Log-likelihood score
        self.lnl = []

        # Prior values
        self.pi_prior = opts.pi_prior
        self.theta_prior = opts.theta_prior

        # Precalculated values
        self._weights = self.Q.max(1)             # Weight assigned to each fragment
        self._total_wt = self._weights.sum()      # Total weight
        self._ambig_wt = self._weights.multiply(self.Y).sum() # Weight of ambig frags
        self._unique_wt = self._weights.multiply(1-self.Y).sum()

        # Weighted prior values
        self._pi_prior_wt = self.pi_prior * self._weights.max()
        self._theta_prior_wt = self.theta_prior * self._weights.max()
        #
        self._pisum0 = self.Q.multiply(1-self.Y).sum(0)

    def estep(self):
        """ Calculate the expected values of z
                E(z[i,j]) = ( pi[j] * theta[j]**Y[i] * Q[i,j] ) /
        """
        # assert len(self.z) == len(self.pi) == len(self.theta)
        _pi = self.pi[-1]
        _theta = self.theta[-1]
        _numerator = self.Q.multiply(csr_matrix(_pi * (_theta ** self.Y)))
        self.z.append(_numerator.norm(1))

    def mstep(self):
        """ Calculate the maximum a posteriori (MAP) estimates for pi and theta

        """
        # assert (len(self.z)-1) == len(self.pi) == len(self.theta)
        # The expected values of z weighted by mapping score
        _weighted = self.z[-1].multiply(self._weights)

        # Estimate theta_hat
        _thetasum = _weighted.multiply(self.Y).sum(0)
        _theta_denom = self._ambig_wt + self._theta_prior_wt * self.K
        _theta_hat = (_thetasum + self._theta_prior_wt) / _theta_denom

        # Estimate pi_hat
        _pisum = self._pisum0 + _thetasum
        _pi_denom = self._ambig_wt + self._unique_wt + self._pi_prior_wt * self.K
        _pi_hat = (_pisum + self._pi_prior_wt) / _pi_denom

        self.theta.append(_theta_hat.A1)
        self.pi.append(_pi_hat.A1)

    def calculate_lnl(self):
        _z, _p, _t = self.z[-1], self.pi[-1], self.theta[-1]
        cur = _z.multiply(self.Q.multiply(_p * _t**self.Y).log1p()).sum()
        self.lnl.append(cur)

    # @profile
    def em(self, use_likelihood=False, loglev=lg.WARNING):
        msg = 'Iteration %d, lnl=%g, diff=%g'
        converged = False
        while not converged:
            self.estep()
            self.mstep()
            self.calculate_lnl()

            iidx = len(self.lnl)
            if iidx <= 1:
                diff_lnl = float('inf')
            else:
                diff_lnl = abs(self.lnl[-1] - self.lnl[-2])
            diff_est = abs(self.pi[-1] - self.pi[-2]).sum()

            lg.log(loglev, msg % (iidx, self.lnl[-1], diff_est))
            if use_likelihood:
                converged = diff_lnl < self.epsilon or iidx >= self.max_iter
            else:
                converged = diff_est < self.epsilon or iidx >= self.max_iter

    def reassign(self, method, thresh=0.9, iteration=-1):
        """ Reassign fragments to expected transcripts

        Running EM finds the expected fragment assignment weights at the MAP
        estimates of pi and theta. This function reassigns all fragments based
        on these assignment weights. A simple heuristic is to assign each
        fragment to the transcript with the highest assignment weight.

        In practice, not all fragments have exactly one best hit. The "method"
        argument defines how we deal with fragments that are not fully resolved
        after EM:
                exclude - reads with > 1 best hits are excluded
                choose  - one of the best hits is randomly chosen
                average - read is evenly divided among best hits
                conf    - only confident reads are reassigned
                unique  - only uniquely aligned reads
        Args:
            method:
            thresh:
            iteration:

        Returns:
            matrix where m[i,j] == 1 iff read i is reassigned to transcript j

        """
        _z = self.z[iteration]
        if method == 'exclude':
            # Identify best hit(s), then exclude rows with >1 best hits
            v = _z.binmax(1)
            return v.multiply(v.sum(1) == 1)
        elif method == 'choose':
            # Identify best hit(s), then randomly choose reassignment
            v = _z.binmax(1)
            return v.choose_random(1)
        elif method == 'average':
            # Identify best hit(s), then divide by row sum
            v = _z.binmax(1)
            return v.norm(1)
        elif method == 'conf':
            # Zero out all values less than threshold
            # If thresh > 0.5 then at most
            v = _z.apply_func(lambda x: x if x >= thresh else 0)
            # Average each row so each sums to 1.
            return v.norm(1)
        elif method == 'unique':
            # Zero all rows that are ambiguous
            return _z.multiply(1 - self.Y).ceil().astype(np.uint8)
        elif method == 'all':
            # Return all nonzero elements
            return _z.apply_func(lambda x: 1 if x > 0 else 0).astype(np.uint8)


class Assigner:
    def __init__(self, annotation, opts):
        self.annotation = annotation
        self.opts = opts

    def get_assign_function(self):
        def _assign_pair_threshold(pair):
            blocks = pair.refblocks
            f = self.annotation.intersect_blocks(pair.ref_name, blocks)
            if not f:
                return self.opts.no_feature_key
            # Calculate the percentage of fragment mapped
            fname, overlap = f.most_common()[0]
            if overlap > pair.alnlen * self.opts.overlap_threshold:
                return fname
            else:
                return self.opts.no_feature_key

        def _assign_pair_intersection_strict(pair):
            pass

        def _assign_pair_union(pair):
            pass

        ''' Return function depending on overlap mode '''
        if self.opts.overlap_mode == 'threshold':
            return _assign_pair_threshold
        elif self.opts.overlap_mode == 'intersection-strict':
            return _assign_pair_intersection_strict
        elif self.opts.overlap_mode == 'union':
            return _assign_pair_union
        else:
            assert False