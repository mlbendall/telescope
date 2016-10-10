# -*- coding: utf-8 -*-
""" Telescope model
"""

import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle

import numpy as np
import scipy.sparse

from sparse_matrix import csr_matrix_plus as csr_matrix

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"


class TelescopeModel:
    '''

    '''
    def __init__(self, read_index, tx_index, data=None, qmat=None):
        ''' Initialize TelescopeModel
        :param read_index: Dictionary mapping read names to row index
        :param tx_index: Dictionary mapping transcript names to column index
        :param data: List of tuples with raw mapping scores
        :param qmat: Sparse matrix with scaled mapping scores
        :return:
        '''
        # read_index is a dictionary mapping read name to row index
        # readnames is a sorted list of read names
        self.read_index = read_index
        self.readnames = [k for k,v in sorted(self.read_index.iteritems(), key=lambda x:x[1])]

        # tx_index is a dictionary mapping transcript name to column index
        # txnames is a sorted list of transcript names
        self.tx_index = tx_index
        self.txnames = [k for k,v in sorted(self.tx_index.iteritems(), key=lambda x:x[1])]

        # shape is the number of reads X number of transcripts
        self.shape = (len(self.read_index), len(self.tx_index))

        # Q[i,] is the scaled mapping scores for read i, where Q[i,j] is the
        # mapping score of read i aligned to transcript j.
        if data is not None:
            # Data provided as a list of tuples:
            # (read_index, transcript_index, alignment_score)
            i,j,d = zip(*data)
            _coo = scipy.sparse.coo_matrix((d,(i,j)), shape=self.shape)
            _raw_scores = csr_matrix(_coo)
            self.Q = _raw_scores.multiply(100.0 / _raw_scores.max()).exp()
        else:
            # Data provided as matrix (loaded from checkpoint)
            assert qmat is not None, "qmat must be provided if data is not"
            self.Q = qmat

        # x[i,] is the transcript indicator for read i, where x[i,j] is the
        # expected value for read i originating from transcript j. The initial
        # estimate of x[i,] (x_init) is the normalized mapping scores:
        # x_init[i,] = Q[i,] / sum(Q[i,])
        self.x_init = self.Q.normr()
        self.x_hat = None

        # Y[i] is the uniqueness indicator for read i, where Y[i]=1 if read i
        # maps uniquely (to only one transcript) and Y[i]=0 otherwise
        self.Y = np.where(self.Q.countr()==1, 1, 0)

        # pi[j] is the proportion of reads that originated from transcript j
        self.pi_0  = None
        self.pi    = None

        # theta[j] is the reassignment parameter representing the proportion
        # of non-unique reads that need to be reassigned to transcript j
        self.theta = None

        # Finally, seed the random number generator for consistent results
        # Function of the first 10 read names
        seed = sum([ord(c) for c in ''.join(self.readnames[:10])][::3])
        np.random.seed(seed)

    def reassign_to_best(self, method, thresh=None, initial=False):
        """ Returns matrix where m[i,j] == 1 iff read i is reassigned to transcript j

            "method" defines how we deal with reads that have more than one
            best transcript:
                exclude - reads with > 1 best hits are excluded
                choose  - one of the best hits is randomly chosen
                average - read is evenly divided among best hits
                conf    - only confident reads are reassigned
                unique  - only uniquely aligned reads
        """
        # Operate on initial or final matrix
        mat = self.x_init if initial else self.x_hat

        # Reads that are ambiguous are not assigned (set to zero)
        if method == 'exclude':
            # The maximum value in each row is set to 1. If the max value
            # appears more than once, all elements equal to the max value are
            # set to 1. Rows with sum greater than 1 have all elements set to 0
            v = mat.maxidxr(choose=False)
            return v.multiply(np.where(v.sumr()>1, 0, 1))

        # Reads that are ambiguous are randomly assigned
        if method == 'choose':
            # The maximum value in each row is set to 1. If the max value
            # appears more than once, one of these elements are randomly chosen
            # and set to 1.
            return mat.maxidxr(choose=True)

        # Reads that are ambiguous are evenly divided among best transcripts
        if method == 'average':
            # The maximum value in each row is set to 1, then divided by the
            # row sum. So, if the max value appears more than once, each
            # element will be set to a fractional value.
            v = mat.maxidxr(choose=False)
            return v.normr()

        # Reads are reassigned if they meet or exceed threshold
        if method == 'conf':
            # Elements in each row that are greater than thresh are set to 1.
            # If thresh > 0.5 then at most 1 element will be equal to 1.
            assert thresh is not None and thresh > 0.5, "Invalid value for thresh: %s" % conf
            f = lambda x: 1 if x >= thresh else 0
            return mat.apply_func(f)

        # Reads that are initially uniquely mapped
        if method == 'unique':
            # All non-zero values are set to 1, then multiplied by uniqueness
            # indicator
            return mat.ceil().multiply(csr_matrix(self.Y[:,None]))


        assert False, "Method is invalid: %s" % method

    def dump(self,fh):
        # Python objects
        pickle.dump([self.read_index, self.tx_index], fh)

        # csr_matrix
        self.Q.dump(fh)

        # Numpy arrays
        if self.pi_0 is None:
            pickle.dump(None, fh)
        else:
            self.pi_0.dump(fh)

        if self.pi is None:
            pickle.dump(None, fh)
        else:
            self.pi.dump(fh)

        if self.theta is None:
            pickle.dump(None, fh)
        else:
            self.theta.dump(fh)

       # csr_matrix
        if self.x_hat is None:
            pickle.dump(None, fh)
        else:
            self.x_hat.dump(fh)

    @classmethod
    def load(cls,fh):
        """ This is an example of loading a TelescopeModel
        with open(opts.generate_filename('checkpoint.pickle'),'r') as fh:
            new_tm = TelescopeModel.load(fh)
            print new_tm.rownames[:5]
            print new_tm.colnames[:5]
            print new_tm.shape
            if new_tm.x_hat is None:
                print "x_hat is none"
            else:
                print new_tm.x_hat
        """
        _read_index, _tx_index = pickle.load(fh)
        _Q = csr_matrix.load(fh)

        obj = cls(_read_index, _tx_index, qmat=_Q)

        obj.pi_0 = np.load(fh)
        obj.pi = np.load(fh)
        obj.theta = np.load(fh)

        obj.x_hat = csr_matrix.load(fh)

        return obj

    def matrix_em(self, opts):
        # Propose initial estimates for pi and theta
        R,T    = self.shape
        self.pi    = np.repeat(1./T, T)
        self.theta = np.repeat(1./T, T)

        # weight of each read is the maximum mapping score (np.ndarray, (R,) )
        _weights = self.Q.maxr()

        # total weight for unique reads: sum(weights * Y)
        _u_total  = _weights.multiply(csr_matrix(self.Y[:,None])).sum()
        # total weight for non-unique reads: sum(weights * (1-Y))
        _nu_total = _weights.multiply(csr_matrix(1 - self.Y[:,None])).sum()

        # weight the prior values by the maximum weight overall
        _pi_prior    = opts.piPrior * _weights.max() #max(weights)
        _theta_prior = opts.thetaPrior * _weights.max() #max(weights)

        # pisum0 is the weighted proportion of unique reads assigned to each
        # genome (np.matrix, 1xG)
        _pisum0 = self.Q.multiply(csr_matrix(self.Y[:,None])).sumc()

        for iter_num in xrange(opts.maxIter):
            #--- Expectation step:
            # delta_hat[i,] is the expected value of x[i,] computed using
            # current estimates for pi and theta.
            _numerator = self.Q.multiply(
                csr_matrix( self.pi * self.theta**((1-self.Y)[:,None]))
            )
            self.x_hat = _numerator.normr()
            # w_hat[i,] is the expected value of x[i,] weighted by mapping score
            # (csr_matrix_plus RxG)
            _w_hat = self.x_hat.multiply(_weights)

            #--- Maximization step
            # thetasum is the weighted proportion of non-unique reads assigned
            # to each genome (np.matrix, 1xG)
            _thetasum = _w_hat.multiply(csr_matrix(1 - self.Y[:,None])).sumc()
            # pisum is the weighted proportion of all reads assigned to each genome
            _pisum = _pisum0 + _thetasum

            # Estimate pi_hat
            _pi_denom = _u_total + _nu_total + _pi_prior * T
            _pi_hat = (_pisum + _pi_prior) / _pi_denom

            # Estimate theta_hat
            _theta_denom = _nu_total + _theta_prior * T
            _theta_hat = (_thetasum + _theta_prior) / _theta_denom

            # Difference between pi and pi_hat
            _pidiff = abs(self.pi - _pi_hat).sum()
            if opts.verbose:
                print >>sys.stderr, "[%d]%g" % (iter_num, _pidiff)

            # Set pi_0 if this is the first iteration
            if iter_num == 0: self.pi_0 = _pi_hat.A1

            self.pi     = _pi_hat.A1
            self.theta  = _theta_hat.A1

            # Perform checkpointing
            if opts.checkpoint:
                if iter_num % opts.checkpoint_interval == 0:
                    if opts.verbose: print >>sys.stderr, "Checkpointing... " ,
                    _fn = opts.generate_filename('checkpoint.%03d.p' % iter_num)
                    with open(_fn,'w') as outh:
                        self.dump(outh)
                    if opts.verbose: print >>sys.stderr, "done."

            # Exit if pi difference is less than threshold
            if _pidiff <= opts.emEpsilon:
                break
