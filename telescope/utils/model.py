__author__ = 'bendall'

import sys
import numpy as np
import scipy.sparse
from sparse_matrix import csr_matrix_plus as csr_matrix

try:
    import cPickle as pickle
except ImportError:
    import pickle

from helpers import phred


def reassign_best(mat, method='exclude'):
    """ Reads are reassigned to the transcript with the highest probability
            "method" defines how we deal with reads that have more than one
            best transcript:
                exclude - reads with > 1 best hits are excluded
                average - the value is 1 / num_best_hits
                choose  - one of the best hits is randomly chosen
    """
    v = mat.maxidxr()
    # Check if every read is assigned once
    if np.all(v.sumr() <= 1):
        return v
    else:
        if method == 'exclude':
            return v.multiply(np.where(v.sumr()>1, 0, 1))
        if method == 'choose':
            return mat.maxidxr(choose=True)
        if method == 'average':
            return v.normr()

def reassign_conf(mat, thresh):
    """ Reads are reassigned to transcript if probability > thresh
    """
    f = lambda x: 1 if x >= thresh else 0
    return mat.apply_func(f)

class TelescopeModel:
    '''

    '''
    report_columns = ['transcript',
                      'final_count', 'final_conf', 'final_prop',
                      'init_aligned', 'init_best','init_unique', 'init_prop',]

    report_formats = ['%s',
                      '%d', '%d', '%.6g',
                      '%d', '%d', '%d', '%.6g',]
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

    def calculate_unique_counts(self):
        ''' Calculates number of uniquely mapping reads for each transcript
                - Multiply Q by Y to set values for non-unique reads to zero,
                  then count the number of nonzero values in each column.
        '''
        return self.Q.multiply(csr_matrix(self.Y[:,None])).countc()

    def calculate_fractional_counts(self):
        ''' Calculates the "fractional count" for each transcript
                - Set nonzero values in x_init to 1, then divide by the row
                  total. Fractional counts are the sums of each column.
        '''
        return self.x_init.ceil().normr().sumc().A1

    def make_report(self, conf_prob ,sortby='final_count'):
        '''

        :param conf_prob:
        :param sortby:
        :param other:
        :return:
        '''
        _header = self.report_columns


        report_data = {}
        report_data['transcript'] = self.txnames
        report_data['final_count'] = reassign_best(self.x_hat, method='exclude').sumc().A1
        report_data['final_conf'] = reassign_conf(self.x_hat, thresh=conf_prob).sumc().A1
        report_data['final_prop'] = self.pi

        # Total number of alignments
        report_data['init_aligned'] = self.x_init.ceil().sumc().A1
        # Number of 'best' alignments
        report_data['init_best'] = reassign_best(self.x_init, method='exclude').sumc().A1
        # Number of unambiguous alignments
        report_data['init_unique'] = self.x_init.ceil().multiply(csr_matrix(self.Y[:,None])).sumc().A1
        # First estimate of pi
        report_data['init_prop']  = self.pi_0

        #v = reassign_best(self.Q, method='choose')
        # report_data['init_unique'] = v.multiply(csr_matrix(self.Y[:,None])).sumc().A1

        #
        #
        #report_data['init_best'] = reassign_best(self.x_init, method='exclude').sumc().A1
        #report_data['init_conf'] = reassign_conf(self.x_init, thresh=conf_prob).sumc().A1

        #report_data['unique_counts'] = v.multiply(csr_matrix(self.Y[:,None])).sumc().A1
        #report_data['best_counts'] = v.sumc().A1
        # report_data['fractional_counts'] = self.x_init.ceil().normr().sumc().A1
        # report_data['total_alignments'] = self.x_init.ceil().sumc().A1

        R,T = self.shape
        header = [h for h in _header if h in report_data]
        _unfmt = [[report_data[h][j] for h in header] for j in range(T)]
        _unfmt.sort(key=lambda x:x[header.index(sortby)], reverse=True)
        # Format the rows
        _rows = [[f % v for f,v in zip(self.report_formats,r)] for r in _unfmt]
        return [header] + _rows

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
