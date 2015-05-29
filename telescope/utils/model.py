__author__ = 'bendall'

import sys
import numpy as np
import scipy.sparse
from sparse_matrix import csr_matrix_plus as csr_matrix

from helpers import phred

def reassign_best(mat):
    """ Reads are reassigned to the transcript with the highest probability
    """
    return mat.maxidxr().normr()

def reassign_conf(mat, thresh=0.99):
    """ Reads are reassigned to transcript if probability > thresh
    """
    f = lambda x: 1 if x >= thresh else 0
    return mat.apply_func(f)

class TelescopeModel:
    def __init__(self, data, row_idx, col_idx):
        self.ridx = row_idx
        self.rownames = [k for k,v in sorted(self.ridx.iteritems(),key=lambda x:x[1])]
        self.cidx = col_idx
        self.colnames = [k for k,v in sorted(self.cidx.iteritems(),key=lambda x:x[1])]
        self.shape = (len(self.ridx),len(self.cidx))

        i,j,d = zip(*data)
        coo = scipy.sparse.coo_matrix((d,(i,j)),shape=self.shape)
        raw_scores = csr_matrix(coo)
        self.Q = raw_scores.multiply(100.0 / raw_scores.max()).exp()
        self.x_init = self.Q.normr()
        self.Y = np.where(self.Q.countr()==1, 1, 0)

        self.pi_0  = None
        self.pi    = None
        self.theta = None
        self.x_hat = None

    def calculate_unique_counts(self):
        # print self.Q.multiply(csr_matrix(self.Y[:,None]))
        return self.Q.multiply(csr_matrix(self.Y[:,None])).countc()

    def calculate_fractional_counts(self):
        return self.x_init.ceil().normr().sumc().A1

    def calculate_weighted_counts(self):
        return self.Q.normr().sumc().A1

    def make_report(self, conf_prob ,sortby='final_best'):
        header = ['genome', 'final_best', 'final_conf', 'final_prop',
                  'init_best', 'init_conf', 'init_prop',
                  'unique_counts', 'weighted_counts','fractional_counts',
                  ]
        report_data = {}
        report_data['genome']   = self.colnames

        report_data['final_best'] = reassign_best(self.x_hat).sumc().A1
        report_data['final_conf'] = reassign_conf(self.x_hat, thresh=conf_prob).sumc().A1
        report_data['final_prop'] = self.pi

        report_data['init_best'] = reassign_best(self.x_init).sumc().A1
        report_data['init_conf'] = reassign_conf(self.x_init, thresh=conf_prob).sumc().A1
        report_data['init_prop']  = self.pi_0

        report_data['unique_counts'] = self.calculate_unique_counts()
        report_data['weighted_counts'] = self.calculate_weighted_counts()
        report_data['fractional_counts'] = self.calculate_fractional_counts()

        R,G = self.shape
        comment = ['# Aligned reads:', str(R), 'Genomes', str(G)]
        header = [h for h in header if h in report_data]
        _rows = [[report_data[h][j] for h in header] for j in range(G)]
        _rows.sort(key=lambda x:x[header.index(sortby)], reverse=True)
        return [comment, header] + _rows

'''
    def make_report(self, l1_thresh=0.5, l2_thresh=0.01):
        # assert ngenomes==len(rdata['final_pi'])
        comment = ['# Aligned reads:', self.shape[0], 'Genomes', self.shape[1]]
        header = ['genome','final_pi','final_best_hit', 'final_best','final_l1','final_l2',
                  'init_pi', 'init_best_hit',  'init_best', 'init_l1', 'init_l2',
                  'fractional_counts','weighted_counts','weighted_raw_counts','unique_counts']
        for row in self.rownames:

  _rows = [[rdata[h][j] for h in header] for j in range(ngenomes)]
  _rows.sort(key=lambda x:x[header.index(sortby)], reverse=True)
  return [comment, header] + _rows
        ret = ''
        mat = self.x_init
        row_max = mat.maxr()
        arr = mat.toarray()
        # maxhits = np.where(mat_exp==row_max)
        # print maxhits
        maxhits = np.where(arr==row_max, 1, 0)
        nmaxhits = maxhits / maxhits.sum(1)[:, np.newaxis]
        best_reads = nmaxhits.sum(0)
        print best_reads

        l1_reads = np.where(arr >= l1_thresh, 1, 0)
        level1 = (1. * l1_reads.sum(0)) / self.shape[0]
        print level1.sum()

        l2_reads = np.where((arr < l1_thresh) & (arr > l2_thresh), 1, 0)
        level2 = (1. * l2_reads.sum(0)) / self.shape[0]
        print level2.sum()

        #l1_reads = np.where(arr >= l1_thresh, 1, 0)
        #level1 = l1_reads.sum(0)

        # num_best = np.unique(maxhits[0].A1, return_counts=True)[1]
        # print num_best
'''

def matrix_em(Q, opts):
  # Initialize model
  R,G   = Q.shape
  pi    = np.repeat(1./G, G)
  theta = np.repeat(1./G, G)

  # Y is the uniqueness indicator. Y_i==1 if read is unique, 0 otherwise (np.ndarray, (R,) )
  Y     = np.where(Q.countr()==1, 1, 0)
  # Y_mat = csr_matrix(np.where(Q.countr()==1, 1, 0)[:,None])

  # weights is the weight of each read (np.ndarray, (R,) )
  weights = Q.maxr() #.todense().A1
  # type(weights)

  # u_total  = sum(weights * Y)
  u_total  = weights.multiply(csr_matrix(Y[:,None])).sum()
  # nu_total = sum(weights * (1-Y))
  nu_total = weights.multiply(csr_matrix(1-Y[:,None])).sum()

  # Prior values weighted by max weight
  pi_prior    = opts.piPrior * weights.max() #max(weights)
  theta_prior = opts.thetaPrior * weights.max() #max(weights)

  # pisum0 is the weighted proportion of unique reads assigned to each genome (np.matrix, 1xG)
  pisum0 = Q.multiply(csr_matrix(Y[:,None])).sumc()

  for iter_num in xrange(opts.maxIter):
    #--- Expectation step:
    # q_hat is the numerator in the expected values of x_i
    # (csr_matrix_plus RxG)
    q_hat = Q.multiply( csr_matrix( pi * theta**((1-Y)[:,None])) )
    delta_hat = q_hat.normr()
    # (csr_matrix_plus RxG)
    w_hat = q_hat.normr().multiply(weights)

    #--- Maximization step
    # thetasum is the weighted proportion of non-unique reads assigned to each genome (np.matrix, 1xG)
    thetasum = w_hat.multiply(csr_matrix(1-Y[:,None])).sumc()
    # pisum is the weighted proportion of all reads assigned to each genome
    pisum = pisum0 + thetasum

    # Estimate pi_hat
    pi_denom = u_total + nu_total + pi_prior * G
    pi_hat = (pisum + pi_prior) / pi_denom

    # Estimate theta_hat
    theta_denom = nu_total + theta_prior * G
    theta_hat = (thetasum + theta_prior) / theta_denom

    cutoff = abs(pi - pi_hat).sum()
    if opts.verbose:
      print >>sys.stderr, "[%d]%g" % (iter_num, cutoff)
    if iter_num == 0: pi_0 = pi_hat.A1
    pi, theta = pi_hat.A1, theta_hat.A1
    if cutoff <= opts.emEpsilon:
      break

  return pi_0, pi, theta, delta_hat