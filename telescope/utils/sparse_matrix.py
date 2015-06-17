__author__ = 'bendall'

import math
import scipy.sparse
import numpy as np

class csr_matrix_plus(scipy.sparse.csr_matrix):

    def __init__(self, *args, **kwargs):
        super(csr_matrix_plus, self).__init__(*args, **kwargs)

    def sumr(self):
        ''' Sum across rows '''
        return self.sum(1)

    def sumc(self):
        ''' Sum across columns '''
        return self.sum(0)

    def maxr(self):
        ''' Returns max value for each row as Rx1 sparse matrix '''
        return type(self)(self.max(1))

    def maxc(self):
        ''' Returns max value for each column as 1xC sparse matrix '''
        return type(self)(self.max(0))

    def argmaxr(self):
        ''' Returns index for max value in each row as numpy array '''
        return np.argmax(self.toarray(),1)

    def argmaxc(self):
        ''' Returns index for max value in each column as numpy array '''
        return np.argmax(self.toarray(),0)

    def countr(self):
        ''' Returns count of nonzero elements in each row as numpy array '''
        return np.bincount(self.nonzero()[0], minlength=self.shape[0])

    def countc(self):
        ''' Returns count of nonzero elements in each column as numpy array '''
        return np.bincount(self.nonzero()[1], minlength=self.shape[1])

    def exp(self):
        ''' Element-wise exp(x), which is the natural exponential function e^x'''
        return type(self)((np.exp(self.data), self.indices, self.indptr), shape=self.shape)

    def expm1(self):
        ''' Element-wise exp(x)-1, provides greater precision '''
        return super(type(self),self).expm1()

    def ceil(self):
        ''' Element-wise '''
        return type(self)((np.ceil(self.data), self.indices, self.indptr), shape=self.shape)

    def normr(self):
        ''' Return copy of matrix normalized by row (row sums == 1)
        :param m: Matrix to be normalized
        :return:  Copy of matrix, normalized
        '''
        assert scipy.sparse.isspmatrix_csr(self), "Matrix is not sparse CSR matrix"
        return type(self)( self.multiply(1./self.sum(1)) )

    def multiply(self,other):
        return type(self)(super(type(self),self).multiply(other))

    def apply_func(self, func):
        vfunc = np.vectorize(func)
        return type(self)((vfunc(self.data), self.indices, self.indptr), shape=self.shape)

    def maxidxr(self):
        newdata = np.copy(self.data)
        for n in xrange(self.shape[0]):
          rowvals = self.data[self.indptr[n]:self.indptr[n+1]]
          newdata[self.indptr[n]:self.indptr[n+1]] = np.where(rowvals==np.max(rowvals), 1, 0)
        return type(self)((newdata, self.indices, self.indptr), shape=self.shape)
    
    def pretty_tsv(self, rownames, colnames):
        ret = [ '\t'.join([''] + colnames) ]
        for i,rn in enumerate(rownames):
            vals = self.getrow(i).toarray()[0]
            ret.append('%s\t%s' % (rn, '\t'.join('%.8g' % f for f in vals)))
        return '\n'.join(ret)

class TelescopeMatrix:
  from copy import deepcopy

  def __init__(self, data, row_idx, col_idx):
    self.ridx = row_idx
    self.rownames = [k for k,v in sorted(self.ridx.iteritems(),key=lambda x:x[1])]
    self.cidx = col_idx
    self.colnames = [k for k,v in sorted(self.cidx.iteritems(),key=lambda x:x[1])]
    self.shape = (len(self.ridx),len(self.cidx))

    self.rows = [dict() for _ in self.ridx]
    for i,j,v in data:
      if j in self.rows[i]:
        print "WARNING: already have value"
      self.rows[i][j] = v

  def _func_row(self,func):
    ret = []
    for d in self.rows:
      vals = d.values()
      vals += [] if len(d) == self.shape[1] else [0]
      ret.append(func(vals))
    return ret

  def count(self):
    return self._func_row(lambda x:sum(v>0 for v in x))

  def max(self):
    return self._func_row(max)

  def min(self):
    return self._func_row(min)

  def sum(self):
    return self._func_row(sum)

  def copy(self):
    return deepcopy(self)

  """
  def mult_row(self, v, inplace=False):
    ''' Multiply each row by corresponding value in v '''
    other = self if inplace else deepcopy(self)
    assert isinstance(v, list), "v must be a list"
    assert len(v) == self.shape[0], "Shape is incorrect %d" % len(v)
    for i,r in enumerate(other.rows):
      for c in r.keys():
        r[c] *= v[i]
    return other

  def mult_item(self, v, inplace=False):
    ''' Multiply each item by v '''
    other = self if inplace else deepcopy(self)
    assert isinstance(v, (int, float)), "v must be int or float"
    for r in other.rows:
      for c in r.keys():
        r[c] *= v
    return other
  """
  def multiply(self, v, inplace=False):
    ''' Determine whether to multiply item or row based on v '''
    other = self if inplace else deepcopy(self)

    if type(v) is list:
      assert len(v) == self.shape[0], "Shape is incorrect %d" % len(v)
      for i,r in enumerate(other.rows):
        for c in r.keys():
          r[c] *= v[i]
    else:
      for r in other.rows:
        for c in r.keys():
          r[c] *= v
    return other

  def add(self, v, inplace=False):
    ''' Determine whether to multiply item or row based on v '''
    other = self if inplace else deepcopy(self)

    if type(v) is list:
      assert len(v) == self.shape[0], "Shape is incorrect %d" % len(v)
      for i,r in enumerate(other.rows):
        for c in r.keys():
          r[c] += v[i]
    else:
      for r in other.rows:
        for c in r.keys():
          r[c] += v
    return other

  def exp(self, inplace=False):
    other = self if inplace else deepcopy(self)
    for r in other.rows:
      for c in r.keys():
        r[c] = math.exp(r[c])
    return other

  def expm1(self, inplace=False):
    other = self if inplace else deepcopy(self)
    for r in other.rows:
      for c in r.keys():
        r[c] = math.expm1(r[c])
    return other

  def normalize(self, inplace=False):
    ''' Normalize matrix by dividing elements in each row by the row sum'''
    other = self if inplace else deepcopy(self)
    _inv = [1./v if v != 0 else 1 for v in other.sum()]
    return other.multiply(_inv, inplace)

  def rescale(self, inplace=False):
    ''' Rescale matrix '''
    other = self if inplace else deepcopy(self)
    # Get max and min values of matrix
    _mmax = max(other.max())
    _mmin = min(other.min())
    if _mmin < 0:
      _scaling_factor = 100.0 / (_mmax - _mmin)
      other = other.add(abs(_mmin), inplace)
    else:
      _scaling_factor = 100.0 / _mmax
    return other.multiply(_scaling_factor).exp()

  def tocounts(self, inplace=False):
    '''
    :param inplace:
    :return:
    '''
    other = self if inplace else deepcopy(self)
    for r in other.rows:
      for c in r.keys():
        if r[c] != 0: r[c] = 1
    return other

  def col_sum(self):
    ret = [0] * self.shape[1]
    for r in self.rows:
      for j,v in r.iteritems():
        ret[j] += v
    return ret

  def col_count(self):
    ret = [0] * self.shape[1]
    for r in self.rows:
      for j,v in r.iteritems():
        ret[j] += 1
    return ret

  def __unicode__(self):
    ret = u'\t%s\n' % u'\t'.join(self.colnames)
    for r,row in enumerate(self.rows):
      ret += u'%s\t' % self.rownames[r]
      ret += u'%s\n' % '\t'.join(str(row[c]) if c in row else '0' for c in range(len(self.colnames)))
    return ret

  def __str__(self):
    return unicode(self).encode('utf-8')

