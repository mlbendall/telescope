# -*- coding: utf-8 -*-
""" Provides sparse matrix classes augmented with additional functions
"""
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import range

import numpy as np
import scipy.sparse

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

def _recip0(v):
    ''' Return the reciprocal of a vector '''
    old_settings = np.seterr(divide='ignore')
    ret = 1. / v
    ret[np.isinf(ret)] = 0
    np.seterr(**old_settings)
    return ret

class csr_matrix_plus(scipy.sparse.csr_matrix):

    def norm(self, axis=None):
        """ Normalize matrix along axis

        Args:
            axis:

        Returns:
        Examples:
            >>> row = np.array([0, 0, 1, 2, 2, 2])
            >>> col = np.array([0, 2, 2, 0, 1, 2])
            >>> data = np.array([1, 2, 3, 4, 5, 6])
            >>> M = csr_matrix_plus((data, (row, col)), shape=(3, 3))
            >>> print(M.norm(1).toarray())
            [[ 0.33333333  0.          0.66666667]
             [ 0.          0.          1.        ]
             [ 0.26666667  0.33333333  0.4       ]]
        """
        # return self._norm_loop(axis)
        return self._norm(axis)

    def _norm(self, axis=None):
        if axis is None:
            return type(self)(self.multiply(1. / self.sum()))
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            return type(self)(self.multiply(_recip0(self.sum(1))))

    def _norm_loop(self, axis=None):
        if axis is None:
            ret = self.copy().astype(np.float)
            ret.data /= sum(ret)
            return ret
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            ret = self.copy().astype(np.float)
            rowiter = zip(ret.indptr[:-1], ret.indptr[1:], ret.sum(1).A1)
            for d_start, d_end, d_sum in rowiter:
                if d_sum != 0:
                    ret.data[d_start:d_end] /= d_sum
            return ret

    def scale(self, axis=None):
        """ Scale matrix so values are between 0 and 1

        Args:
            axis:

        Returns:
        Examples:
            >>> M = csr_matrix_plus([[10, 0, 20],[0, 0, 30],[40, 50, 60]])
            >>> print(M.scale().toarray())
            [[ 0.1  0.   0.2]
             [ 0.   0.   0.3]
             [ 0.4  0.5  1. ]]
            >>> print(M.scale(1).toarray())
            [[ 0.5  0.   1. ]
             [ 0.   0.   1. ]
             [ 0.4  0.5  1. ]]
        """
        return self._scale(axis)

    def _scale(self, axis=None):
        if axis is None:
            return type(self)(self.multiply(1. / self.max()))
            # ret = self.copy().astype(np.float)
            # return ret / ret.max()
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            return type(self)(self.multiply(_recip0(self.max(1).toarray())))

    def binmax(self, axis=None):
        """ Set max values to 1 and others to 0

        Args:
            axis:

        Returns:
        Examples:
            >>> M = csr_matrix_plus([[6, 0, 2],[0, 0, 3],[4, 5, 6]])
            >>> print(M.binmax().toarray())
            [[1 0 0]
             [0 0 0]
             [0 0 1]]
            >>> print(M.binmax(1).toarray())
            [[1 0 0]
             [0 0 1]
             [0 0 1]]
        """
        if axis is None:
            raise NotImplementedError
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            _data = np.zeros(self.data.shape, dtype=np.int8)
            rit = zip(self.indptr[:-1], self.indptr[1:], self.max(1).toarray())
            for d_start, d_end, d_max in rit:
                _data[d_start:d_end] = (self.data[d_start:d_end] == d_max)
            ret = type(self)((_data, self.indices.copy(), self.indptr.copy()),
                              shape=self.shape)
            ret.eliminate_zeros()
            return ret

    def count(self, axis=None):
        if axis is None:
            raise NotImplementedError
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            ret = self.indptr[1:] - self.indptr[:-1]
            return np.array(ret, ndmin=2).T

    def choose_random(self, axis=None):
        if axis is None:
            raise NotImplementedError
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            ret = self.copy()
            for d_start, d_end in zip(ret.indptr[:-1], ret.indptr[1:]):
                if d_end - d_start > 1:
                    chosen = np.random.choice(range(d_start, d_end))
                    for j in range(d_start, d_end):
                        if j != chosen:
                            ret.data[j] = 0
            ret.eliminate_zeros()
            return ret

    def check_equal(self, other):
        if self.shape != other.shape:
            return False
        return (self != other).nnz == 0

    def apply_func(self, func):
        ret = self.copy()
        ret.data = np.fromiter((func(v) for v in self.data),
                               self.data.dtype, count=len(self.data))
        return ret

    def save(self, filename):
        np.savez(filename, data=self.data, indices=self.indices,
                 indptr=self.indptr, shape=self.shape)
    @classmethod
    def load(cls, filename):
        loader = np.load(filename)
        return cls((loader['data'], loader['indices'], loader['indptr']),
                   shape = loader['shape'])
