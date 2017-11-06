# -*- coding: utf-8 -*-
""" Provides sparse matrix classes augmented with additional functions
"""
from __future__ import division

from future import standard_library
standard_library.install_aliases()
from builtins import range

try:
    import pickle as pickle
except ImportError:
    import pickle

import numpy as np
import scipy.sparse

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

class csr_matrix_plus(scipy.sparse.csr_matrix):

    def oldnorm(self, axis=None):
        _data = np.array(self.data, dtype=float, copy=True)
        if axis is None:
            raise NotImplementedError
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            _rowmax = self.max(1).data
            for i in range(self.shape[0]):
                _data[self.indptr[i]:self.indptr[i + 1]] /= _rowmax[i]

        return type(self)((_data, self.indices.copy(), self.indptr.copy()))

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
        if axis is None:
            ret = self.copy().astype(np.float)
            return ret / ret.sum()
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
        if axis is None:
            ret = self.copy().astype(np.float)
            return ret / ret.max()
        elif axis == 0:
            raise NotImplementedError
        elif axis == 1:
            ret = self.copy().astype(np.float)
            rowiter = zip(ret.indptr[:-1], ret.indptr[1:], ret.max(1).data)
            for d_start, d_end, d_max in rowiter:
                ret.data[d_start:d_end] /= d_max
            return ret

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
        ret = self.scale(axis).floor().astype(np.uint8)
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

    def checkequal(self, other):
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
