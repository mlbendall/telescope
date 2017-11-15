# -*- coding: utf-8 -*-

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"

from tempfile import TemporaryFile

from telescope.utils.sparse_plus import csr_matrix_plus

def sparse_equal(m1, m2):
    if m1.shape != m2.shape:
        return False
    return (m1!=m2).nnz == 0

def test_mplus_identity():
    m1 = csr_matrix_plus([[1, 0, 2],[0, 0, 3],[4, 5, 6 ]])
    assert m1[0, 0] == 1
    assert m1[0, 2] == 2
    assert m1[1, 2] == 3
    assert m1[2, 0] == 4
    assert m1[2, 1] == 5
    assert m1[2, 2] == 6

def test_mplus_norm():
    m1     = csr_matrix_plus([[        1,        0,        2],
                              [        0,        0,        3],
                              [        4,        5,        6]]
                             )
    a_none = csr_matrix_plus([[  (1./21),        0,  (2./21)],
                              [        0,        0,  (3./21)],
                              [  (4./21),  (5./21),  (6./21)]]
                             )
    assert sparse_equal(m1.norm(), a_none)

def test_mplus_norm_row():
    m1     = csr_matrix_plus([[        1,        0,        2],
                              [        0,        0,        3],
                              [        4,        5,        6]]
                             )
    a_row  = csr_matrix_plus([[   (1./3),        0,   (2./3)],
                              [        0,        0,       1.],
                              [  (4./15),  (5./15),  (6./15)]]
                             )
    assert sparse_equal(m1.norm(1), a_row), '\n{}\n{}'.format(m1.norm(1), a_row)

def test_mplus_norm_row_withzero():
    m1     = csr_matrix_plus([[        1,        0,        2],
                              [        0,        0,        0],
                              [        4,        5,        6]]
                             )
    a_row  = csr_matrix_plus([[   (1./3),        0,   (2./3)],
                              [        0,        0,        0],
                              [  (4./15),  (5./15),  (6./15)]]
                             )
    assert sparse_equal(m1.norm(1), a_row), '\n{}\n{}'.format(m1.norm(1), a_row)

def test_mplus_save_load():
    m1     = csr_matrix_plus([[        1,        0,        2],
                              [        0,        0,        3],
                              [        4,        5,        6]]
                             )
    outfile = TemporaryFile()
    m1.save(outfile)
    outfile.seek(0)
    m2 = csr_matrix_plus.load(outfile)
    assert sparse_equal(m1, m2)
