# -*- coding: utf-8 -*-
from builtins import object

import os
from nose import with_setup
from nose.tools import assert_equals
from random import randrange

from telescope.tests import TEST_DATA_DIR
from telescope.utils.annotation import get_annotation_class

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


class TestAnnotationIntervalTree(object):

    @classmethod
    def setup_class(cls):
        print("setup_class() before any methods in this class")
        cls.AnnotationClass = get_annotation_class('intervaltree')

    @classmethod
    def teardown_class(cls):
        print("teardown_class() after any methods in this class")
        del cls.AnnotationClass

    def setup(self):
        print("TestAnnotationIntervalTree:setup() before each test method")
        self.gtffile = os.path.join(TEST_DATA_DIR, 'annotation_test.2.gtf')
        self.A = self.AnnotationClass(self.gtffile, 'locus')

    def teardown(self):
        print("TestAnnotationIntervalTree:teardown() after each test method")
        del self.A

    def test_correct_type(self):
        assert type(self.A) is get_annotation_class('intervaltree')

    def test_annot_created(self):
        print(type(self.A))
        assert_equals(self.A.key, 'locus')

    def test_annot_treesize(self):
        assert_equals(len(self.A.itree['chr1']), 3)
        assert_equals(len(self.A.itree['chr2']), 4)
        assert_equals(len(self.A.itree['chr3']), 2)

    def test_empty_lookups(self):
        assert not self.A.intersect_blocks('chr1', [(1, 9999)])
        assert not self.A.intersect_blocks('chr1', [(20001, 39999)])
        assert not self.A.intersect_blocks('chr1', [(50001, 79999)])
        assert not self.A.intersect_blocks('chr1', [(90001, 90001)])
        assert not self.A.intersect_blocks('chr1', [(190000, 590000)])
        assert not self.A.intersect_blocks('chr2', [(1, 9999)])
        assert not self.A.intersect_blocks('chr3', [(1, 9999)])
        assert not self.A.intersect_blocks('chr4', [(1, 1000000000)])
        assert not self.A.intersect_blocks('chrX', [(1, 1000000000)])

    def test_simple_lookups(self):
        lines = (l.strip('\n').split('\t') for l in open(self.gtffile, 'rU'))
        for l in lines:
            iv = (int(l[3]), int(l[4]))
            loc = l[8].split('"')[1]
            r = self.A.intersect_blocks(l[0], [iv])
            assert loc in r
            assert (r[loc] - 1) == (iv[1] - iv[0]), '{} not equal to {}'.format(r[loc], iv[1] - iv[0])

    def test_overlap_lookups(self):
        assert self.A.intersect_blocks('chr1', [(1, 10000)])['locus1'] == 1
        assert self.A.intersect_blocks('chr2', [(1, 10000)])['locus4'] == 1
        assert self.A.intersect_blocks('chr3', [(1, 10000)])['locus7'] == 1
        r = self.A.intersect_blocks('chr1', [(19990, 40000)])
        assert r['locus1'] == 11 and r['locus2'] == 1
        r = self.A.intersect_blocks('chr2', [(44990, 46010)])
        assert r['locus5'] == 22
        r = self.A.intersect_blocks('chr3', [(44990, 46010)])
        assert r['locus8'] == 1021

    def test_subregion_chrom(self):
        sA = self.A.subregion('chr3')
        assert not sA.intersect_blocks('chr1', [(1, 10000)])
        assert not sA.intersect_blocks('chr2', [(1, 10000)])
        assert sA.intersect_blocks('chr3', [(1, 10000)])['locus7'] == 1
        r = sA.intersect_blocks('chr3', [(44990, 46010)])
        assert r['locus8'] == 1021

    def test_subregion_reg(self):
        sA = self.A.subregion('chr3', 30000, 50000)
        assert not sA.intersect_blocks('chr1', [(1, 10000)])
        assert not sA.intersect_blocks('chr2', [(1, 10000)])
        assert not sA.intersect_blocks('chr3', [(1, 10000)])
        assert sA.intersect_blocks('chr3', [(40000, 45000)])['locus8'] == 5001
        r = sA.intersect_blocks('chr3', [(44990, 46010)])
        assert r['locus8'] == 1021
