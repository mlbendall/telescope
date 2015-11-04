__author__ = 'bendall'

import os

# Import path to test data
from telescope.tests import TEST_DATA_DIR

# Import classes for testing
from telescope.utils.annotation_parsers import _AnnotationBisect

class TestAnnotationBisect:

    @classmethod
    def setup_class(cls):
        _gtffile = os.path.join(TEST_DATA_DIR,'annotation_test.1.gtf')
        cls.annotation = _AnnotationBisect(_gtffile)

    @classmethod
    def teardown_class(cls):
        del cls.annotation

    def test_create_annotation(self):
        assert self.annotation.key == 'locus'
        assert self.annotation._locus_lookup['locus1'][0] == ('chr1', 10000, 20000)

    def test_lookup_chr1_9999(self):
        assert self.annotation.lookup('chr1', 9999) is None

    def test_lookup_chr1_10000(self):
        assert self.annotation.lookup('chr1', 10000) == 'locus1'

    def test_lookup_chr1_20001(self):
        assert self.annotation.lookup('chr1', 15000) == 'locus1'

    def test_lookup_chr1_20000(self):
        assert self.annotation.lookup('chr1', 20000) == 'locus1'

    def test_lookup_chr1_20001(self):
        assert self.annotation.lookup('chr1', 20001) is None

    def test_lookup_chr2_39999(self):
        assert self.annotation.lookup('chr2', 39999) is None

    def test_lookup_chr2_40000(self):
        assert self.annotation.lookup('chr2', 40000) == 'locus5'

    def test_lookup_chr2_45000(self):
        assert self.annotation.lookup('chr2', 45000) == 'locus5'

    def test_lookup_chr2_45001(self):
        assert self.annotation.lookup('chr2', 45001) is None

    def test_lookup_chr2_45500(self):
        assert self.annotation.lookup('chr2', 45500) is None

    def test_lookup_chr2_45999(self):
        assert self.annotation.lookup('chr2', 45999) is None

    def test_lookup_chr2_46000(self):
        assert self.annotation.lookup('chr2', 46000) == 'locus5'

    def test_lookup_chr2_48000(self):
        assert self.annotation.lookup('chr2', 48000) == 'locus5'

    def test_lookup_chr2_51000(self):
        assert self.annotation.lookup('chr2', 51000) == 'locus5'

    def test_lookup_chr2_51001(self):
        assert self.annotation.lookup('chr2', 51001) is None

    # Test interval lookup


# Import classes for testing
from telescope.utils.annotation_parsers import _AnnotationIntervalTree
from intervaltree import Interval, IntervalTree

class TestAnnotationIntervalTree:

    @classmethod
    def setup_class(cls):
        _gtffile = os.path.join(TEST_DATA_DIR, 'annotation_test.2.gtf')
        cls.annotation = _AnnotationIntervalTree(_gtffile, )

    @classmethod
    def teardown_class(cls):
        del cls.annotation

    def test_create_annotation(self):
        assert self.annotation.key == 'locus'
        assert isinstance(self.annotation.itree['chr1'], IntervalTree)

    def test_tree_size(self):
        assert len(self.annotation.itree['chr3']) == 2

    def test_lookup_chr1_9999(self):
        assert self.annotation.lookup('chr1', 9999) is None

    def test_lookup_chr1_10000(self):
        assert self.annotation.lookup('chr1', 10000) == 'locus1'

    def test_lookup_chr1_20001(self):
        assert self.annotation.lookup('chr1', 15000) == 'locus1'

    def test_lookup_chr1_20000(self):
        assert self.annotation.lookup('chr1', 20000) == 'locus1'

    def test_lookup_chr1_20001(self):
        assert self.annotation.lookup('chr1', 20001) is None

    def test_lookup_chr2_39999(self):
        assert self.annotation.lookup('chr2', 39999) is None

    def test_lookup_chr2_40000(self):
        assert self.annotation.lookup('chr2', 40000) == 'locus5'

    def test_lookup_chr2_45000(self):
        assert self.annotation.lookup('chr2', 45000) == 'locus5'

    def test_lookup_chr2_45001(self):
        assert self.annotation.lookup('chr2', 45001) is None

    def test_lookup_chr2_45500(self):
        assert self.annotation.lookup('chr2', 45500) is None

    def test_lookup_chr2_45999(self):
        assert self.annotation.lookup('chr2', 45999) is None

    def test_lookup_chr2_46000(self):
        assert self.annotation.lookup('chr2', 46000) == 'locus5'

    def test_lookup_chr2_48000(self):
        assert self.annotation.lookup('chr2', 48000) == 'locus5'

    def test_lookup_chr2_51000(self):
        assert self.annotation.lookup('chr2', 51000) == 'locus5'

    def test_lookup_chr2_51001(self):
        assert self.annotation.lookup('chr2', 51001) is None

    def test_lookup_interval_chr3_100_260(self):
        # Interval is completely outside
        assert self.annotation.lookup_interval('chr3', 100, 260) is None

    def test_lookup_interval_chr3_10100_10260(self):
        # Interval is completely inside
        assert self.annotation.lookup_interval('chr3', 10100, 10260) == 'locus7'

    def test_lookup_interval_chr3_19900_20100(self):
        # Interval overlaps 50%
        assert self.annotation.lookup_interval('chr3', 19900, 20100) == 'locus7'

    def test_lookup_interval_chr3_19995_20195(self):
        # Interval overlap is not sufficent
        assert self.annotation.lookup_interval('chr3', 19995, 20195) is None

    def test_lookup_interval_chr3_44893_20195(self):
        # Interval overlaps two annotations from same locus
        # These should be merged
        assert self.annotation.lookup_interval('chr3', 44893, 45093) == 'locus8'