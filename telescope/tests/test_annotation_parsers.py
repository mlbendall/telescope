__author__ = 'bendall'

import os

# Import path to test data
from telescope.tests import TEST_DATA_DIR

# Import classes for testing
from telescope.utils.annotation_parsers import _AnnotationBisect

class TestAnnotationBisect:

    @classmethod
    def setup_class(cls):
        _gtffile = os.path.join(TEST_DATA_DIR,'annotation_test.gtf')
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