# -*- coding: utf-8 -*-
""" Parse SAM/BAM alignment files

"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
from builtins import *

from .helpers import merge_blocks

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2017 Matthew L. Bendall"


class AlignedPair(object):
    def __init__(self, r1, r2 = None):
        self.r1 = r1
        self.r2 = r2
        # set properties
        self.numreads = 1 if r2 is None else 2
        self.is_paired = r2 is not None
        self.is_unmapped = r1.is_unmapped

    def write(self, outfile):
        """ Write AlignedPair to file

        :param outfile:
        :return:
        """
        ret = outfile.write(self.r1)
        if self.r2:
            ret += outfile.write(self.r2)
        return ret

    def set_tag(self, tag, value, value_type=None, replace=True):
        self.r1.set_tag(tag, value, value_type, replace)
        if self.r2:
            self.r2.set_tag(tag, value, value_type, replace)

    def set_mapq(self, value):
        self.r1.mapping_quality = value
        if self.r2:
            self.r2.mapping_quality = value

    def set_flag(self, b):
        self.r1.flag = (self.r1.flag | b)
        if self.r2:
            self.r2.flag = (self.r2.flag | b)
        assert (self.r1.flag & b) == b

    def unset_flag(self, b):
        self.r1.flag = (self.r1.flag ^ (self.r1.flag & b))
        if self.r2:
            self.r2.flag = (self.r2.flag ^ (self.r2.flag & b))
        assert (self.r1.flag & b) == 0

    @property
    def ref_name(self):
        return self.r1.reference_name

    @property
    def query_id(self):
        if self.r2 is None:
                if self.r1.is_read2:
                    return self.r1.query_name + '/2'
                else:
                    return self.r1.query_name + '/1'
        else:
            return self.r1.query_name

    @property
    def refblocks(self):
        if self.r2 is None:
            return merge_blocks(self.r1.get_blocks(), 1)
        else:
            blocks = self.r1.get_blocks() + self.r2.get_blocks()
            return merge_blocks(blocks, 1)
    @property
    def alnlen(self):
        return sum(b[1]-b[0] for b in self.refblocks)

    @property
    def alnscore(self):
        if self.r2 is None:
            return self.r1.get_tag('AS')
        else:
            return self.r1.get_tag('AS') + self.r2.get_tag('AS')
