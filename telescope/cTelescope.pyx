# -*- coding: utf-8 -*-
from telescope.utils.helpers import merge_blocks

cdef class AlignedPair:

    def __init__(self, AlignedSegment r1, AlignedSegment r2 = None):
        self.r1 = r1
        self.r2 = r2

    cpdef int write(self, AlignmentFile outfile):
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

    @property
    def numreads(self):
        return 1 if self.r2 is None else 2

    @property
    def is_paired(self):
        return self.r2 is not None

    @property
    def is_unmapped(self):
       return self.r1.is_unmapped

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


def readkey(aln):
    ''' Key for read '''
    return (aln.query_name, aln.is_read1,
            aln.reference_id, aln.reference_start,
            aln.next_reference_id, aln.next_reference_start)


def matekey(aln):
    ''' Key for mate '''
    return (aln.query_name, not aln.is_read1,
            aln.next_reference_id, aln.next_reference_start,
            aln.reference_id, aln.reference_start)


def pair_alignments(alniter):
    readcache = {}
    for aln in alniter:
        if not aln.is_paired:
            yield AlignedPair(aln)
        else:
            mate = readcache.pop(matekey(aln), None)
            if mate is not None:  # Mate found in cache
                if aln.is_read1:
                    yield AlignedPair(aln, mate)
                else:
                    yield AlignedPair(mate, aln)
            else:  # Mate not found in cache
                readcache[readkey(aln)] = aln
    # Yield the remaining reads in the cache as unpaired
    for aln in readcache.values():
        yield AlignedPair(aln)

def organize_bundle(alns):
    if not alns[0].is_paired:
        assert all(not a.is_paired for a in alns), 'mismatch pair flag'
        return [AlignedPair(a) for a in alns], None
    else:
        if len(alns) == 2 and alns[0].is_unmapped and alns[1].is_unmapped:
            # Unmapped fragment
            return [AlignedPair(alns[0], alns[1])], None
        elif alns[0].is_proper_pair:
            return list(pair_alignments(alns)), None
        else:
            ret1, ret2 = list(), list()
            for aln in alns:
                if aln.is_read1:
                    ret1 += [AlignedPair(aln)]
                else:
                    assert aln.is_read2
                    ret2 += [AlignedPair(aln)]
            return ret1, ret2


def fetch_bundle(samfile, **kwargs):
    """ Iterate over alignment over reads with same ID """
    samiter = samfile.fetch(**kwargs)
    bundle = [ next(samiter) ]
    for aln in samiter:
        if aln.query_name == bundle[0].query_name:
            bundle.append(aln)
        else:
            yield bundle
            bundle = [aln]
    yield bundle


def fetch_fragments(samfile, **kwargs):
    biter = fetch_bundle(samfile, **kwargs)
    for bundle in biter:
        f1, f2 = organize_bundle(bundle)
        yield f1
        if f2 is not None:
            yield f2

