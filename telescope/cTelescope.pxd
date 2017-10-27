from pysam.libcalignedsegment cimport AlignedSegment

cdef class AlignedPair:
    cdef AlignedSegment r1
    cdef AlignedSegment r2

