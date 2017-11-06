from pysam.libcalignedsegment cimport AlignedSegment
from pysam.libcalignedsegment cimport AlignmentFile

from pysam.libchtslib cimport bam_destroy1

cdef class AlignedPair:
    cdef public AlignedSegment r1
    cdef public AlignedSegment r2

    cpdef int write(self, AlignmentFile outfile)
