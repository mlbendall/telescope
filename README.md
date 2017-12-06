Telescope
========

###### *Single locus resolution of* **T***ransposable* **ELE***ment expression.*

Affiliations:

+ [Computational Biology Institute](http://cbi.gwu.edu)
+ [Department of Microbiology, Immunology, & Tropical Medicine](https://smhs.gwu.edu/microbiology/)

## Setup

Telescope depends on numpy (>=1.13), scipy (>=0.19.0), intervaltree, and pysam (>=0.12).
If the dependencies are not met they will be automatically installed. 

```bash
pip install git+git://github.com/mlbendall/telescope.git
```

## `telescope assign`

The `telescope assign` program finds overlapping reads between an alignment
(SAM/BAM) and an annotation (GTF) then reassigns reads using a statistical
model.

## `telescope resume`

The `telescope resume` program loads the checkpoint from a previous run and 
reassigns reads using a statistical model.

## Appendix

### SAM tags used by telescope

+ `ZF:Z` Assigned Feature - The name of the feature that alignment is assigned to.
+ `ZT:Z` Telescope tag - A value of `PRI` indicates that this alignment is the
     best hit for the feature and is used in the likelihood calculations. 
     Otherwise the value will be `SEC`, meaning that another alignment to the
     same feature has a higher score.
+ `ZB:Z` Best Feature = The name(s) of the highest scoring feature(s) for the fragment.          
+ `YC:Z` Specifies color for alignment as R,G,B.
UCSC sanctioned tag, see documentation
[here.](http://genome.ucsc.edu/goldenpath/help/hgBamTrackHelp.html)
+ `XP:Z` Alignment probability - estimated posterior probability for this alignment.

### Version History

#### v0.5.4.1
  + Fixed bug where random seed is out of range
  
#### v0.5.4
  + Added MIT license
  + Changes to logging/reporting
  
#### v0.5.3
  + Improvements to `telescope resume`

#### v0.5.2
  + Implemented checkpoint and `telescope resume`

#### v0.5.1
  + Refactoring Telescope class with TelescopeLikelihood
  + Improved memory usage

#### v0.4.2
  +  Subcommand option parsing class
  +  Cython for alignment parsing
  +  HTSeq as alternate alignment parser

#### v0.3.2
  +  Python3 compatibility  

#### v0.3
  + Implemented IntervalTree for Annotation data structure
  + Added support for annotation files where a locus may be non-contiguous.
  + Overlapping annotations with the same key value (locus) are merged
  + User can set minimum overlap criteria for assigning read to locus, default = 0.1

#### v0.2

  + Implemented checkpointing
  + Output tables as pickled objects
  + Changes to report format and output  