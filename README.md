Telescope
========

######*Single locus resolution of* **T***ransposable* **ELE***ment expression.*

Affiliations:

+ [Computational Biology Institute](http://cbi.gwu.edu)
+ [Department of Microbiology, Immunology, & Tropical Medicine](https://smhs.gwu.edu/microbiology/)

## Setup

Telescope depends on numpy (>=1.7), scipy (>=0.17.0), intervaltree, and pysam (>=0.8.2.1).
If the dependencies are not met they will be automatically installed. 

```bash
pip install git+git://github.com/mlbendall/telescope.git
```

### SAM tags used by telescope

#### telescope tag

+ `YC:Z` Specifies color for alignment as R,G,B.
UCSC sanctioned tag, see documentation
[here.](http://genome.ucsc.edu/goldenpath/help/hgBamTrackHelp.html)

+ `XC:i` Alignment Count - total number of alignments for this read.

+ `XT:Z` Alignment transcript - name of the transcript containing this alignment

+ `ZC:i` Best Count - number of "best" alignments sharing the same maximum score.

+ `ZS:i` Best Score - maximum alignment score for this read.

+ `ZT:Z` Best transcript - name(s) of the transcript(s) containing the best alignment(s)

#### telescope id

+ `YC:Z` Specifies color for alignment as R,G,B.
UCSC sanctioned tag, see documentation
[here.](http://genome.ucsc.edu/goldenpath/help/hgBamTrackHelp.html)

  + **vermilion** - best location, high confidence read (default is greater than 90%) 
  + **yellow** - best location, low confidence read
  + **teal** - best location, but read assigned to multiple genomes
  + **light green** -  not the best location
  + **white** - alternate, the read aligned multiple times to the same feature

+ `XP:Z` Alignment probability - estimated probability for this alignment

+ `XT:Z` Alignment transcript - name of the transcript containing this alignment

+ `ZT:Z` Best transcript - name(s) of the transcript(s) containing the best alignment(s)

#### telescope load

Load a telescope checkpoint or output matrix and create tables with parameter values. Possible
parameters for `--outparam` include:

  + **pi** - Final genome proportion
  + **pi_0** - Initial genome proportion (before EM)
  + **theta** - Reassignment parameter
  + **x_hat** - Genome indicator, probability of read originating from genome i
  + **x_init** - Initial genome indicator (before EM)
  + **Q** - Rescaled, normalized mapping scores

### Version History

#### v0.3

  + Implemented IntervalTree for Annotation data structure
  + Added support for annotation files where a locus may be non-contiguous.
  + Overlapping annotations with the same key value (locus) are merged
  + User can set minimum overlap criteria for assigning read to locus, default = 0.1

#### v0.2

  + Implemented checkpointing
  + Output tables as pickled objects
  + Changes to report format and output  