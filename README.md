Telescope (private repo)
========

######*Single locus resolution of* **T***ransposable* **ELE***ment expression.*

Affiliations:

+ [Computational Biology Institute](http://cbi.gwu.edu)
+ [Department of Microbiology, Immunology, & Tropical Medicine](https://smhs.gwu.edu/microbiology/)

### SAM tags used by telescope

#### telescope tag

+ `YC:Z` Specifies color for alignment as R,G,B.
UCSC sanctioned tag, see documentation
[here.](http://genome.ucsc.edu/goldenpath/help/hgBamTrackHelp.html)

+ `XC:i` Alignment Count - total number of alignments for this read.

+ `XF:Z` Alignment Feature - name of the feature containing this alignment

+ `ZC:i` Best Count - number of "best" alignments sharing the same maximum score.

+ `ZS:i` Best Score - maximum alignment score for this read.

+ `ZF:Z` Best Feature - name(s) of the feature(s) containing the best alignment(s)

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

+ `XF:Z` Alignment Feature - name of the feature containing this alignment

+ `ZF:Z` Best Feature - name(s) of the feature(s) containing the best alignment(s)

#### telescope load

Load a telescope checkpoint or output matrix and create tables with parameter values. Possible
parameters for `--outparam` include:

  + **pi** - Final genome proportion
  + **pi_0** - Initial genome proportion (before EM)
  + **theta** - Reassignment parameter
  + **x_hat** - Genome indicator, probability of read originating from genome i
  + **x_init** - Initial genome indicator (before EM)
  + **Q** - Rescaled, normalized mapping scores


