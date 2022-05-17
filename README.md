Stellarscope [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/telescope/README.html)
========

###### *Single locus resolution of* **T***ransposable* **ELE***ment expression.*

**Affiliations:**

+ [Computational Biology Institute](http://cbi.gwu.edu) at George Washington University
+ [Weill Cornell Medicine Division of Infectious Diseases](https://medicine.weill.cornell.edu/divisions-programs/infectious-diseases)

**Table of Contents:**

* [Installation](#installation)
* [Usage](#usage)
  * [`stellarscope assign`](#stellarscope-assign)
  * [`stellarscope cellsort`](#stellarscope-cellsort)
  * [`stellarscope merge`](#stellarscope-merge)
  * [`stellarscope resume`](#stellarscope-resume)
* [Output](#Output)
  * [Telescope report](#telescope-report)
  * [Updated SAM file](#updated-sam-file)
  * [Feature list (stellarscope only)](#feature-list)
  * [Barcode list (stellarscope only)](#barcode-list)
* [Version History](#version-history)

## Installation

**Recommended:**

Install Stellarscope using [bioconda](https://bioconda.github.io):

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/telescope/README.html) 

```bash
conda install -c bioconda stellarscope
```

See [Getting Started](https://bioconda.github.io/user/install.html) for
instructions on setting up bioconda.


**Alternative:**

Use conda package manager to install dependencies, then 
use `pip` to install Stellarscope.

The following has been testing using miniconda3 on macOS and Linux (CentOS 7):

```bash
conda create -n stellarscope_env python=3.6 future pyyaml cython=0.29.7 \
  numpy=1.16.3 pandas=1.1.3 scipy=1.2.1 pysam=0.15.2 htslib=1.9 intervaltree=3.0.2

conda activate stellarscope_env
pip install git+git://github.com/nixonlab/stellarscope.git
stellarscope assign -h
```

## Testing

A BAM file (`alignment.bam`) and annotation (`annotation.gtf`) are included in
the telescope package for testing. The files are installed in the `data` 
directory of the package root. We've included a subcommand, `stellarscope test`,
to generate an example command line with the correct paths. 
For example, to generate an example command:

```
stellarscope test
```

The command can be executed using `eval`:

```
eval $(stellarscope test)
```

The expected output to STDOUT includes the final log-likelihood, which was 
`95252.596293` in our tests. The test also outputs a report,
`stellarscope-stellarscope_report.tsv`, which can be compared to the report 
included in the `data` directory. NOTE: The precise values may be 
platform-dependent due to differences in floating point precision.

## Usage

### `stellarscope assign`

The `stellarscope assign` program finds overlapping reads between an alignment
(SAM/BAM) and an annotation (GTF) then reassigns reads using a statistical
model. This algorithm enables locus-specific quantification of transposable
element expression.

#### Basic usage

Basic usage requires a file containing read alignments to the genome and an 
annotation file with the transposable element gene model. 
To obtain single-cell TE counts from a BAM/SAM file:

```
stellarscope assign [samfile] [gtffile]
```

The alignment file must be in SAM or BAM format must be collated so that all 
alignments for a read pair appear sequentially in the file. Fragments should be
permitted to map to multiple locations (i.e. `-k` option in `bowtie2`).

The annotation file must be in GTF format and indicate the genomic regions that
represent transposable element transcripts. The transcripts are permitted to be
disjoint in order to exclude insertions of other element types. A collection of
valid transposable element gene models are available for download at 
[mlbendall/telescope_annotation_db](https://github.com/mlbendall/telescope_annotation_db).

#### Advanced usage

```
Input Options:

  samfile               Path to alignment file. Alignment file can be in SAM
                        or BAM format. File must be collated so that all
                        alignments for a read pair appear sequentially in the
                        file.
  gtffile               Path to annotation file (GTF format)
  --attribute ATTRIBUTE
                        GTF attribute that defines a transposable element
                        locus. GTF features that share the same value for
                        --attribute will be considered as part of the same
                        locus. (default: locus)
  --no_feature_key NO_FEATURE_KEY
                        Used internally to represent alignments. Must be
                        different from all other feature names. (default:
                        __no_feature)
  --ncpu NCPU           Number of cores to use. (Multiple cores not supported
                        yet). (default: 1)
  --tempdir TEMPDIR     Path to temporary directory. Temporary files will be
                        stored here. Default uses python tempfile package to
                        create the temporary directory. (default: None)
  --barcodefile BARCODEFILE 
                        Path to file containing cellular barcodes for which TE read 
                        counts should be returned. Optional but recommended to provide
                        a file with barcodes with sufficient numbers of reads. 
                        (default: None)
  --celltypefile CELLTYPEFILE 
                        Path to whitespace-delimited file containing cellular barcodes 
                        and each barcode's associated cell type or cluster. Only used
                        when pooling_mode='celltype'. (default: None)                   

Reporting Options:

  --quiet               Silence (most) output. (default: False)
  --debug               Print debug messages. (default: False)
  --logfile LOGFILE     Log output to this file. (default: None)
  --outdir OUTDIR       Output directory. (default: .)
  --exp_tag EXP_TAG     Experiment tag (default: telescope)
  --updated_sam         Generate an updated alignment file. (default: False)
  
Run Modes:

  --reassign_mode {exclude,choose,average,conf,unique}
                        Reassignment mode. After EM is complete, each fragment
                        is reassigned according to the expected value of its
                        membership weights. The reassignment method is the
                        method for resolving the "best" reassignment for
                        fragments that have multiple possible reassignments.
                        Available modes are: "exclude" - fragments with
                        multiple best assignments are excluded from the final
                        counts; "choose" - the best assignment is randomly
                        chosen from among the set of best assignments;
                        "average" - the fragment is divided evenly among the
                        best assignments; "conf" - only assignments that
                        exceed a certain threshold (see --conf_prob) are
                        accepted; "unique" - only uniquely aligned reads are
                        included. NOTE: Results using all assignment modes are
                        included in the Telescope report by default. This
                        argument determines what mode will be used for the
                        "final counts" column. (default: exclude)
  --conf_prob CONF_PROB
                        Minimum probability for high confidence assignment.
                        (default: 0.9)
  --overlap_mode {threshold,intersection-strict,union}
                        Overlap mode. The method used to determine whether a
                        fragment overlaps feature. (default: threshold)
  --overlap_threshold OVERLAP_THRESHOLD
                        Fraction of fragment that must be contained within a
                        feature to be assigned to that locus. Ignored if
                        --overlap_method is not "threshold". (default: 0.2)
  --annotation_class {intervaltree,htseq}
                        Annotation class to use for finding overlaps. Both
                        htseq and intervaltree appear to yield identical
                        results. Performance differences are TBD. (default:
                        intervaltree)                    
  --stranded_mode {None, RF, R, FR, F}
                        Options for considering feature strand when assigning reads. 
                        If None, for each feature in the annotation, returns counts 
                        for the positive strand and negative strand. If not None, 
                        this argument specifies the orientation of paired end reads 
                        (RF - read 1 reverse strand, read 2 forward strand) and
                        single end reads (F - forward strand) with respect to the 
                        generating transcript. (default: None)
  --pooling_mode {pseudobulk, individual, celltype} 
                        The population of cells from which reads will be used to estimate 
                        locus expression levels in EM. The "pseudobulk" option pools all reads 
                        from all cells, "individual" only uses reads from each individual cell, 
                        and "celltype" uses reads from specified subsets of cells, provided in 
                        a file (--celltypefile). (default: pseudobulk)
  --use_every_reassign_mode 
                        Whether to output count matrices using every reassign mode. 
                        If specified, six output count matrices will be generated, 
                        corresponding to the six possible reassignment methods (all, exclude, 
                        choose, average, conf, unique). (default: False)
  --barcode_tag 
                        The name of the field in the BAM/SAM 
                        file containing the barcode for each read. (default: CB)
  --umi_tag 
                        The name of the field in the BAM/SAM 
                        file containing the unique molecular identifier (UMI) for each read. (default: UB)
Model Parameters:

  --pi_prior PI_PRIOR   Prior on π. Equivalent to adding n unique reads.
                        (default: 0)
  --theta_prior THETA_PRIOR
                        Prior on θ. Equivalent to adding n non-unique reads.
                        (default: 200000)
  --em_epsilon EM_EPSILON
                        EM Algorithm Epsilon cutoff (default: 1e-7)
  --max_iter MAX_ITER   EM Algorithm maximum iterations (default: 100)
  --use_likelihood      Use difference in log-likelihood as convergence
                        criteria. (default: False)
  --skip_em             Exits after loading alignment and saving checkpoint
                        file. (default: False)
```


### `stellarscope resume`

The `stellarscope resume` program loads the checkpoint from a previous run and 
reassigns reads using a statistical model.

#### Basic usage

Basic usage requires a checkpoint file created by an earlier run of 
`stellarscope assign`. Useful if the run fails after the initial load:

```
stellarscope resume [checkpoint]
```

#### Advanced usage

Options are available for tuning the EM optimization, similar to 
`stellarscope assign`.

```
Input Options:

  checkpoint            Path to checkpoint file.

Reporting Options:

  --quiet               Silence (most) output. (default: False)
  --debug               Print debug messages. (default: False)
  --logfile LOGFILE     Log output to this file. (default: None)
  --outdir OUTDIR       Output directory. (default: .)
  --exp_tag EXP_TAG     Experiment tag (default: telescope)
  --updated_sam         Generate an updated alignment file. (default: False)
  
Run Modes:

  --reassign_mode {exclude, choose, average, conf, unique}
                        Reassignment mode. After EM is complete, each fragment
                        is reassigned according to the expected value of its
                        membership weights. The reassignment method is the
                        method for resolving the "best" reassignment for
                        fragments that have multiple possible reassignments.
                        Available modes are: "exclude" - fragments with
                        multiple best assignments are excluded from the final
                        counts; "choose" - the best assignment is randomly
                        chosen from among the set of best assignments;
                        "average" - the fragment is divided evenly among the
                        best assignments; "conf" - only assignments that
                        exceed a certain threshold (see --conf_prob) are
                        accepted; "unique" - only uniquely aligned reads are
                        included. NOTE: Results using all assignment modes are
                        included in the Telescope report by default. This
                        argument determines what mode will be used for the
                        "final counts" column. (default: exclude)
  --conf_prob CONF_PROB
                        Minimum probability for high confidence assignment.
                        (default: 0.9)
  --overlap_mode {threshold, intersection-strict, union}
                        Overlap mode. The method used to determine whether a
                        fragment overlaps feature. (default: threshold)
  --overlap_threshold OVERLAP_THRESHOLD
                        Fraction of fragment that must be contained within a
                        feature to be assigned to that locus. Ignored if
                        --overlap_method is not "threshold". (default: 0.2)
  --annotation_class {intervaltree, htseq}
                        Annotation class to use for finding overlaps. Both
                        htseq and intervaltree appear to yield identical
                        results. Performance differences are TBD. (default:
                        intervaltree)                    
  --stranded_mode {None, RF, R, FR, F}
                        Options for considering feature strand when assigning reads. 
                        If None, for each feature in the annotation, returns counts 
                        for the positive strand and negative strand. If not None, 
                        this argument specifies the orientation of paired end reads 
                        (RF - read 1 reverse strand, read 2 forward strand) and
                        single end reads (F - forward strand) with respect to the 
                        generating transcript. (default: None)
  --pooling_mode {pseudobulk, individual, celltype} 
                        The population of cells from which reads will be used to estimate 
                        locus expression levels in EM. The "pseudobulk" option pools all reads 
                        from all cells, "individual" only uses reads from each individual cell, 
                        and "celltype" uses reads from specified subsets of cells, provided in 
                        a file (--celltypefile). (default: pseudobulk)
  --use_every_reassign_mode
                        Whether to output count matrices using every reassign mode. 
                        If specified, six output count matrices will be generated, 
                        corresponding to the six possible reassignment methods (all, exclude, 
                        choose, average, conf, unique). (default: False)
  --barcode_tag
                        The name of the field in the BAM/SAM 
                        file containing the barcode for each read. (default: CB)
  --umi_tag
                        The name of the field in the BAM/SAM 
                        file containing the unique molecular identifier (UMI) for each read. (default: UB)
Model Parameters:

  --pi_prior PI_PRIOR   Prior on π. Equivalent to adding n unique reads.
                        (default: 0)
  --theta_prior THETA_PRIOR
                        Prior on θ. Equivalent to adding n non-unique reads.
                        (default: 200000)
  --em_epsilon EM_EPSILON
                        EM Algorithm Epsilon cutoff (default: 1e-7)
  --max_iter MAX_ITER   EM Algorithm maximum iterations (default: 100)
  --use_likelihood      Use difference in log-likelihood as convergence
                        criteria. (default: False)
  --skip_em             Exits after loading alignment and saving checkpoint
                        file. (default: False)
```
                        
## Output

Stellarscope has three main output files: the transcript counts estimated via EM (`stellarscope-TE_counts.tsv`), 
a statistical report of the run containing model parameters and additional information
(`stellarscope-stats_report.tsv`), and an updated SAM file (optional). 
The count file is most important for downstream differential
expression analysis. The updated SAM file is useful for downstream locus-specific analyses. 

If the user would like the tool to output count matrices generated
via each of the six assignment methods, they can use the `--use_every_reassign_mode`
option (`stellarscope assign [samfile] [gtffile] --use_every_reassign_mode`).

### Stellarscope statistics report

In addition to outputting transcript counts,
Stellarscope (`stellarscope assign`) provides a more detailed 
statistical report of each read assignment run. 
The first line in the  report is a comment (starting with a “#”) that
contains information about the run such as the number of fragments processed,
number of mapped fragments, number of uniquely and ambiguously mapped 
fragments, and number of fragments mapping to the annotation. The total number
of mapped fragments may be useful for normalization. 

The rest of the report is a table with expression values for 
individual transposable element locations, as well as estimated and initial model parameters.
Comparing the results from different assignment methods may shed light on the 
model's behaviour. The columns of the table are: 

+ `transcript` - Transcript ID, by default from "locus" field. See --attribute argument to use a different attribute.
+ `transcript_length` - Approximate length of transcript. This is calculated from the annotation, not the data, and is equal to the spanning length of the annotation minus any non-model regions.
+ `final_prop` - Final proportion of fragments represented by transcript. This is the final estimate of the π parameter.
+ `init_prop` - Initial proportion of fragments represented by transcript.

### Updated SAM file

The updated SAM file contains those fragments that has at least 1 initial 
alignment to a transposable element. The final assignment and probabilities are
encoded in the SAM tags:

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