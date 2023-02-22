# Version History

---

All notable changes to this project will be documented in this file.

## [Unreleased](https://github.com/mlbendall/telescope/tree/main)

### Added

- `CHANGELOG.md` for documenting changes

- `DEVELOPER.md` developer notes

- `environment.yml` conda environment specification

- Included [versioneer](https://github.com/python-versioneer/python-versioneer)
  for managing version strings through git tags

### Changed
- Depends on python >= 3.7, ensure dict objects maintain insertion-order.
  See [What’s New In Python 3.7](https://docs.python.org/3/whatsnew/3.7.html)
- Removed `bulk` and `sc` subcommands since the single-cell 
  [stellarscope](https://github.com/nixonlab/stellarscope) project has moved
  into its own repository. 

### Fixed

- Error where numpy.int is deprecated

----

_The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html)._


## Previous versions

Prior to this CHANGELOG, changes were tracked in [README.md](./README.md).
The previous version history is relocated here:

### v1.0.3.1

  + Checks that GTF feature type is "exon"
  + Skips GTF lines missing attribute (with warning)

### v1.0.3
  + Added cimport statements to calignment.pyx (MacOS bug fix)
  + Fixed warning about deprecated PyYAML yaml.load
  + Compatibility with intervaltree v3.0.2

### v1.0.2
  + Temporary files are written as BAM

### v1.0.1
  + Changed default `theta_prior` to 200,000

### v1.0
  + Removed dependency on `git`
  + Release version

### v0.6.5
  + Support for sorted BAM files
  + Parallel reading of sorted BAM files
  + Improved performance of EM
  + Improved memory-efficiency of spare matrix

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
