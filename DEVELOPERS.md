# Developers

Conda environment:

```bash
mamba create -n stellarscope_dev\
 python=3.9.7\
 future= 0.18.2\
 cython=0.29.24\
 pyyaml=5.4.1\
 numpy=1.21.2\
 scipy=1.7.1\
 htslib=1.13\
 pysam=0.17.0\
 HTSeq=0.13.5\
 intervaltree=3.0.2\
 pyranges=0.0.111


conda activate stellarscope_dev
pip install -e .
```

## Extending

Each subcommand is defined in a top-level module named as 
`[basecommand]_[subcommand].py`. The module will contain a subclass of
`utils.SubcommandOptions` that will 