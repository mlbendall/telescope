# -*- coding: utf-8 -*-
""" Setup telescope-ngs package

"""
from __future__ import print_function

import os
from distutils.core import setup
from setuptools import Extension
from setuptools import find_packages

from telescope._version import VERSION

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

USE_CYTHON = False

CONDA_PREFIX = os.environ.get("CONDA_PREFIX", None)
HTSLIB_INCLUDE_DIR = os.environ.get("HTSLIB_INCLUDE_DIR", None)

if HTSLIB_INCLUDE_DIR is None:
    if CONDA_PREFIX is not None:
        HTSLIB_INCLUDE_DIR = os.path.join(CONDA_PREFIX, 'include')

ext = '.pyx' if USE_CYTHON else '.c'
extensions = [
    Extension("telescope.cTelescope",
              ["telescope/cTelescope"+ext],
              include_dirs=[HTSLIB_INCLUDE_DIR, ]
              ),
]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)



setup(
    name='telescope-ngs',
    version=VERSION.split('+')[0],
    packages=find_packages(),

    install_requires=[
        'numpy>=1.13.0',
        'scipy>=0.19.0',
        'intervaltree',
        'pysam>=0.12',
        'pyyaml',
        'future',
        'cython',
    ],

    # Runnable scripts
    entry_points={
        'console_scripts': [
            'telescope=telescope.__main__:main',
        ],
    },

    # cython
    ext_modules=extensions,

    # metadata for upload to PyPI
    author='Matthew L. Bendall',
    author_email='bendall@gwu.edu',
    description='Single locus resolution of Transposable ELEment expression using next-generation sequencing.',
    license='',
    keywords='',
    url='https://github.com/mlbendall/telescope',

    zip_safe=False
)
