# -*- coding: utf-8 -*-
""" Setup telescope-ngs package

"""
from setuptools import setup, find_packages

from telescope._version import VERSION

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

setup(
    name='telescope-ngs',
    version=VERSION.split('+')[0],
    packages=find_packages(),

    install_requires=[
        'numpy>=1.7.0',
        'scipy>=0.17.0',
        'intervaltree',
        'pysam>=0.8.2.1',
    ],

    # Runnable scripts
    entry_points={
        'console_scripts': [
            'telescope=telescope.__main__:main',
        ],
    },

    # metadata for upload to PyPI
    author='Matthew L. Bendall',
    author_email='bendall@gwu.edu',
    description='Single locus resolution of Transposable ELEment expression using next-generation sequencing.',
    license='',
    keywords='',
    url='https://github.com/mlbendall/telescope',

    zip_safe=False
)
