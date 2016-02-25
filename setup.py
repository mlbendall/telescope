from setuptools import setup, find_packages

import re
VERSIONFILE="telescope/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name = 'telescope-ngs',
    version = verstr,
    packages = find_packages(),
    scripts = ['bin/telescope'],
    
    # Project uses reStructuredText, so ensure that the docutils get
    # installed or upgraded on the target machine
    # install_requires = ['docutils>=0.3'],

    #package_data = {
    #    # If any package contains *.txt or *.rst files, include them:
    #    '': ['*.txt', '*.rst'],
    #    # And include any *.msg files found in the 'hello' package, too:
    #    'hello': ['*.msg'],
    #},

    # metadata for upload to PyPI
    author='Matthew L. Bendall',
    author_email='bendall@gwu.edu',
    description='Single locus resolution of Transposable ELEment expression using next-generation sequencing.',
    license='',
    keywords='',
    url='https://github.com/mlbendall/telescope',

    zip_safe=False
)
