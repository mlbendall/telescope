# -*- coding: utf-8 -*-
""" Set the version information for the package

"""

import os
import subprocess

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2016 Matthew L. Bendall"

VERSION = '0.3.1'

# Add git hash to version number (if possible)
wd = os.getcwd()
try:
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'],
                                       stderr = subprocess.STDOUT)
    VERSION = '%s+dev%s' % (VERSION, git_hash.strip())
except subprocess.CalledProcessError as e:
    # The intended behavior for installed versions is to return the version
    # without a local identifier (i.e. devdev89c2fe5). Since installed
    # versions are not git repositories, the call to "git rev-parse" will fail.
    # Simply ignore the error to get the desired VERSION string.
    # print >> sys.stderr, 'Returncode: %s\n%s' % (e.returncode, e.message)
    pass
finally:
    os.chdir(wd)
