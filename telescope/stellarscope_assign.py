# -*- coding: utf-8 -*-
""" Stellarscope assign

"""

from __future__ import print_function
from __future__ import absolute_import
from builtins import super

__author__ = 'Matthew L. Bendall'
__copyright__ = "Copyright (C) 2021 Matthew L. Bendall"

class StellarscopeAssignOptions(utils.SubcommandOptions):
    """

    """
    OPTS = pkgutil.get_data('telescope', 'cmdopts/stellarscope_assign.yaml')

    def __init__(self, args, sc):
        super().__init__(args)

        if hasattr(self, 'tempdir') and self.tempdir is None:
            if hasattr(self, 'ncpu') and self.ncpu > 1:
                self.tempdir = tempfile.mkdtemp()
                atexit.register(shutil.rmtree, self.tempdir)

    def outfile_path(self, suffix):
        basename = '%s-%s' % (self.exp_tag, suffix)
        return os.path.join(self.outdir, basename)

