from __future__ import print_function
from __future__ import absolute_import

from . import utils
from .utils.cellmerge import CellMerge
import logging as lg

class CellMergeOptions(utils.SubcommandOptions):

    OPTS = """
    - Input Options:
        - file_of_filenames:
            type: str
            default: sam_filenames.fofn
            help: Path to file of filenames "sam_filenames.fofn" outputted by telescope cellsplit
        - cluster_file:
            default: None
            help: Path to a comma-delimited file mapping cell barcodes (first column) to clusters (second column). Optional.
    - Reporting Options:
        - quiet:
            action: store_true
            help: Silence (most) output.
        - debug:
            action: store_true
            help: Print debug messages.
        - logfile:
            type: argparse.FileType('r')
            help: Log output to this file.
        - outdir:
            default: .
            help: Output directory.
        - outfile_name:
            default: merged.bam
            help: Name of outputted merged bam files.
    """

def run(args):
    """

    Args:
        args:

    Returns:

    """
    opts = CellMergeOptions(args)

    utils.configure_logging(opts)
    lg.info('\n{}\n'.format(opts))

    lg.info('Merging alignment files...')

    merger = CellMerge(opts)
    merger()
