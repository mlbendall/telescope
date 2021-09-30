from __future__ import print_function
from __future__ import absolute_import
from builtins import super

from . import utils
import logging as lg

class CellSplitOptions(utils.SubcommandOptions):

    OPTS = """
    - Input Options:
        - samfile:
            positional: True
            help: Path to alignment file. Alignment file can be in SAM or BAM
                  format. File must be collated so that all alignments for a
                  read pair appear sequentially in the file.
        - barcodefile:
            positional: True
            help: Path to barcode file, containing barcodes from which reads should be kept.
        - barcode_tag:
            type: str
            default: CB
            help: Name of the field in the BAM/SAM file containing the barcode for each read.
        - umi_tag:
            type: str
            default: UB
            help: Name of the field in the BAM/SAM file containing the unique molecular identifier (UMI) for each read.
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
    """

def run(args):
    """

    Args:
        args:

    Returns:

    """
    opts = CellSplitOptions(args)

    utils.configure_logging(opts)
    lg.info('\n{}\n'.format(opts))

    lg.info('Splitting alignment file by cell...')

    splitter = utils.cellsplit.CellSplit(opts)
    splitter()
