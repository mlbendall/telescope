import sys
import os
import pysam
import pandas as pd
import logging as lg
from collections import defaultdict, Counter

class CellSplit:

    def __init__(self, opts):
        self.opts = opts

    def __call__(self):

        bam_in = pysam.AlignmentFile(self.opts.samfile)  # creates AlignmentFile object
        bam_header = str(bam_in.header).strip()  # get the header for the large bamfile
        file_handles = {}  # dictionaries of filehandles
        reads_per_umis = defaultdict(set)  # umis counter per cell
        reads_per_barcode = Counter()  # count how many barcodes
        data_name = '.'.join(os.path.basename(self.opts.samfile).split('.')[0:-1])
        extension = os.path.basename(self.opts.samfile).split('.')[-1]
        selected_barcodes = set(pd.read_csv(self.opts.barcodefile, header=None)[0].tolist())  # read the selected barcodes and turn them into a list
        barcode_alignments = defaultdict(list)

        if not os.path.exists(self.opts.outdir):  # Make directories if they dont exists
            os.makedirs(self.opts.outdir)  # if not, create corresponding directories

        for read in bam_in.fetch(until_eof=True):  # For each read in bamfile

            if read.has_tag(self.opts.barcode_tag) and read.has_tag(self.opts.umi_tag):  # if the read has the selected barcode
                cbc = read.get_tag(self.opts.barcode_tag)  # get the barcode
                umi_code = read.get_tag(self.opts.umi_tag)  # get the umicode

                if cbc.split('-')[0] in selected_barcodes:  # if the read is in the filtered barcodes file
                    reads_per_barcode[cbc] += 1  # counter for reads per cellbarcode
                    reads_per_umis[cbc].add(umi_code)  # store the umi for the cellbarcode

                    if cbc not in file_handles:  # if its not already created, create file handle
                        file_handle = os.path.join(self.opts.outdir, 'cbc_%s.%s.sam' % (cbc, data_name))  # create the file handle
                        file_handles[cbc] = file_handle  # add the filehandle to dictionary of filehandles

                    barcode_alignments[cbc].append(read.to_string())

                    if sum(reads_per_barcode.values()) % 2.5e6 == 0 and not self.opts.quiet:
                        lg.info(f'...{sum(reads_per_barcode.values()) / 1e6:.1f}M alignments processed')

        # write cell-level bam files
        for cbc, file in file_handles.items():
            with open(file, 'w') as f:
                f.write(bam_header)
                for read in barcode_alignments[cbc]:
                    f.write(read)

        ##########
        # Create File of Filenames
        fofn_path = os.path.join(self.opts.outdir, 'sam_filenames.fofn')
        fofn = open(fofn_path, 'w')

        # Print report to the log to either file/sterr
        if self.opts.logfile is None:
            outh = sys.stderr  # print log to stderr
        else:
            outh = open(self.opts.logfile, 'w')  # if log exists, open it

        # Print log header
        print('barcode\treads_per_barcode\tumis_per_cell\tfile_path', file=outh)
        for barcode in file_handles.keys():  # log(3) : cell_barcode    number_of_reads    path_to_bam
            n_umis_per_cell = len(set(reads_per_umis[barcode]))
            print(barcode + '\t' + str(reads_per_barcode[barcode]) + '\t' +
                  str(n_umis_per_cell) + '\t' + str(file_handles[barcode]),
                  file=outh)
            sam_filename = str(file_handles[barcode]).split('/')[-1]
            print(sam_filename, file=fofn)  # print the filename

        # #Close the log
        fofn.close()
        if self.opts.logfile is not None:
            outh.close()