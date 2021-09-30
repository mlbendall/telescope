import pandas as pd
from subprocess import call

class CellMerge:

    def __init__(self, opts):
        self.opts = opts

    def __call__(self):

        # read in file of file names, map barcodes to filenames
        with open(self.opts.file_of_filenames, 'r') as f:
            barcode_filenames = {}
            for line in f.readlines():
                barcode, filename = line.split('\t')
                barcode_filenames[barcode] = filename

        # if clusters file is specified, loop through clusters and concatenate bam files within cluster
        if self.opts.cluster_file is not None:
            clusters = pd.read_csv(self.opts.cluster_file)
            for cluster, df in clusters.groupby(by=clusters.columns.iloc[1]):
                cluster_barcodes = df.iloc[:,0]
                cluster_files = ' '.join([barcode_filenames.get(bc) for bc in
                                          cluster_barcodes if barcode_filenames.get(bc)]
                                         )
                out = self.opts.outfile_name
                outfile_name = out[:out.rfind('.')] + f'cluster{cluster}' + '.bam'
                call(f'samtools cat -o {outfile_name} {cluster_files}', shell=True)
        # otherwise just merge all bam files
        else:
            files = ' '.join(barcode_filenames.values())
            call(f'samtools cat -o {self.opts.outfile_name} {files}', shell=True)
