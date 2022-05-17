import pandas as pd
import scipy
from scipy import io
from collections import defaultdict, OrderedDict

from .model import Telescope


class Stellarscope(Telescope):

    def __init__(self, opts):

        super().__init__(opts)
        self.single_cell = True
        self.mapped_read_barcodes = {}  # Dictionary for storing alignment ids mapped to barcodes
        self.mapped_read_umis = {}  # Dictionary for storing alignment ids mapped to UMIs
        self.barcode_read_indices = defaultdict(list)  # Dictionary for storing barcodes mapped to read indices
        self.barcode_umis = defaultdict(list)  # Dictionary for storing barcodes mapped to UMIs

    def output_report(self, tl, stats_filename, counts_filename, barcodes_filename, features_filename):

        _rmethod, _rprob = self.opts.reassign_mode, self.opts.conf_prob
        _fnames = sorted(self.feat_index, key=self.feat_index.get)
        _flens = self.feature_length

        ''' Only output stats file for pseudobulk '''
        if self.opts.pooling_mode == 'pseudobulk':
            _stats_rounding = pd.Series(
                [2, 3, 2, 3],
                index=['final_conf',
                       'final_prop',
                       'init_best_avg',
                       'init_prop']
            )

            # Report information for run statistics
            _stats_report0 = {
                'transcript': _fnames,  # transcript
                'transcript_length': [_flens[f] for f in _fnames],  # tx_len
                'final_prop': tl.pi,  # final_prop
                'init_prop': tl.pi_init  # init_prop
            }

            # Convert report into data frame
            _stats_report = pd.DataFrame(_stats_report0)

            # Sort the report by transcript proportion
            _stats_report.sort_values('final_prop', ascending=False, inplace=True)

            # Round decimal values
            _stats_report = _stats_report.round(_stats_rounding)

            # Run info line
            _comment = ["## RunInfo", ]
            _comment += ['{}:{}'.format(*tup) for tup in self.run_info.items()]

            with open(stats_filename, 'w') as outh:
                outh.write('\t'.join(_comment) + '\n')
                _stats_report.to_csv(outh, sep='\t', index=False)

        ''' Aggregate fragment assignments by cell using each of the 6 assignment methods'''
        _methods = ['conf', 'all', 'unique', 'exclude', 'choose', 'average']
        _allbc = self.barcodes
        _bcidx = OrderedDict(
            {bcode: rows for bcode, rows in self.barcode_read_indices.items() if len(rows) > 0}
        )
        _bcumi = OrderedDict(
            {bcode: umis for bcode, umis in self.barcode_umis.items() if len(_bcidx[bcode]) > 0}
        )

        ''' Write cell barcodes and feature names to a text file '''
        pd.Series(_allbc).to_csv(barcodes_filename, sep='\t', index=False, header=False)
        pd.Series(_fnames).to_csv(features_filename, sep='\t', index=False, header=False)

        for _method in _methods:

            if _method != _rmethod and not self.opts.use_every_reassign_mode:
                continue

            if self.opts.use_every_reassign_mode:
                counts_outfile = counts_filename[:counts_filename.rfind('.')] + '_' + _method + '.mtx'
            else:
                counts_outfile = counts_filename

            _assignments = tl.reassign(_method, _rprob)
            _cell_count_matrix = scipy.sparse.dok_matrix((len(_allbc), _assignments.shape[1]))

            for i, _bcode in enumerate(_allbc):
                ''' If the barcode has reads that map to the annotation, sum the barcode's reads '''
                if _bcode in _bcidx:
                    _rows = _bcidx[_bcode]
                    _umis = _bcumi[_bcode]
                    _cell_assignments = _assignments[_rows, :]
                    _cell_assignment_matrix = scipy.sparse.lil_matrix(_cell_assignments)
                    _cell_final_assignments = _cell_assignments.argmax(axis=1)
                    _umi_assignments = pd.Series(
                        [(umi, assignment) for umi, assignment in zip(_umis, _cell_final_assignments.A1)]
                    )
                    _duplicate_umi_mask = _umi_assignments.duplicated(keep='first').values
                    _cell_assignment_matrix[_duplicate_umi_mask, :] = 0
                    _cell_count_matrix[i, :] = _cell_assignment_matrix.sum(0).A1
                else:
                    _cell_count_matrix[i, :] = 0

            io.mmwrite(counts_outfile, _cell_count_matrix)  # include nofeat
