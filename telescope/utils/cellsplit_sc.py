#! /usr/bin/env python

import sys
import os
import pysam
import argparse
from collections import defaultdict, Counter

############################################
# Check if running as snakemake script
try:
    snakemake
    snakemake_flag = True
except:
    snakemake_flag = False
############################################
def split_bam_by_cell_barcode(bamfile, selected_barcodes_file, dest, log, barcode_tag, umicode_tag):

    import pandas as pd

    bam_in = pysam.AlignmentFile(bamfile) #creates AlignmentFile object
    bam_header = str(bam_in.header).strip() #get the header for the large bamfile
    file_handles = {} #dictionaries of filehandles
    reads_per_umis = defaultdict(set) #umis counter per cell
    reads_per_barcode = Counter() #count how many barcodes
    # data_name,extension = os.path.basename(bamfile).split('.')
    data_name = '.'.join(os.path.basename(bamfile).split('.')[0:-1])
    extension = os.path.basename(bamfile).split('.')[-1]
    selected_barcodes = set(pd.read_csv(selected_barcodes_file, header=None)[0].tolist()) #read the selected barcodes and turn them into a list

    if not os.path.exists(dest): #Make directories if they dont exists
        os.makedirs(dest) #if not, create corresponding directories


    for read in bam_in.fetch(until_eof=True):#For each read in bamfile
        if(read.has_tag(barcode_tag) and read.has_tag(umicode_tag)): # if the read has the selected barcode
            cbc = read.get_tag(barcode_tag) #get the barcode
            umi_code = read.get_tag(umicode_tag) #get the umicode

            if(cbc.split('-')[0] in selected_barcodes): #if the read is in the filtered barcodes file
                reads_per_barcode[cbc] += 1 #counter for reads per cellbarcode
                reads_per_umis[cbc].add(umi_code) #store the umi for the cellbarcode

                if(cbc not in file_handles):# if its not already created, create file handle
                    file_handle = os.path.join(dest, 'cbc_%s.%s.sam' % (cbc, data_name)) #create the file handle
                    file_handles[cbc] = file_handle #add the filehandle to dictionary of filehandles
                    #open the file and write the header
                    with open(file_handle, 'w') as f: #create and open the file
                        print(bam_header, file=f) #append the header
                        print(read.to_string(), file=f) #append the read
                        f.close() #close the file
                else:
                    with open(file_handles[cbc], 'a') as f: #if it already exists
                        print(read.to_string(), file=f) #print the aignment
                        f.close() #close the file


    ##########
    #Create File of Filenames
    # fofn_path = '/'.join(os.path.abspath(log).split('/')[0:-1])+'/sam_filenames.fofn'
    fofn_path = os.path.join(dest,'sam_filenames.fofn')
    # print(fofn_path)
    fofn = open(fofn_path, 'w')
    # exit()

    # Print report to the log to either file/sterr
    if log is None:
        outh = sys.stderr #print log to stderr
    else:
        outh = open(log, 'w') #if log exists, open it
    # Print log header
    print('barcode\treads_per_barcode\tumis_per_cell\tfile_path', file=outh)
    for barcode in file_handles.keys(): # log(3) : cell_barcode    number_of_reads    path_to_bam
        # print(barcode)
        n_umis_per_cell = len(set(reads_per_umis[barcode]))
        print(barcode+'\t'+str(reads_per_barcode[barcode])+'\t'+str(n_umis_per_cell)+'\t'+str(file_handles[barcode]),file=outh)
        sam_filename = str(file_handles[barcode]).split('/')[-1]
        print(sam_filename,file=fofn) #print the filename
    #     # print('%s\t%d\t%d' % (barcode, float(n_reads_per_cell[barcode]), paths[barcode]), file=outh)
    # #Close the log
    fofn.close()
    if log is not None:
        outh.close()
############################################
# If the snakemake flag is true / run with snakemake
if snakemake_flag:

    if snakemake.params and 'barcode_tag' in snakemake.params:
        barcode_tag = snakemake.params.barcode_tag
    if snakemake.params and 'umicode_tag' in snakemake.params:
        umicode_tag = snakemake.params.umicode_tag
    #bamfile, selected_barcodes_file, dest, log, barcode_tag, umicode_tag
    split_bam_by_cell_barcode(snakemake.input.bam, snakemake.input.tsv, snakemake.output[0], snakemake.output[1],barcode_tag,umicode_tag)



# If not runned by snakemake

#If not snakemake get parameters from command line
elif __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='''Split a BAM file by single cell barcode. It is assmued that the cell
                       has barcode in the 'BC' tag '''
    )

    parser.add_argument("--bam",
        help="Path to BAM"
    )

    parser.add_argument("--brcds",
        help="Path to filtered barcodes from starsolo"
    )
    parser.add_argument("--dest",
        help='''Destination path. A BAM file for each cell barcode is created within the
                destination path (i.e. cbcXXXXXX). Default is in the same directory as
                the source BAM.'''
    )

    parser.add_argument("--log",
        help='Filename to write barcode count information. Default is STDERR'
    )

    parser.add_argument("--barcode_tag", default='CB',
        help='Delimiter to split cel_barcode . Default is "CB"'
    )

    parser.add_argument("--umicode_tag", default='UB',
        help='Delimiter to get the umi_code. Default is "UB"'
    )

    args = parser.parse_args()
    split_bam_by_cell_barcode(args.bam, args.brcds,args.dest, args.log,args.barcode_tag,args.umicode_tag)
