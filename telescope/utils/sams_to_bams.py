#! /usr/bin/env python

import os
import argparse
############################################
# Check if running as snakemake script
try:
    snakemake
    snakemake_flag = True
except:
    snakemake_flag = False
############################################
def sams2bams(sam_folder):

    sam_files = [_ for _ in os.listdir(sam_folder) if _.endswith('.sam')] #get all filenames

    remove_sam_fofn = 'rm '+sam_folder+'*.fofn'
    os.system(remove_sam_fofn) #remove sam file of filenames
    bam_fofn = sam_folder+'bam_filenames.fofn' #path to fofn
    bam_fofn_fh = open(bam_fofn, "w") #open the file of filenames

    for sam_file in sam_files:

        data_name = '.'.join(sam_file.split('.')[0:-1]) #get the dataname w/o extension
        bam_file = data_name+'.bam' #add the bam extension
        sam_path = sam_folder+sam_file #generate the path
        bam_path = sam_folder+bam_file #generate the path

        cmd = 'cat '+ sam_path+' | samtools view -u -S | samtools sort -n > '+bam_path #conversion command
        remove_cmd = 'rm '+sam_path #removal command

        os.system(cmd) #run the conversion
        os.system(remove_cmd) #remove sam file
        print(bam_file, file = bam_fofn_fh)

    bam_fofn_fh.close() #close the fofn file
############################################
if snakemake_flag:
    sams2bams(snakemake.input)

elif __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''Converts all sam files of a folder to bam files (*deletes sam files)''')
    parser.add_argument("--sam_folder","-sf",help="Path to folder w/sam files")
    args = parser.parse_args()
    sams2bams(args.sam_folder)
