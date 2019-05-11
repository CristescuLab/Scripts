"""
This script will take a set of folders with the results of a blast
and a blast_processing.py results and concatenate them into a single
 table with samples in the column.

 Usage:

 python join_processed_blasts.py <folder_pattern> <outputfile_suffix> <final_ouput_name>

for example if you have a bunch of folders that start with output and within each one
you have two parsing of blasts named test.98_98_number_of_reads_in_species.tsv, AND
test.90_95_number_of_reads_in_species.tsv and if you do:

python join_processed_blasts.py output 98_98 final.tsv

This will create a final.tsv folder, concatenating all 98_98 files within the folder that
contain output in their name

"""
import pandas as pd
from glob import glob
import sys, os

prefix=sys.argv[1]
suffix=sys.argv[2]
files=glob(os.path.join('*%s*' % prefix, '*%s*_number_of_reads_in_species.tsv' % suffix))
df = []
for f in files:
    name = os.path.split(f)[1].split('.')[0]
    df.append(pd.read_table(f, sep='\t', index_col=0, header=None, names=[name]))
df = pd.concat(df,axis=1, sort=False)
df.to_csv(sys.argv[3], sep='\t')
