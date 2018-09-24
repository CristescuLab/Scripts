import pandas as pd
from glob import glob
import sys

prefix=sys.argv[1]
files=glob('%s*/*_number_of_reads_in_species.tsv' % prefix)
df = []
for f in files:
    name = f.split('.')[0]
    df.append(pd.read_table(f, sep='\t', index_col=0, header=None, names=[name]))
df = pd.concat(df,axis=1)
df.to_csv(sys.argv[2], sep='\t')