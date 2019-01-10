import pandas as pd
import sys
from glob import glob

def getsp(x):
    res = []
    for y in x:
        h = y[1:-1].split("'")
        res.extend([z for z in h if z.strip()])
    return res

pattern = sys.argv[1]
# patter is how to look for the number_unique_species. It should match the ls
# unix command.
gl = glob(pattern)
summ = []
for fn in gl:
    sample = fn.split('_')[0]
    df = pd.read_table(fn, sep='\t')
    ntax = df['No. unique taxa'].sum()
    nreads = df.qseqid.nunique()
    arr = repr(getsp(df['Unique taxa'].unique()))
    summ.append({'Sample': sample, 'Count unique reads': nreads,
                 'Sum of unique taxa': ntax, 'Species in sample': arr})

pd.DataFrame(summ).to_csv('summary.tsv', sep='\t', index=False)