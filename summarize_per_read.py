import pandas as pd
import sys
from glob import glob


def getsp(x):
    res = []
    for y in x:
        h = y[1:-1].split("'")
        res.extend([z for z in h if z.strip()])
    return res


if len(sys.argv) < 2 or '-h' in sys.argv or '--help' in sys.arv:
    print('Usage:')
    print("python summarize_per_read.py '*_output_12S_Numver_unique_species_"
          "per_read.tsv' 12S")
    print('**NOTE: pay imporant attention to the quotes over the pattern')
    sys.exit()
pattern = sys.argv[1]
prefix = sys.argv[2]
# patter is how to look for the number_unique_species. It should match the ls
# unix command.
gl = glob(pattern)
summ = []
order = ['Sample', 'Count unique reads', 'Sum of unique taxa',
         'Species in sample']
for fn in gl:
    sample = fn.split('_')[0]
    df = pd.read_table(fn, sep='\t')
    ntax = df['No. unique taxa'].sum()
    nreads = df.qseqid.nunique()
    arr = repr(getsp(df['Unique taxa'].unique()))
    summ.append({'Sample': sample, 'Count unique reads': nreads,
                 'Sum of unique taxa': ntax, 'Species in sample': arr})

pd.DataFrame(summ).to_csv('%s_summary.tsv' % prefix, sep='\t', index=False,
                          columns=order)
