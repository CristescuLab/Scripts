"""
Count taxonomic assignemntes on mbasta
Usage:
python count_basta.py <basta output> <outfile>
"""
import sys
import pandas as pd
from collections import Counter

df = pd.read_table(sys.argv[1], sep='\t', header=None, names=['OTU', 'LCA',
                                                              'best_hit'])
tax = df.LCA.str.split(';', expand=True).reindex(columns=range(7))
with open(sys.argv[2], 'w') as out:
    for i in range(7):
        out.write('RANK %d\n--------\n' % i)
        col = tax.iloc[:, i]
        c = Counter(col)
        line=''
        for t, count in c.most_common():
            if t == '':
                t='Unassigned'
            line += '%s: %d\n' % (t, count)
        out.write('%s\n\n' % line)


