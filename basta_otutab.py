from subprocess import run
import pandas as pd
import argparse
import os
# import matplotlib.pyplot as plt
#
# plt.style.use('ggplot')


def mod_zotu(x):
    return x[1:].capitalize()


def process_basta(df):
    taxonlevels = 'Kingdom Phylum Class Order Family Genus Species'.split()
    ndf = df.tax.str.split(';', expand=True).reindex(columns=range(7))
    ndf = ndf.rename(columns=dict(zip(range(7),taxonlevels)))
    ndf['#OTU ID'] = df.Zotu
    ndf['Zotu'] = df.Zotu
    return ndf


def mergethem(otu, basta):
    cols = basta.columns.tolist()
    cols.pop(cols.index('#OTU ID'))
    cols += otu.columns.tolist()
    m = otu.merge(basta, on='#OTU ID', how='outer')
    for row in m.itertuples():
        if pd.isnull(row[-1]):
            continue
        else:
            try:
                assert row[-1] == row[1]
            except:
                print(row)
                raise
    return m.reindex(columns=cols)


def run_basta(blast, outpref, pid, qcov=95, qcol=5, config='config.txt'):
    awk = "<(awk ' $%d > %f ' %s)" % (qcol, qcov, blast)
    args = ['basta', 'sequence', awk, '%s.out' % outpref,  'gb', '-v',
            '%s.verbose' % outpref, '-e',  '0.0001', '-i', str(pid), '-b', 'T',
            '-x', 'T', '-c', config, '-m', '1']
    st = run(' '.join(args), shell=True, executable='/bin/bash')
    return '%s.out' % outpref


def process_one(args):
    outpref = '%s_p%d_q%d' % (args.outprefix, args.pid, args.qcov)
    if not os.path.isfile(args.config):
        with open(args.config, 'w') as c:
            c.write("query_id\t0\nsubject_id\t1\n")
            c.write("align_length\t6\nevalue\t3\npident\t2")
    if args.merge_only is None:
        basta = run_basta(args.blast, outpref, args.pid, qcov=args.qcov,
                          qcol=args.qcol)
    else:
        basta = args.merge_only
    if not args.basta_only:
        basta = pd.read_csv(basta, sep='\t', names=['Zotu', 'tax', 'best'],
                            header=None)
        basta = process_basta(basta)
        otutab = pd.read_csv(args.otutab, sep='\t')
        result = mergethem(otutab, basta)
        result.to_csv('%s_mergedotutab.tsv' % outpref, sep='\t',
                      index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('blast', help='blast file with format 6')
    parser.add_argument('otutab', help='File with the otu table to merge')
    parser.add_argument('-c', '--config', help='File with the mapping of '
                                               'positions and  headers for ' 
                                               'basta', default='config.txt')
    parser.add_argument('-p', '--pid', help='percent identity to pass to basta',
                        type=int)
    parser.add_argument('-b', '--basta_only',
                        help='Run basta without merging with OTU table',
                        default=False, action='store_true')
    parser.add_argument('-o', '--outprefix', help='Prefix for outputs')
    parser.add_argument('-q', '--qcov', help='Minumum query coverage',
                        default=95, type=float)
    parser.add_argument('-r', '--qcol', help='column in blast with query cov',
                        default=5, type=int)
    parser.add_argument('-m', '--merge_only',
                        help='Basta file to be  merged with OTU table. This '
                             'turns off basta execution', default=None)



    args = parser.parse_args()
    process_one(args)



#res97.to_csv('basta97_mergedotutab.tsv', sep='\t', index=False)
# plotcols = 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Zotu,#OTU ID' \
#                ''.split(',')
# a = res90.reindex(columns=plotcols).count()
# b = res97.reindex(columns=plotcols).count()
# c = pd.concat([a,b], axis=1).rename(columns={0:90, 1:97})
# c.rename(columns={'Zotu': 'Any Basta', '#OTU ID': 'OTU'})
# plt.figure()
# c.plot.bar()
# plt.ylabel('Number of reads with assignment')
# plt.tight_layout()
# plt.savefig('test.pdf')