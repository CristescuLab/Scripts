import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('ggplot')

def mod_zotu(x):
    return x[1:].capitalize()

def process_basta(df):
    taxonlevels = 'Kingdom Phylum Class Order Family Genus Species'.split()
    ndf = df.tax.str.split(';', expand=True).reindex(columns=range(7))
    ndf = ndf.rename(columns=dict(zip(range(7),taxonlevels)))
    ndf['#OTU ID'] = df.Zotu.apply(mod_zotu)
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
                assert row[-1] == 'Z%s' % row[1].lower()
            except:
                print(row)
                raise
    return m.reindex(columns=cols)

basta90 = pd.read_table('COIzotus_8_90_95.basta', sep='\t', names=[
    'Zotu', 'tax', 'best'], header=None)
basta90 = process_basta(basta90)
basta97 = pd.read_table('COIzotus_8_97_95.basta', sep='\t', names=[
    'Zotu', 'tax', 'best'], header=None)
basta97 = process_basta(basta97)
otutab = pd.read_table('allCOI.otutab.txt', sep='\t')

res90 = mergethem(otutab, basta90)
res90.to_csv('basta90_mergedotutab.tsv', sep='\t', index=False)
res97 = mergethem(otutab, basta97)
res97.to_csv('basta97_mergedotutab.tsv', sep='\t', index=False)

plotcols = 'Kingdom,Phylum,Class,Order,Family,Genus,Species,Zotu,#OTU ID' \
           ''.split(',')
a = res90.reindex(columns=plotcols).count()
b = res97.reindex(columns=plotcols).count()
c = pd.concat([a,b], axis=1).rename(columns={0:90, 1:97})
c.rename(columns={'Zotu': 'Any Basta', '#OTU ID': 'OTU'})
plt.figure()
c.plot.bar()
plt.ylabel('Number of reads with assignment')
plt.tight_layout()
plt.savefig('test.pdf')