"""
**blast_processing.py
** Copyright (C) 2018  Jose Sergio Hleap

Parse a Blast output on format 6 with the following columns:
'qseqid sseqid pident evalue qcovs qlen length staxid stitle'. It assumes that
the stitle contains the lineage or the first two fields are species. If
lineage is in stitle, it will assume 7 taxon levels:
'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

E-mail: jose.hleaplozano@mcgill.ca

Python modules:
1. pandas
2. matplotlib
"""
import csv
import optparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.style.use('ggplot')
import pandas as pd
import numpy as np
import re
from subprocess import check_output
from io import BytesIO

kings = ['Bacteria', ' ;', 'Eukaryota', 'Archaea']
SIX = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
plast_names = 'qseqid sseqid pident length nb_misses nb_gaps qstart qend ' \
              'sstart send e-evalue bit_score qlen query_frame ' \
              'query_translated, qcovs query_gaps subject_length subject_frame ' \
              'subject_translated sovs,subject_gaps'.split(',')
#TODO: check fix for staxid
names = 'qseqid sseqid pident evalue qcovs qlen length stitle'
names = names.split()
# line for taxonkit run
taxonkit = "grep -v '#' %s | cut -f 8 | sort -u| taxonkit lineage| taxonkit reformat"


def get_sps_coi(line):
    """
    Helper when using COI markers from BOLD systems
    :param line: line from blast output
    :return: species string
    """
    return line.split('|')[1]


def get_sps(line):
    """
    Process each blast line
    :param line: One row/line of the blast dataframe
    :return: species string
    """
    if '|' in line:
        return ' '.join(line.split('|')[2].split()[:2])
    elif 'Fish' in line:
       return line.split()[0]
    else:
        # try to assess if there is an accession number before species
        line = line.split()
        if bool(re.search(r'\d', line[0])):
            idx = 1
        else:
            idx = 0
        #return ' '.join(line.split()[1:3])
        return ' '.join(line[idx:3])


def get_lineages(fn):
    """
    If leneages are not in stitle compute them using the field 8 as staxid

    :param fn: blast hits file name
    :return: dataframe with lineages
    """
    o = check_output(taxonkit % fn, shell=True)
    df = pd.read_table(BytesIO(o), header=None, names=['staxid', '_', 'lineage'])
    return df.reindex(columns=['staxid', 'lineage'])


def split_acc_lineage(x):
    # if it does not have a kingdom and there is no space before the ;
    if x.startswith(';'):
        return x
    # check with kigndom
    for i in kings:
        if x.find(i) != -1:
            return x[x.find(i):].strip()


def parse_blast(fn, names, filters={}, top_n_hits=None, output_filtered=False,
                coi=False):
    """
    Parse a blast file, and filter it if required

    :param coi: Use COI fromat to parse species
    :param output_filtered: Output the filtered dataframe
    :param top_n_hits: Number of top blast hits to retain
    :param filters: Doctionary with the filters to include
    :param fn: Filename of the blast to parse
    :return: Dataframe with the information
    """
    fnc = {True: get_sps_coi, False: get_sps}
    df = pd.read_table(fn, sep='\t', header=None, comment='#',
                       quoting=csv.QUOTE_NONE, encoding='utf-8')
    if df.shape[1] > 6:
        by = 'evalue pident qcovs qlen length'.split()
        asc = [True, False, False, False, False]
        if df.shape[1] == 9:
            names = ['qseqid', 'sseqid', 'pident', 'evalue', 'qcovs', 'qlen',
                     'length', 'staxid', 'stitle']
    else:
        names = 'qseqid sseqid pident evalue qcovs stitle'.split()
        by = 'evalue pident qcovs'.split()
        asc = [True, False, False]
    df.rename(columns=dict(zip(range(len(names)), names)), inplace=True)
    if filters:
        query = ' & '.join(
            ['(%s > %d)' % (k, v) if k != 'evalue' else '(%s < %e)' % (k, v)
             for k, v in filters.items()])
        df = df.query(query)
    if top_n_hits is not None:
        args = dict(by=by, ascending=asc)
        df = df.sort_values(**args).groupby('qseqid').head(top_n_hits)
    if df.stitle.str.count(';').mean() > 3:
        # Lineage present, incorporate it
        ndf = df.stitle.apply(lambda x: x[x.find(' ')+1:].strip())
        ndf = ndf.str.split(';', expand=True)
        # Assume 7 level taxonomy
        ndf.rename(columns=dict(zip(range(7), SIX)), inplace=True)
        # Join the dataframes
        df = pd.concat([df, ndf], axis=1)
    elif df.shape[1] == 9 and not pd.isnull(df.staxid).all():
        # if not taxonomic info in stitle but staxid is present, run taxonkit
        lin = get_lineages(fn)
        df = df.merge(lin, on='staxid', how='left')
        df.rename(columns={'stitle':'stitle_old'}, inplace=True)
        df.rename(columns={'lineage': 'stitle'}, inplace=True)
        ndf = df.stitle.apply(split_acc_lineage)
            #lambda x: x.strip().split()[0] if " " in x.strip() else x.strip().split()[0])
        ndf = ndf.str.split(';', expand=True)
        # Assume 7 level taxonomy
        ndf.rename(columns=dict(zip(range(7), SIX)), inplace=True)
        # Join the dataframes
        df = pd.concat([df, ndf], axis=1)
    else:
        # Assume that species is in the first two fields of stitle
        # def get_sps(x): return ' '.join(x.strip().split()[:2])
        df.loc[:, 'species'] = df.stitle.apply(fnc[coi])
    if output_filtered:
        outfn = fn[:fn.rfind('.')]
        df.to_csv('%s_filtered.tsv' %  outfn, sep='\t', index=False, header=False)
    print(df.head())
    print(df.columns)
    return df


def report_any(x, taxlevel):
    """
    For the taxonomic level that is being analyzed, fill it with the lowest rank if
    taxlevel is empty
    :param x: the row (series) where this is being applied
    :param taxlevel: taxonomic level of interest
    :return: new df
    """
    if pd.isnull(x.loc[taxlevel]) or x.loc[taxlevel] == '':
        for i in SIX[SIX.index(taxlevel)+1:]:
            if pd.isnull(x.loc[i]):
                continue
            else:
                x.loc[taxlevel] = x.loc[i]
                return x
        for j in reversed(SIX[:SIX.index(taxlevel)]):
            if pd.isnull(x.loc[j]):
                continue
            else:
                x.loc[taxlevel] = x.loc[j]
                return x
    else:
        return x

    x.loc[taxlevel] = "No taxonomic information"
    return x


def get_reads_per_group(df, prefix, taxlevel='species', min_reads=10):
    """
    Get the number of reads per taxonomic level and the number of unique taxa
    per read

    :param min_reads: Minumum number of reads to retain group
    :param df: filtered blast dataframe
    :param taxlevel: taxonomic
    :return:
    """
    # Get number of reads per taxonomic group
    # Get empty taxlevels
    df = df.apply(report_any, args=(taxlevel,), axis=1)
    cou = df.groupby([taxlevel])['qseqid'].nunique()
    cou.to_csv('%s_number_of_reads_in_%s.tsv' % (prefix, taxlevel), sep='\t')
    # Get number of unique species per read
    re = pd.concat([df.groupby('qseqid')[taxlevel].nunique().rename(
        'No. unique taxa'), df.groupby('qseqid')[taxlevel].unique().rename(
        'Unique taxa')], axis=1).sort_values(by='No. unique taxa', ascending=False)
    re.to_csv('%s_Number_unique_%s_per_read.tsv' % (prefix, taxlevel), sep='\t')
    # List number of unique species above the min_reads
    sps = cou[cou > min_reads].index.unique().to_series()
    sps.to_csv('%s_List_unique_%s.txt' % (prefix, taxlevel), header=False, index=False)


def plot_tax(df, n, taxlevel='species', tax_for_pattern=None, pattern=None,
             suffix=None, min_reads=10):
    """
    Create a diversity bar chart
    :param min_reads: Minimum reads per taxonomic level to retain it.
    :param df: Filtered blast datafrae
    :param n: number of records to retain
    :param taxlevel: Taxonomic level to display
    :param tax_for_pattern: Taxonomic level to restrict by pattern
    :param pattern: pattern to restrict the plotting
    :param suffix: Output suffix (before extension)
    """
    prefix = 'top%dhits_%s' % (n, taxlevel)
    if pattern is not None:
        df = df[df.loc[:, tax_for_pattern].str.contains(pattern)]
        prefix += '_%s' % pattern
    if suffix is not None:
        prefix += '_%s' % suffix
    cols = [taxlevel, 'qseqid']
    toplot = df.reindex(columns=cols)
    #ntaxa = toplot.loc[:, taxlevel].nunique()
    #color = matplotlib.cm.inferno_r(np.linspace(.0, 1., ntaxa))
    toplot = toplot.groupby([taxlevel]).nunique()
    toplot = toplot[toplot.qseqid > min_reads]
    fig, ax = plt.subplots()
    toplot.plot(kind='barh', stacked=False, ax=ax, fontsize=12, legend=False)
                # color=color)
    ax.set_xlabel("Number of reads", fontsize=12)
    ax.set_ylabel(taxlevel, fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)
    plt.tight_layout()
    plt.savefig('%s.pdf' % prefix)
    plt.close()


def main(blast_file, prefix, names, pident=None, eval=None, qcov=None, qlen=None, length=None,
         output_filtered=False, taxlevel='species', min_reads=0, plot=False,
         tax_for_pattern=None, pattern=None, suffix_for_plot=None, ntop=None,
         use_coi=False):
    """
    Execute the code

    :param names: Column names of the blast table
    :param blast_file: File with blast results
    :param pident: Minimum percent identity to retain
    :param eval: Maximum evalue to retain
    :param qcov: Minimum Query coverage to retain
    :param qlen: Minimum query length to retain
    :param length: Minimum alignment length to retain
    :param output_filtered: Output the filtered dataframe to file
    :param taxlevel: Taxonomic level to display
    :param min_reads: Minimum reads per taxonomic level to retain it.
    :param plot: Plot a barchart of number of reads per taxonomic level
    :param tax_for_pattern: Taxonomic level to restrict by pattern
    :param pattern: pattern to restrict the plotting
    :param suffix_for_plot: Output suffix (before extension)
    :param ntop: Number of blast top hits to retain
    :param use_coi: Use COI formatting to parse species
    """
    filters = dict(zip('pident evalue qcovs qlen length'.split(),
                       [pident,   eval,   qcov,   qlen,   length]))
    filters = {k: v for k, v in filters.items() if v is not None}
    kwargs = dict(filters=filters, output_filtered=output_filtered,
                  top_n_hits=ntop, coi=use_coi)
    df = parse_blast(blast_file, names, **kwargs)
    get_reads_per_group(df, prefix, taxlevel=taxlevel, min_reads=min_reads)
    if plot:
        plot_tax(df, ntop, taxlevel=taxlevel, tax_for_pattern=tax_for_pattern,
                 pattern=pattern, suffix=suffix_for_plot, min_reads=min_reads)


if __name__ == '__main__':
    intro='''%prog [options] blast_file prefix
    If in Compute Canada (Graham/Cedar), use:
    module load scipy-stack/2018b'''
    opts = optparse.OptionParser(usage='%s' % intro)
    opts.add_option('--output_filtered', '-o', action='store_true',
                    default=False, help=('Output a TSV with the filtered table'
                                         ' [default: %default]'))

    opts.add_option('--pident', '-p', action='store', type=float, default=None,
                    help='Minimum percent identity [default: %default]')
    opts.add_option('--eval', '-e', action='store', type=float, default=None,
                    help='Maximum evalue [default: %default]')
    opts.add_option('--qcov', '-q', action='store', type=float, default=None,
                    help='Minimum query coverage [default: %default]')
    opts.add_option('--qlen', '-Q', action='store', type=int, default=None,
                    help='Minimum query length [default: %default]')
    opts.add_option('--length', '-l', action='store', type=int, default=None,
                    help='Minimum alignment length [default: %default]')
    opts.add_option('--taxlevel', '-t', action='store', default='species',
                    help='Taxonomic level to display [default: %default]')
    opts.add_option('--min_reads', '-r', action='store', type=int, default=0,
                    help=('Minimum number of reads to retain group '
                          '[default: %default]'))
    opts.add_option('--plot', '-P', action='store_true', default=False,
                    help='Make a barchart with your group [default: %default]')
    opts.add_option('--tax_for_pattern', '-a', default=None,
                    help=('Parental taxonomic level to subset based on pattern'
                          ' [default: %default]'))
    opts.add_option('--pattern', '-b', default=None,
                    help=('Pattern to subset the tax_for_pattern with '
                          '[default: %default]'))
    opts.add_option('--suffix_for_plot', '-s', default=None,
                    help=('Suffix for plot (before extension) '
                          '[default: %default]'))
    opts.add_option('--ntop', '-n', action='store', type=int, default=None,
                    help='Number of hits per query [default: %default]')
    opts.add_option('--use_coi', '-c', action='store_true', default=False,
                    help=('If no special formating in the database and using'
                          ' COI from bold-like DBs [default: %default]'))

    opt, arg = opts.parse_args()
    main(arg[0], arg[1], names, opt.pident, opt.eval, opt.qcov, opt.qlen, opt.length,
         opt.output_filtered, opt.taxlevel, opt.min_reads, opt.plot,
         opt.tax_for_pattern, opt.pattern, opt.suffix_for_plot, opt.ntop,
         opt.use_coi)
