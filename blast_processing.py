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
from dask import dataframe as dd
import os
import psutil
import re
from subprocess import run, PIPE
from io import BytesIO
from joblib import Parallel, delayed

kings = ['Bacteria', ' ;', 'Eukaryota', 'Archaea']
SIX = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# line for taxonkit run
taxonkit = "grep -v '#' %s | cut -f 8 | sort -u| taxonkit lineage| taxonkit " \
           "reformat"
taxonkit2 = "taxonkit name2taxid| taxonkit lineage -i 2 | taxonkit reformat " \
            "-i 3"


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
    if '([' in line:
        # is from mitofish
        sp = ' '.join(line.split()[:2])
    elif '|' in line:
        sp = ' '.join(line.split('|')[2].split()[:2])
    elif 'Fish' in line:
       sp = line.split()[0]
    else:
        # try to assess if there is an accession number before species
        line = line.split()
        if bool(re.search(r'\d', line[0])):
            idx = 1
        else:
            idx = 0
        sp = ' '.join(line[idx: idx + 2])
    if 'sp.' in sp:
        return sp.split()[0]
    else:
        return sp


def taxon2exe(sp):
    st = run(taxonkit2, input=sp.encode('utf-8'), shell=True, stdout=PIPE)
    return pd.read_table(BytesIO(st.stdout), names=['species', 'staxid', '_',
                                                    'lineage'])


def get_lineages(fn, typeof=1, cpus=-1):
    """
    If lineages are not in stitle compute them using the field 8 as staxid

    :param fn: blast hits file name
    :return: dataframe with lineages
    """
    if isinstance(typeof, list):
        # Assume you have species and want lineages
        dfs = Parallel(n_jobs=cpus, prefer="threads")(delayed(taxon2exe)(sp)
                                                      for sp in set(typeof))
        df = pd.concat(dfs).reindex(columns=['species', 'staxid', 'lineage'])
    else:
        # assume that staxid is in the 8th column of the hits file
        o = run(taxonkit % fn, shell=True, stdout=PIPE).stdout
        df = pd.read_table(BytesIO(o), header=None,
                           names=['staxid', '_', 'lineage']).reindex(
            columns=['staxid', 'lineage'])

    return df


def split_acc_lineage(x):
    try:
        if pd.isnull(x):
            return
    except:
        print(x, type(x))
        raise
    # if it does not have a kingdom and there is no space before the ;
    try:
        if x.startswith(';'):
            return x
    except AttributeError:
        print(x)
        raise AttributeError
    # check with kigndom
    for i in kings:
        if x.find(i) != -1:
            return x[x.find(i):].strip()

def parse_blast(fn, filters={}, top_n_hits=None, output_filtered=False,
                coi=False, same_blast=None, cpus=-1, taxlevel='species',
                names='qseqid sseqid pident evalue qcovs qlen length staxid '
                      'stitle'):
    """
    Parse a blast file, and filter it if required

    :param coi: Use COI fromat to parse species
    :param output_filtered: Output the filtered dataframe
    :param top_n_hits: Number of top blast hits to retain
    :param filters: Doctionary with the filters to include
    :param fn: Filename of the blast to parse
    :return: Dataframe with the information
    """
    # check how big the file is
    size = os.stat(fn).st_size
    avail = psutil.virtual_memory().available
    kwargs = dict(sep='\t', comment='#', encoding='utf-8',
                  quoting=csv.QUOTE_NONE)
    if size < avail:
        module = pd
    else:
        module = dd
        kwargs['blocksize'] = 5E6
    names = names.split()
    sortable = 'evalue pident qcovs qlen length'.split()
    orientation = [True, False, False, False, False]
    sorts = dict(zip(sortable, orientation))
    fnc = {True: get_sps_coi, False: get_sps}
    if same_blast is not None:
        df = module.read_table(same_blast, **kwargs)
    else:
        df = module.read_table(fn, names=names, **kwargs)

    by = list(set(df.columns).intersection(sortable))
    asc = [sorts[x] for x in by]
    flipped = ['evalue', 'maxqlen']
    if filters:
        query = ' & '.join(['(%s > %d)' % (k, v) if k not in flipped else
                            '(%s < %e)' % (k, v) for k, v in filters.items()])
        if 'maxqlen' in query:
            query = query.replace('maxqlen', 'qlen')
        df = df.query(query)
    if top_n_hits is not None:
        args = dict(by=by, ascending=asc)
        df = df.sort_values(**args).groupby('qseqid').head(top_n_hits)

    if ('staxid' in df.columns) and not df.head(100).staxid.isnull().all()  \
        and (same_blast is None):
        # if not taxonomic info in stitle but staxid is present, run taxonkit
        lin = get_lineages(fn, typeof=df.staxid)
        df = df.merge(lin, on='staxid', how='left')
        df.to_csv('df_before_split.tsv', sep='\t', index=False)
        if 'stitle' in df.columns:
            df.rename(columns={'stitle': 'stitle_old'}, inplace=True)
            df.rename(columns={'lineage': 'stitle'}, inplace=True)
        df.to_csv('df_before_split_acc_lineage.tsv', sep='\t', index=False)
        ndf = df.stitle.apply(split_acc_lineage)
        ndf = ndf.str.split(';', expand=True)
        # Assume 7 level taxonomy
        ndf.rename(columns=dict(zip(range(7), SIX)), inplace=True)
        # Join the dataframes
        df = pd.concat([df, ndf], axis=1)
    elif (df.head(100).stitle.str.count(';').mean() > 3) and (
            same_blast is None):
        # Lineage present, incorporate it
        ndf = df.stitle.apply(lambda x: x[x.find(' ')+1:].strip())
        ndf = ndf.str.split(';', expand=True)
        # Assume 7 level taxonomy
        ndf.rename(columns=dict(zip(range(7), SIX)), inplace=True)
        # Join the dataframes
        df = pd.concat([df, ndf], axis=1)
    elif same_blast is None:
        # Assume that species is in the first two fields of stitle
        df['species'] = df.stitle.apply(fnc[coi])
        if taxlevel != 'species':
            # avoid computation if
            lin = get_lineages('', typeof=df.species.tolist(), cpus=cpus)
            df = df.merge(lin, on='species', how='left')
            ndf = df.stitle.apply(split_acc_lineage)
            ndf = ndf.str.split(';', expand=True)
            df = pd.concat([df, ndf], axis=1)

    if output_filtered:
        df.to_csv('%s_filtered.tsv' % output_filtered, sep='\t', index=False,
                  header=True)
        df.reindex(columns=by + [taxlevel]).groupby(taxlevel).describe(
        ).to_csv('%s_filtered_stats.tsv' % output_filtered, sep='\t')

    return df


def report_any(x, taxlevel, names):
    """
    For the taxonomic level that is being analyzed, fill it with the lowest rank if
    taxlevel is empty
    :param x: the row (series) where this is being applied
    :param taxlevel: taxonomic level of interest
    :return: new df
    """
    names = names.split() if isinstance(names, str) else names
    taxinfo = [None if pd.isnull(i) or (i == '') else i for i in
               x[x.index.difference(names)]]
    if pd.isnull(taxinfo).all():
        x.loc[taxlevel] = "No taxonomic information"
        return x
    elif pd.isnull(x.loc[taxlevel]) or x.loc[taxlevel] == '':
        for i in SIX[SIX.index(taxlevel)+1:]:
            if pd.isnull(x.loc[i]) or x.loc[i] == '':
                continue
            else:
                x.loc[taxlevel] = x.loc[i] + '_from_%s' % i
                return x
        for j in reversed(SIX[:SIX.index(taxlevel)]):
            if pd.isnull(x.loc[j]):
                continue
            else:
                x.loc[taxlevel] = x.loc[j] + '_from_%s' % j
                return x
    else:
        return x


def get_reads_per_group(df, prefix, taxlevel='species', min_reads=10, names=[]
                        ):
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
    df = df.apply(report_any, args=(taxlevel, names,), axis=1)
    cou = df.groupby([taxlevel])['qseqid'].nunique()
    if 'size' in df.columns:
        size = df.groupby([taxlevel])['size'].sum()
        size.name = 'Total'
        cou.name = 'Unique'
        cou = pd.concat((cou, size), axis=1)
    cou.to_csv('%s_number_of_reads_in_%s.tsv' % (prefix, taxlevel), sep='\t')
    # Get number of unique species per read
    re = pd.concat([df.groupby('qseqid')[taxlevel].nunique().rename(
        'No. unique taxa'), df.groupby('qseqid')[taxlevel].unique().rename(
        'Unique taxa')], axis=1).sort_values(by='No. unique taxa',
                                             ascending=False)
    re.to_csv('%s_Number_unique_%s_per_read.tsv' % (prefix, taxlevel),
              sep='\t')
    # List number of unique species above the min_reads
    sps = cou[cou > min_reads].index.unique().to_series()
    sps.to_csv('%s_List_unique_%s.txt' % (prefix, taxlevel), header=False,
               index=False)
    return df


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


def parse_dedup(fn):
    gr = []
    with open(fn) as dedup:
        for line in dedup:
            if line.startswith('>'):
                size = int(line[line.find('size=')+5: line.rfind(';')])
                read = line.strip().split()[0][1:]
                gr.append({'qseqid':read, 'size': size})
    return pd.DataFrame(gr)


def main(blast_file, prefix, names, pident=None, evalue=None, query_len=None,
         query_coverage=None,length=None, output_filtered=False, min_reads=0,
         taxon_level='species', plot=False, tax_for_pattern=None, pattern=None,
         suffix_for_plot=None, n_top=None, use_coi=False, report_dedup=None,
         same_blast=None, cpus=-1, max_qlen=None):
    """
    Execute the code

    :param prefix: Prefix for outputs
    :param names: Column names of the blast table
    :param blast_file: File with blast results
    :param pident: Minimum percent identity to retain
    :param evalue: Maximum evalue to retain
    :param query_coverage: Minimum Query coverage to retain
    :param query_len: Minimum query length to retain
    :param length: Minimum alignment length to retain
    :param output_filtered: Output the filtered dataframe to file
    :param taxon_level: Taxonomic level to display
    :param min_reads: Minimum reads per taxonomic level to retain it.
    :param plot: Plot a barchart of number of reads per taxonomic level
    :param tax_for_pattern: Taxonomic level to restrict by pattern
    :param pattern: pattern to restrict the plotting
    :param suffix_for_plot: Output suffix (before extension)
    :param n_top: Number of blast top hits to retain
    :param use_coi: Use COI formatting to parse species
    """
    filters = dict(zip('pident evalue qcovs qlen length maxqlen'.split(),
                       [pident, evalue, query_coverage, query_len, length,
                        max_qlen]))
    filters = {k: v for k, v in filters.items() if v is not None}
    if output_filtered:
        output_filtered = prefix
    kwargs = dict(filters=filters, output_filtered=output_filtered, cpus=cpus,
                  top_n_hits=n_top, coi=use_coi, same_blast=same_blast,
                  names=names)
    df = parse_blast(blast_file, **kwargs)
    if report_dedup is not None:
        dedup = parse_dedup(report_dedup)
        df = df.merge(dedup, on='qseqid')
    df = get_reads_per_group(df, prefix, taxlevel=taxon_level, names=names,
                             min_reads=min_reads)
    if plot:
        plot_tax(df, n_top, taxlevel=taxon_level, pattern=pattern,
                 tax_for_pattern=tax_for_pattern, suffix=suffix_for_plot,
                 min_reads=min_reads)


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
    opts.add_option('--maxqlen', '-R', action='store', type=int, default=None,
                    help='Maximum query length [default: %default]')
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
    opts.add_option('--report_dedup', '-d', action='store', default=None,
                    help=('Use the dereplication information from this file '
                          '(assumes usearch output)[default: %default]'))
    opts.add_option('--same_blast', '-S', action='store', default=None,
                    help=('Use ta previous filtered file and redo filtering '
                          'with different (more stringent) parameters. Pass '
                          'the filename or None [default: %default]'))
    opts.add_option('--cpus', '-C', action='store', default=-1, type=int,
                    help='number of cpus to use in the run [default: %default]'
                    )
    opts.add_option('--colnames', '-H', action='store',
                    default='qseqid sseqid pident evalue qcovs qlen length '
                            'staxid stitle', help='number of cpus to use in '
                                                  'the run [default: %default]'
                    )

    opt, arg = opts.parse_args()
    main(arg[0], arg[1], opt.colnames, pident=opt.pident, evalue=opt.eval,
         query_coverage=opt.qcov, query_len=opt.qlen, max_qlen=opt.maxqlen,
         length=opt.length, output_filtered=opt.output_filtered,
         taxon_level=opt.taxlevel, min_reads=opt.min_reads, plot=opt.plot,
         pattern=opt.pattern, tax_for_pattern=opt.tax_for_pattern,
         n_top=opt.ntop, cpus=opt.cpus, use_coi=opt.use_coi,
         suffix_for_plot=opt.suffix_for_plot, report_dedup=opt.report_dedup,
         same_blast=opt.same_blast)
