"""
This script takes a list of denoised fasta files and the zotus table per each
and merge the zotus based on matches
"""
import sys
from io import BytesIO
import pandas as pd
from joblib import Parallel, delayed
from subprocess import run, PIPE, CalledProcessError
import shelve
import dill
from glob import glob
from functools import reduce
from itertools import cycle
from tqdm import tqdm
import os
import numpy as np


def parse_fasta(files, fn2):
    """
    Parse a fastas file and put it in a shelve DB and rename it with the prefix
    of the filename. It will also write a combined fasta with the renamed
    entries

    :param files: list of fasta files
    :param fn2: name of the shelve
    :return: shelve db name
    """
    with shelve.open(fn2) as dic, open('%s.fasta' % fn2, 'w') as out:
        for filename in files:
            prefix = filename[:filename.find('.fa')]
            name = None
            seq = ''
            with open(filename) as F:
                for line in tqdm(F, desc="Parsing %s" % filename):
                    if line.startswith('>'):
                        if name is not None:
                            dic[name] = seq
                            out.write('%s\n%s' % (name, seq))
                        seq = ''
                        name = '%s_%s' % (line.strip(), prefix)
                    else:
                        seq += line
                if name not in dic:
                    dic[name] = seq
                    out.write('%s\n%s' % (name, seq))
    return fn2


def stdin_run(args, inpt, **kwargs):
    """
    Run programs in the shell

    :param args: list of arguments to run
    :param inpt: standard input to provide
    :param kwargs: any key word arguments for subprocess.run
    :return: stdout
    """
    inpt = inpt.encode('utf-8') if isinstance(inpt, str) else input
    exe = run(args, input=inpt, stderr=PIPE, stdout=PIPE, **kwargs)
    try:
        exe.check_returncode()
    except CalledProcessError:
        raise Exception(exe.stdout, exe.stderr)
    return exe.stdout


def iter_fasta(fastas, done=[]):
    """
    Iterator over sequences in a fasta database (shelve)

    :param fn: Name of the shleve database
    :param done: list of excusions (useful if relauching)
    """
    for header, sequence in tqdm(fastas.items(), desc="Yielding sequences"):
        if header not in done:
            done.append(header)
            yield '%s\n%s' % (header, sequence)


def parallel_blast(db, query, evalue=1E-50, p_id=100, cpus=-1, out='hit.hits'):
    """
    Run blast in parallel threads

    :param db: Database to query
    :param query: name of the shelve database where sequences are
    :param evalue: evalue for blast run
    :param p_id: Percent identity of blast run
    :param mts: Max target sequences to retrive in the blast run
    :param cpus: Number of cpus to use
    :param out: Outfilename
    :return: dataframe with blast
    """
    outfmt_str = 'qseqid sseqid pident evalue qcovs qlen length'
    args = ['blastn', '-db', db, '-query', '-', '-evalue', str(evalue),
            '-perc_identity', str(p_id), '-outfmt', "6 %s" % outfmt_str]
    with shelve.open(query) as query:
        blasts = Parallel(n_jobs=cpus, prefer='threads')(
            delayed(stdin_run)(args, inp) for inp in iter_fasta(query))
    blasts = [pd.read_table(BytesIO(x), header=None, names=outfmt_str.split())
              for x in blasts]
    blasts = pd.concat(blasts, sort=False)
    # filter blast so that qlen length are the same
    blasts = blasts[blasts.qlen == blasts.length]
    # filter out the blasts that are the same qseqid and sseqid
    blasts.to_csv('%s.hits' % out, sep='\t', index=False)
    return blasts


def process_lane_and_otus(string, tables, new_otu):
    """
    Get the apropriate row for a given otu

    :param string: undescore-delimited
    :param tables: dictionary with the zotus table
    :return: the apropriate row
    """
    bl = string.split('_')
    lane = bl[1]
    otu = bl[0]
    tab = tables[lane]
    df = tab[tab.loc[:, '#OTU ID'] == otu]
    df.loc[:, '#OTU ID'] = new_otu
    if df.empty:
        print(string, 'empty')
        with open('Not_in_zotutab.txt', 'a') as outf:
            outf.write('%s\n' % string)
    return df


def single_execution(outprefix, files, cpus, tables):
    # parse fastas
    if not os.path.isfile('%s.shelve' % outprefix):
        fn2 = parse_fasta(files, '%s.shelve' % outprefix)
    else:
        fn2 = '%s.shelve' % outprefix
    # create a single blast database with the combined fasta
    to_db = '%s.fasta' % fn2
    db = '%s.db' % outprefix
    mkbl = ['makeblastdb', '-in', to_db, '-dbtype', 'nucl', '-parse_seqids',
            '-hash_index', '-out', db]
    if not os.path.isfile('%s.nhr' % db):
        run(mkbl)
    # blast all the sequences against this database
    if not os.path.isfile('%s.hits' % outprefix):
        blast = parallel_blast(db, fn2, out=outprefix, cpus=cpus)
    else:
        blast = pd.read_table('%s.hits' % outprefix, sep='\t')
    # group the blast by qseqid and sseq id
    blast = blast.reset_index(drop=True)
    grps = blast.qseqid.unique()
    df = blast.copy()
    # Go over each group, retrieved the hits and the respective zotus table
    fasta = ''
    zotus = []
    done = []
    count=0
    for g in tqdm(grps, desc="Loping over groups", total=len(grps)):
        if g not in done:
            notu = 'Otu%d' % count
            boolean = (df.qseqid.str.contains(g) | df.sseqid.str.contains(g))
            d = df[boolean]
            df = df[~boolean]
            matches = d.loc[:, ['qseqid', 'sseqid']].unstack().unique()
            dfs = [process_lane_and_otus(match, tables, notu) for match in
                   set(matches)]
            done += matches.tolist()
            zotus.append(reduce(lambda x, y: x.merge(y, on='#OTU ID',
                                                     how='outer'), dfs))
            with shelve.open(fn2) as fas, open('%s.mapping' % outprefix,'w') as o:
                fasta += '>%s\n%s' % (notu, fas['>%s' % g])
                o.write('%s\t%s\n' % (notu, '\t'.join(matches.astype(str))))
            count += 1

    new_zotus = pd.concat(zotus, sort=True, join='outer',ignore_index=True
                          ).reset_index(drop=True).fillna(0)
    return new_zotus, fasta


def main(outprefix, fasta_suffix='fasta', zotu_table_suffix='txt', cpus=-1):
    # Get fastas in the current working directory
    files = glob('*.%s' % fasta_suffix)
    table_names = glob('*.%s' % zotu_table_suffix)
    tables = {'%ss' % t[:t.find('tab')]: pd.read_table(t, sep='\t')
              for t in table_names}
    new_zotus, fasta = single_execution(outprefix, files, cpus, tables)
    new_zotus.to_csv('%s.zotus' % outprefix, sep='\t', index=False)
    with open('%s.fas' % outprefix, 'w') as f:
        f.write(fasta)




if __name__ == '__main__':
    # Usage:
    # python mergezotus.py outprefix fasta_ext zotu_table_ext n_cpus
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))