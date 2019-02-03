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
    df = tab[tab.loc[:,'#OTU ID'] == otu]
    df.loc[:, '#OTU ID'] = new_otu
    if df.empty:
        print(string, 'empty')
        with open('Not_in_zotutab.txt', 'a') as outf:
            outf.write('%s\n' % string)
    return df


def rename_fasta(dbname, mapping, outfn):
    """
    Loop over the shelve and write a consolidated deduplicated fasta

    :param dbname: selve database with the merge fastas
    :param mapping: dictionary with name of the sequence and the respective otu
    """
    if not os.path.isfile(outfn):
        done = []
        count = max([int(x[3:]) for x in set(mapping.values())])
        with shelve.open(dbname) as fas, open(outfn, 'w') as out:
            for header, sequence in fas.items():
                name = header.strip()[1:]
                try:
                    otu = mapping[name]
                except KeyError:
                    # is a singleton
                    count += 1
                    otu = 'OTU%d' % count
                    print(name, 'is a singleton not in mapping. Assigning it '
                                'to new OTU', otu)
                if otu not in done:
                    newseq = '>%s\n%s' % (otu, sequence)
                    out.write(newseq)
                    done.append(otu)


def loop_zotus(df, count, tables):
    new_otu = 'OTU%d' % count
    matches = df.sseqid.unique().tolist() + df.qseqid.unique().tolist()
    # map the group and matches to the new otu
    mapping = dict(zip(set(matches), cycle([new_otu])))
    dfs = [process_lane_and_otus(match, tables, new_otu) for match in
           set(matches)]
    # merge zotus
    d = reduce(lambda x, y: x.merge(y, on='#OTU ID', how='outer'), dfs)
    return d, mapping


def loop_zotus2(df, count, tables):
    new_otu = 'OTU%d' % count
    matches = df.sseqid.unique().tolist()
    # map the group and matches to the new otu
    mapping = dict(zip(set(matches), cycle([new_otu])))
    dfs = [process_lane_and_otus(match, tables, new_otu) for match in
           set(matches)]
    # merge zotus
    d = reduce(lambda x, y: x.merge(y, on='#OTU ID', how='outer'), dfs)
    return d, mapping


def single_execution(outprefix, files, cpus, tables, second=False):
    # parse fastas
    fnc = loop_zotus2 if second else loop_zotus
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
        blast = parallel_blast(db, fn2, out=outprefix)
    else:
        blast = pd.read_table('%s.hits' % outprefix, sep='\t')
    # group the blast by qseqid and sseq id
    blast = blast.reset_index(drop=True)
    grpq = blast.groupby('qseqid')
    print(grpq.size())
    # Go over each group, retrieved the hits and the respective zotus table
    if os.path.isfile('par.dump'):
        with open('par.dump', 'rb') as dump:
            par = dill.load(dump)
    else:
        par = Parallel(n_jobs=cpus, prefer='threads')(
            delayed(fnc)(df[1], count, tables) for count, df in
            tqdm(enumerate(grpq), desc="Loping over groups"))
    with open('par.dump', 'wb') as p:
        dill.dump(par, p)
    new_zotus, mapping = zip(*par)
    mapping = {k: v for d in mapping for k, v in d.items()}
    new_zotus = pd.concat(new_zotus, sort=True, join='outer',
                          ignore_index=True).reset_index(drop=True).fillna(0)
    return new_zotus, mapping, fn2


def main(outprefix, fasta_suffix='fasta', zotu_table_suffix='txt', cpus=-1):
    # Get fastas in the current working directory
    files = glob('*.%s' % fasta_suffix)
    table_names = glob('*.%s' % zotu_table_suffix)
    tables = {'%ss' % t[:t.find('tab')]: pd.read_table(t, sep='\t')
              for t in table_names}
    # execute first pass
    new_zotus, mapping, fn2 = single_execution(outprefix, files, cpus, tables)
    outfas = '%s.fas' % outprefix
    outprefix2 = '%s2' % outprefix
    tables.update({outprefix2: new_zotus})
    rename_fasta(fn2, mapping, outfas)
    # execute second pass
    new_zotus, mapping, fn3 = single_execution(outprefix2, [outfas], cpus,
                                               tables)
    rename_fasta(fn3, mapping, outfas)
    # Second blast to make sure not duplicates in result

    # fn3 = parse_fasta([outfas], '%s.shelve' % outprefix2)
    # db = '%s.db' % outprefix2
    # mkbl = ['makeblastdb', '-in', outfas, '-dbtype', 'nucl', '-parse_seqids',
    #         '-hash_index', '-out', db]
    # run(mkbl)
    # blast2 = parallel_blast(db, fn3, out=outprefix2)
    # blast2 = blast2.reset_index(drop=True)
    # grpq2 = blast2.groupby('qseqid')
    # if os.path.isfile('par.dump'):
    #     with open('par.dump', 'rb') as dump:
    #         par = dill.load(dump)
    # else:
    #     par2 = Parallel(n_jobs=cpus, prefer='threads')(
    #     delayed(loop_zotus)(df[1], count, tables) for count, df in
    #     tqdm(enumerate(grpq2), desc="Loping over groups, second pass"))
    # with open('par.dump', 'wb') as p:
    #     dill.dump(par, p)
    #new_zotus, mapping = zip(*par2)
    #mapping = {k: v for d in mapping for k, v in d.items()}
    new_zotus = pd.concat(new_zotus, sort=True, join='outer',
                          ignore_index=True).reset_index(drop=True).fillna(0)

    new_zotus.to_csv('%s.zotus' % outprefix, sep='\t', index=False)


if __name__ == '__main__':
    # Usage:
    # python mergezotus.py outprefix fasta_ext zotu_table_ext n_cpus
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))







