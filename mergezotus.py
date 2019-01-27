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
    blasts = pd.concat(blasts)
    # filter blast so that qlen length are the same
    blasts = blasts[blasts.qlen == blasts.length]
    # filter out the blasts that are the same qseqid and sseqid
    blasts = blasts[blasts.qseqid != blasts.sseqid].reset_index()
    blasts.to_csv('%s.hits' % out, sep='\t', index=False)
    return blasts


def process_lane_and_otus(string, tables):
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
    df = tab[tab['#OTU ID'].isin([otu])]
    print(df)
    assert not df.qseqid.isin(df.sseqid)
    return df


def rename_fasta(dbname, mapping, outfn):
    """
    Loop over the shelve and write a consolidated deduplicated fasta

    :param dbname: selve database with the merge fastas
    :param mapping: dictionary with name of the sequence and the respective otu
    """
    done=[]
    with shelve.open(dbname) as fas, open(outfn, 'w') as out:
        for header, sequence in fas.items():
            name = header.strip()[1:]
            otu = mapping[name]
            if otu not in done:
                newseq = '>%s\n%s' % (otu, sequence)
                out.write(newseq)
                done.append(otu)


def main(outprefix, fasta_suffix='fasta', zotu_table_suffix='txt', cpus=-1):
    # Get fastas in the current working directory
    files = glob('*.%s' % fasta_suffix)
    table_names = glob('*.%s' % zotu_table_suffix)
    tables = {'%ss' % t[:t.find('tab')]: pd.read_table(t, sep='\t')
              for t in table_names}
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
        blast = parallel_blast(db, fn2, out=outprefix)
    else:
        blast = pd.read_table('%s.hits' % outprefix, sep='\t')
        blast = blast[blast.qseqid != blast.sseqid].reset_index()
    # get all the zotus with name label
    zotus = blast.qseqid.unique().tolist()
    # group the blast by qseqid and sseq id
    grpq = blast.groupby('qseqid')

    # initialize a dataframe with the columns in the zotu tables
    new_zotus = pd.DataFrame(columns=tables[list(tables.keys())[0]].columns)
    count=0
    done = []
    mapping = {}
    # Go over each group, retrieved the hits and the respective zotus table
    for g in zotus:#tqdm(zotus, total=len(zotus), desc="Loping over groups"):
        print('Performing g')
        if g not in done:
            count += 1
            new_otu = 'OTU%d' % count
            done.append(g)
            df = grpq.get_group(g)
            matches = df.sseqid.unique().tolist() + [g]
            # map the group and matches to the new otu
            mapping.update(dict(zip(matches, cycle([new_otu]))))
            done.extend(matches)
            dfs = Parallel(n_jobs=cpus, prefer='threads')(
                delayed(process_lane_and_otus)(match, tables) for match in
                tqdm(matches, desc="Matches in group %s" % g))
            print('TYPE DFs', type(dfs))
            d = dfs[0]
            for i in range(1,len(dfs)):
                d += dfs[i]
            # try:
            #     d = reduce(lambda x, y: x.add(y, fill_value=0), dfs)
            # except TypeError:
            #     print(df)
            #     raise
            # change the name of the OTU
            d.reset_index().loc[count-1, '#OTU ID'] = new_otu
            # append to the dataframe
            new_zotus = new_zotus.append(d)
        else:
            print(g, 'in', done)
    rename_fasta(fn2, mapping, '%s.fas' % outprefix)

if __name__ == '__main__':
    # Usage:
    # python mergezotus.py outprefix fasta_ext zotu_table_ext n_cpus
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]))







