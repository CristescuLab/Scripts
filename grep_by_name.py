#!/usr/bin/env python3
from itertools import zip_longest
import gzip
import shelve
import sys
import os
from tqdm import tqdm


def parse_fasta(filename, fn2=None):
    """
    Parse a single fasta file and put it in a shelve DB

    :param filename: fasta file
    :param fn2: name of the shelve
    :param rename: boolean of whether to append prefix to each sequence
    :return: shelve db name
    """

    if isinstance(filename, shelve.DbfilenameShelf):
        Warning('THIS IS ALREADY A SHELF!!!')
        return filename
    elif isinstance(filename, bytes):
        fnc = BytesIO
        fn = 'shelve' if fn2 is None else fn2
    else:
        fnc = open
        fn = filename[:filename.find('.fa')]
    db = '%s.shelve' % fn
    if not os.path.isfile(db):
        dic = shelve.open(db)
        name = None
        seq = ''
        with fnc(filename) as F:
            for line in tqdm(F, desc="Creating shelve"):
                line = line.decode('utf-8') if isinstance(line, bytes) else line
                if line.startswith('>'):
                    if name is not None:
                        dic[name] = seq
                    seq = ''
                    name = line.strip()
                else:
                    seq += line
            if name not in dic:
                dic[name] = seq
    return db


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def fastq_yielder(filename):
    with gzip.open(filename, 'r') as infile:
        for entry in grouper(infile, 4, ''):
            yield entry


def create_shelve(filename, prefix):
    dbname = '%s.shelve' % prefix
    if not os.path.isfile('%s.shelve' % prefix):
        with shelve.open(dbname) as d:
            for sequence in tqdm(fastq_yielder(filename),
                                 desc="Creating shelve"):
                header = sequence[0].strip().split()[0][1:].decode('utf-8')
                d[header] = b''.join(sequence)
    return dbname


def write_fasta(d, name, handle):
    try:
        line = d[name.strip()]
    except KeyError:
        name = '>%s' % name.strip()
        line = d[name]
    handle.write('%s\n%s\n' % (name, line))


def write_fq(d, name, handle):
    handle.write(d[name.strip()])


def main(intersect, filename, prefix, file_type='fq'):
    if file_type == 'fq':
        dbname = create_shelve(filename, prefix)
        outfn = '%s.fastq.gz' % prefix
        fnc = gzip.open
        outfnc = write_fq
    else:
        dbname = parse_fasta(filename, prefix)
        outfn = '%s.fasta' % prefix
        fnc = open
        outfnc = write_fasta

    with open(intersect) as infile, fnc(outfn, 'w') as fq, shelve.open(
            dbname) as d:
        for name in tqdm(infile, desc='Grepping'):
            outfnc(d, name, fq)


if __name__ == '__main__':
    # Usage:
    # python grep_by_name.py intersectfile fqfile prefix <fq|fasta>
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])