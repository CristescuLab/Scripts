#!/usr/bin/env python3
from itertools import zip_longest
import gzip
import shelve
import sys
import os
from tqdm import tqdm


def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def fastq_yielder(filename):
    with gzip.open(filename, 'r') as infile:
        for entry in grouper(infile, 4, ''):
            yield entry


def create_shelve(filename, prefix):
    if not os.path.isfile('%s.shelve' % prefix):
        with shelve.open('%s.shelve' % prefix) as d:
            for sequence in tqdm(fastq_yielder(filename),
                                 desc="Creating shelve"):
                header = sequence[0].strip().split()[0][1:].decode('utf-8')
                d[header] = b''.join(sequence)


def main(intersect, filename, prefix):
    create_shelve(filename, prefix)
    with open(intersect) as infile, gzip.open('%s.fastq.gz' % prefix,
                                              'w') as fq, shelve.open(
            '%s.shelve' % prefix) as d:
        for name in tqdm(infile, desc='Grepping'):
            fq.write(d[name.strip()])


if __name__ == '__main__':
    # Usage:
    # python grep_by_name.py intersectfile fqfile prefix
    main(sys.argv[1], sys.argv[2], sys.argv[3])