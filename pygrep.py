from itertools import zip_longest
from joblib import Parallel, delayed
from typing import Tuple
from tqdm import tqdm
import shelve
import gzip
import os

class FastQ:
    dictionary = None
    tipo = None
    prefix = None

    def __init__(self, filename: str, cpus: int = -1) -> None:
        self.filename = filename
        self.cpus = cpus
        self.parser()

    @property
    def filename(self):
        return self.__filename

    @filename.setter
    def filename(self, filename):
        if 'gz' in filename:
            self.open = gzip.open
            self.tipo = 'gz'
            self.prefix = filename[:filename.find('.fastq.gz')]
        else:
            self.open = open
            self.prefix = filename[:filename.rfind('.')]
        self.__filename = filename

    def yielder(self):
        with self.open(self.filename, 'r') as infile:
            args = [iter(infile)] * 4
            for entry in zip_longest(*args, fillvalue=''):
                if self.tipo is None:
                    yield entry
                else:
                    yield entry.decode('utf-8')

    @staticmethod
    def to_dict(seq: str) -> Tuple[str, str]:
        name = seq[0].strip().split()[0][1:]
        return name, '\n'.join(seq)

    def parser(self):
        if not os.path.isfile('%s.shelve' % self.prefix):
            d = Parallel(n_jobs=self.cpus)(delayed(self.to_dict)(seq)
                                           for seq in tqdm(self.yielder(),
                                                           desc="Parsing file")
                                           )
            self.dictionary = dict(d)
            with shelve.open('%s.shelve' % self.prefix) as dic:
                dic.update(d)
        else:
            with shelve.open('%s.shelve' % self.prefix) as dic:
                self.dictionary = dict(dic)


def parse_pattern(pattern):
    if os.path.isfile(pattern):
        with open(pattern) as p:
            return p.read().strip().split('\n')
    else:
        return [pattern]


def print_it(key, dictionary):
    print(dictionary[key])


def main(filename: str, cpus: int, pattern: str, inverse: bool):
    fastq = FastQ(filename=filename, cpus=cpus)
    seqs = fastq.dictionary
    pattern = parse_pattern(pattern)
    if inverse:
        subset = set(seqs.keys()).difference(pattern)
    else:
        subset = set(seqs.keys()).intersection(pattern)

    _ = Parallel(n_jobs=cpus)(delayed(print_it)(key, seqs)
                              for key in tqdm(subset, total=len(subset),
                                              desc='Subsetting'))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('fastq', help="FastQ file with sequences")
    parser.add_argument('-c', '--cpus', help="Number of cpus to use", type=int,
                        default=-1)
    parser.add_argument('-v', '--inverse_match', action='store_true',
                        help="Return all matches not found in pattern")
    parser.add_argument('-p', '--pattern', action='store', required=True,
                        help="Pattern string or file to find")
    args = parser.parse_args()
    main(filename=args.fastq, cpus=args.cpus, inverse=args.inverse_match,
         pattern=args.pattern)