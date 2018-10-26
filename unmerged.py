"""
**unmergerged.py
** Copyright (C) 2018  Jose Sergio Hleap

Given a pairend sequencing, if the amplicon is too long for the reads,
blast each one separately and report the matches nhgc

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
"""
from io import StringIO
import pandas as pd
from multiprocessing import Pool
from joblib import Parallel, delayed
from subprocess import Popen, PIPE
import optparse
from itertools import zip_longest
import gzip
import csv
import os


names = 'qseqid sseqid pident evalue qcovs qlen length staxid'.split()

def iterzip(a, b, kwargs=None):
    """
    zip two iterators (assumes the same length)
    :param a: iterator one
    :param b: iterator 2
    :return: tuple of elements
    """
    for e_a in a:
        e_b = next(b)
        assert e_a['id'].strip().split()[0] == e_b['id'].strip().split()[0]
        yield e_a, e_b, kwargs


def grouper(iterable, n, fillvalue=None):
    """
    Collect data into fixed-length chunks or blocks
    From itertools recipe

    :param iterable: container to be iterate over
    :param n: Number of items to return
    :param fillvalue: Padding
    :return: zip_longest iterator
    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


def read_fastq(fn):
    """
    Iterator over fastq files

    :param fn: filename
    :return: dictionary
    """
    fnc = gzip.open if 'gz' in fn else open
    fields = ['id', 'seq', 'plus', 'qual']
    with fnc(fn) as f:
        for lines in grouper(f, 4, ''):
            assert len(lines) == 4
            yield dict(zip(fields,lines))


def blast_helper(idx, r):
    r1, r2, kwargs = r
    db, evalue, p_id, tgts = kwargs
    blast_line = ['blastn', '-db', db, '-evalue', str(evalue),
                  '-perc_identity', str(p_id), '-max_target_seqs', str(tgts),
                  '-outfmt', '6 qseqid sseqid pident evalue qcovs qlen length '
                             'staxid']
    seq1 = '>read%d\n%s' % (idx, r1['seq'].decode('utf-8').strip())
    seq2 = '>read%d\n%s' % (idx, r2['seq'].decode('utf-8').strip())
    mapreads = '%s\t%s\tread%d' % (r1['id'].strip(), r2['id'].strip(), idx)
    # run blasts
    p1 = Popen(blast_line, stdin=PIPE, stdout=PIPE, stderr=PIPE,
              env=os.environ.copy(), encoding='utf8', cwd=os.getcwd())
    o1, e1 = p1.communicate(r1['seq'].decode('utf-8').strip())
    p2 = Popen(blast_line, stdin=PIPE, stdout=PIPE, stderr=PIPE,
              env=os.environ.copy(), encoding='utf8', cwd=os.getcwd())
    o2, e2 = p2.communicate(r2['seq'].decode('utf-8').strip())
    # read outputs
    df1 = pd.read_table(StringIO(o1), sep='\t', header=None, names=names,
                        comment='#', quoting=csv.QUOTE_NONE, encoding='utf-8')
    df2 = pd.read_table(StringIO(o2), sep='\t', header=None, names=names,
                        comment='#', quoting=csv.QUOTE_NONE, encoding='utf-8')
    # check taxids
    b = df2.staxid.isin(df1.staxid)
    c = df1.staxid.isin(df2.staxid)
    if b.any():
        concordant = df1[c].merge(df2[b], on=['qseqid', 'staxid'],
                                  suffixes=('_R1', '_R2'),)
        discordant = None
        c_appended = '>seq%s_appended\n%s%s' % (idx, r1['seq'].decode(
            'utf-8').strip(), r2['seq'].decode('utf-8').strip())
        d_appended = None
    else:
        df1['qseqid'] = df1.qseqid.apply(lambda x: '%s_R1' % x)
        df2['qseqid'] = df2.qseqid.apply(lambda x: '%s_R2' % x)
        discordant = pd.concat([df1, df2])
        concordant = None
        c_appended = None
        d_appended = (seq1, seq2)
    return concordant, discordant, mapreads, c_appended, d_appended


def blast_pair(prefix, db, evalue=1E-10, p_id=95, tgts=50, cpus=-1):
    cpus = os.cpu_count() if cpus == -1 else cpus
    fwd = '%s_R1.fastq.gz' % prefix
    rev = '%s_R2.fastq.gz' % prefix
    kwargs = (db, evalue, p_id, tgts)
    it = iterzip(read_fastq(fwd), read_fastq(rev), kwargs)
    with Pool(processes=cpus) as pool:
        out = pool.starmap(blast_helper, enumerate(it))
    #out = Parallel(n_jobs=cpus)(delayed(blast_helper)(
    #    r[0], r[1], idx, db, evalue, p_id, tgts) for idx, r in enumerate(it))
    concordant, discordant, mapreads, c_appended, d_appended  = zip(*out)
    pd.concat(concordant).to_csv('%s_concordantblast.tsv' % prefix, sep='\t',
                                 index=False)
    pd.concat(concordant).to_csv('%s_discordantblast.tsv' % prefix, sep='\t',
                                 index=False)
    fwd, rev = zip(*[x for x in d_appended if x is not None])
    with open('%s_mapreads.txt' % prefix, 'w') as m, \
            open('%s_appended.fasta' % prefix, 'w') as a, \
            open('%s_discordant_R1.fasta' % prefix, 'w') as r1, \
            open('%s_discordant_R2.fasta' % prefix, 'w') as r2:
        m.write('R1\tR2\tread\n%s' % '\n'.join(mapreads))
        a.write('\n'.join([x for x in c_appended if x is not None]))
        r1.write('\n'.join(fwd))
        r2.write('\n'.join(rev))


if __name__ == '__main__':
    intro = '''python %prog fastq_prefix db [options]'''

    opts = optparse.OptionParser(usage='%s' % intro)
    opts.add_option('--pident', '-p', action='store', type=float, default=95,
                    help='Minimum percent identity [default: %default]')
    opts.add_option('--eval', '-e', action='store', type=float, default=1E-10,
                    help='Maximum evalue [default: %default]')
    opts.add_option('--max_target_seqs', '-m', action='store', type=int,
                    default=50, help='Maximum number of targets '
                                     '[default: %default]')
    opts.add_option('--threads', '-t', action='store', type=int,
                    default=-1, help='threads to use [default: %default]')

    opt, arg = opts.parse_args()
    blast_pair(arg[0], arg[1], evalue=opt.eval, p_id=opt.pident,
               tgts=opt.max_target_seqs, cpus=opt.threads)