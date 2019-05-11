"""
**pythonblast.py
** Copyright (C) 2018  Jose Sergio Hleap

Python wrapper for multithreaded blast one sequence at a time, it uses shelve,
a persistent dictionary-like object, to avoid overusing memory as well.

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
from io import BytesIO
import pandas as pd
from joblib import Parallel, delayed
from subprocess import run, PIPE, CalledProcessError
import shelve
import dill


def parse_fasta(filename, fn2=None, rename=False):
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
    dic = shelve.open(db)
    name = None
    seq = ''
    with fnc(filename) as F:
        for line in F:
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


def stdin_run(args, inpt, **kwargs):
    file_name = None
    inpt = inpt.encode('utf-8') if isinstance(inpt, str) else input
    if isinstance(input, tuple):
        file_name, inpt = inpt
    exe = run(args, input=inpt, stderr=PIPE, stdout=PIPE, **kwargs)
    try:
        exe.check_returncode()
    except CalledProcessError:
        raise Exception(exe.stdout, exe.stderr)
    return exe.stdout, file_name


def iter_fasta(fn, done=[]):
    fastas = shelve.open(fn)
    for header, sequence in fastas.items():
        if header not in done:
            done.append(header)
            yield '%s\n%s' % (header, sequence)


def main(db, query, evalue, p_id,  max_target_seqs=50, cpus=-1, out='hit.hits',
         outfmt_str='qseqid sseqid pident evalue qcovs qlen length staxid '
                    'stitle'):
    fasta = parse_fasta(query)
    args = ['blastn', '-db', db, '-query', '-', '-evalue', evalue,
            '-perc_identity', p_id, '-outfmt', '"6 %s"' % outfmt_str,
            '-max_target_seqs', max_target_seqs]
    blasts = Parallel(n_jobs=cpus)(
        delayed(stdin_run)(args, inp) for inp in iter_fasta(fasta))
    blasts = [pd.read_table(BytesIO(x), header=None, names=outfmt_str.split())
              for x in blasts]
    blasts = pd.concat(blasts)
    blasts.to_csv(out, sep='\t', index=False, header=False)
    return blasts


if __name__ == '__main__':
    import optparse
    opts = optparse.OptionParser(usage='%prog [options] fasta db')
    opts.add_option('--output', '-o', default='hit.hits',
                    help='Output filename [default: %default]')
    opts.add_option('--pident', '-p', action='store', type=float, default=None,
                    help='Minimum percent identity [default: %default]')
    opts.add_option('--eval', '-e', action='store', type=float, default=None,
                    help='Maximum evalue [default: %default]')
    opts.add_option('--max_targets', '-m', action='store', type=int,
                    default=50, help='Max target sequences [default: %default]'
                    )
    opts.add_option('--cpus', '-c', action='store', default=-1, type=int,
                    help='number of cpus to use in the run [default: %default]'
                    )
    opts.add_option('--outfmt', '-f', action='store',
                    default='qseqid sseqid pident evalue qcovs qlen length '
                            'staxid stitle', help='outfmt [default: %default]'
                    )
    opt, arg = opts.parse_args()
    main(arg[1], arg[0], opt.evalue, opt.pident, opt.max_targets, opt.cpus,
         opt.output, opt.outfmt)
