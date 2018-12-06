"""
**pythonblast.py
** Copyright (C) 2018  Jose Sergio Hleap

Python wrapper for multithreaded blast one sequence at a time

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
from multiprocessing import Pool
from joblib import Parallel, delayed
from subprocess import run, PIPE, CalledProcessError
import shelve

def parse_fasta(filename, fn2=None):
    """
    Parse a single fasta file

    :param str filename: name of the fasta file
    """
    if isinstance(filename, shelve.DbfilenameShelf):
        Warning('THIS IS ALREADY A SHELF!!!')
        return filename
    elif isinstance(filename, bytes):
        fnc = BytesIO
        fn = 'shelve' if fn2 is None else fn2
    else:
        fnc = open
        fn = filename[:filename.find('.')]
    dic = shelve.open('%s.shelve' % fn)
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
    return dic


def stdin_run(args, input, **kwargs):
    file_name = None
    input = input.encode('utf-8') if isinstance(input, str) else input
    if isinstance(input, tuple):
        file_name, input = input
    exe = run(args, input=input, stderr=PIPE, stdout=PIPE,
              **kwargs)
    try:
        exe.check_returncode()
    except CalledProcessError:
        raise Exception(exe.stdout, exe.stderr)
    return exe.stdout, file_name


def iter_fasta(fn):
    fastas = parse_fasta(fn)
    for header, sequence in fastas.items():
            yield '%s\n%s' % (header, sequence)


def main(db, query, evalue, p_id,  max_target_seqs=50, cpus=-1,
         outfmt_str='qseqid sseqid pident evalue qcovs qlen length staxid '
                    'stitle'):
    fasta = parse_fasta(query)
    args = ['blastn', '-db', db, '-', query, '-evalue', evalue,
            '-perc_identity', p_id, '-outfmt', '"6 %s"' % outfmt_str,
            '-max_target_seqs', max_target_seqs]
    blasts = Parallel(n_jobs=cpus)(
        delayed(stdin_run)(args, input) for input in iter_fasta(fasta))
    blasts = [pd.read_table(BytesIO(x), header=None, names=outfmt_str.split())
              for x in blasts]
    blasts = pd.concat(blasts)
    blasts.to_csv()