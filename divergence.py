import math
from subprocess import check_output
import sys
import pandas as pd


def calculate_jukes_cantor(sequence1, sequence2):
    '''
    Calculate the Jukes-Cantor distance between the two provided aligned
    sequences.

    Based on pseudocode provided on page 116 of the textbook.
    Taken from https://github.com/Xyaneon/jukes-cantor/blob/master/jc_algorithm.py
    '''
    # Initialization
    difference_counter = 0
    length_counter = 0
    # Step 1: Count differences between sequences, ignoring gaps
    for i in range(min(len(sequence1), len(sequence2))):
        if sequence1[i] != '-' and sequence2[i] != '-':
            length_counter += 1
            if sequence1[i] != sequence2[i]:
                difference_counter += 1
    # Step 2: Calculate and return results
    substitutions_fraction = float(difference_counter) / float(length_counter)
    jukes = -3.0 / 4.0 * math.log(1 - (4.0 / 3.0 * substitutions_fraction))
    return jukes

def parse_pair(o):
    q = o[o.find('\n'):o.rfind('>')].strip().replace('\n', '')
    a = o.split('\n>')[1]
    s = a[a.find('\n'):].strip().replace('\n', '')
    return q, s


def Parse_Fasta(filename):
    """
    Parse a single fasta file

    :param str filename: name of the fasta file
    """
    dic = {}
    name = None
    seq = ''
    with open(filename) as F:
        for line in F:
            if line.startswith('>'):
                if name is not None:
                    dic[name] = seq
                seq = ''
                name = line.strip()
            else:
                seq += line
        if not name in dic:
            dic[name] = seq
    return dic


table = sys.argv[1]
fasta = sys.argv[2]
df = pd.read_table(table, sep='\t')
d = Parse_Fasta(fasta)

for t in df.itertuples():
    try:
        quer = t.qaccver
        subj = t.saccver
        qseq = d['>%s' % quer]
        sseq = d['>%s' % subj]
        with open('temp.fas', 'w') as f:
            f.write('>%s\n%s>%s\n%s' % (quer, qseq, subj, sseq))
        o = check_output(['muscle', '-in', 'temp.fas'])
        dist = calculate_jukes_cantor(*parse_pair(o))
        with open('jc_difs.tsv', 'a') as f:
            f.write('%s\t%s\t%.3f\n' % (quer, subj, dist*100))
    except KeyError:
        continue



