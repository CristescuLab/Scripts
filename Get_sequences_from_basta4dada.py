from tqdm import tqdm
import pandas as pd
import shelve
import sys
import os

def parse_fasta(filename):
    """
    Parse a fasta file and put it in a shelve DB and rename it with the prefix
    of the filename. It will also write a combined fasta with the renamed
    entries

    :param fil: name of fasta files
    :param fn2: name of the shelve
    :return: shelve db name
    """
    fn2 = '%s.shelve' % filename[:filename.rfind('.fasta')]
    with shelve.open(fn2) as dic:
        name = None
        seq = ''
        with open(filename) as F:
            for line in tqdm(F, desc="Parsing %s" % filename):
                if line.startswith('>'):
                    if name is not None:
                        dic[name] = seq
                    seq = ''
                    name = '%s' % (line.strip())
                else:
                    seq += line
            if name not in dic:
                dic[name] = seq
    return fn2


def main(basta, fasta, otutable,  outdir):
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    df = pd.read_csv(otutable, sep='\t')
    samples = df.index.tolist()
    basta = pd.read_csv(basta, sep='\t', header=None)
    basta['assigned'] = basta[1].apply(
        lambda x: x.strip().strip(';').split(';')[-1])
    fas = parse_fasta(fasta)
    with shelve.open(fas) as dic:
        for row in samples:
            seqs = df.loc[row,:][df.loc[row,:] > 0]
            sdf = basta[basta[0].isin(seqs.index.tolist())]
            for name, d in tqdm(sdf.groupby('assigned')):
                path = os.path.join(outdir, row, name.replace(' ', '_'))
                if d.empty:
                    continue
                if not os.path.isdir(os.path.dirname(path)):
                    os.mkdir(os.path.dirname(path))
                with open('%s.fas' % path, 'w') as out:
                    q = d[0].unique().tolist()
                    for seq in q:
                        out.write('>%s\n%s' % (seq, dic['>%s' % seq]))


if __name__ == '__main__':
    # Usage:
    # python Get_sequences_from_basta4dada.py bastaout allfasta otutable outdir
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


