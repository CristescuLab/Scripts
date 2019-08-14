from tqdm import tqdm
import pandas as pd
import shelve
import sys

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


def main(filtered, fasta, taxa):
    taxa = taxa.lower()
    prefix = fasta[:fasta.rfind('.fasta')]
    fas = parse_fasta(fasta)
    df = pd.read_csv(filtered, sep='\t')
    with shelve.open(fas) as dic:
        seqs = [x[1:] for x in list(dic.keys())]
        sdf = df[df.qseqid.isin(seqs)]
        for name, d in tqdm(sdf.groupby(taxa), desc="Getting %s" % taxa):
            if d.empty:
                continue
            with open('%s_%s.fas' % (prefix, name.replace(' ', '_')), 'w'
                       ) as out:
                q = d.qseqid.unique().tolist()
                for seq in q:
                    out.write('>%s\n%s' % (seq, dic['>%s' % seq]))


if __name__ == '__main__':
    # Usage:
    # python get_sequences.py filtered_blast fasta taxa
    main(sys.argv[1], sys.argv[2], sys.argv[3])