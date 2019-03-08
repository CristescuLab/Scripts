import sys, os
import pandas as pd
from joblib import Parallel, delayed
import dill

def process_info(line, filter_exp):
    # get each field
    out = {}
    bl = line.split(';')
    for inf in bl:
        k, v = inf.split('=')
        out[k] = v
    df = pd.DataFrame([out], dtype=float)
    try:
        d = df.query(filter_exp)
        return d.shape[0] > 0
    except pd.core.computation.expr.UndefinedVariableError:
        return False


def single(line, filter_exp):
    if line.startswith('#'):
        return line
    else:
        bline = line.strip().split()
        if process_info(bline[7], filter_exp):
            return line


def main(fn, outfn, filter_exp):
    if not os.path.isfile('lines.pckl'):
        with open(fn) as fil, open('lines.pckl', 'wb') as pckl:
            outs = Parallel(n_jobs=-1, prefer='threads')(delayed(single)(
                line, filter_exp) for line in fil)
            dill.dump(outs, pckl)
    else:
        with open('lines.pckl', 'rb') as pckl:
            outs = dill.load(pckl)

    with open(outfn, 'w') as out:
        out.write('\n'.join(x.strip() for x in outs if x is not None))


if __name__ == '__main__':
    # Usage:
    #
    # python filter_variants.py vcf_file outfilename expression
    # stringent expression:
    # AF < 0.99 & SOR > 1.0 & ReadPosRankSum < 1.5 & ReadPosRankSum > 0.0 & \
    # MQRankSum > -0.5 & MQRankSum < 1.5 & MQ > 59.9 & MQ < 60.1 & FS < 6.0 & \
    # QD > 10.0

    main(sys.argv[1], sys.argv[2], sys.argv[3])
