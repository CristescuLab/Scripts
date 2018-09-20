from joblib import Parallel, delayed
from itertools import zip_longest
from datetime import datetime
from tqdm import tqdm
import requests
import json
import time
import sys

SIX = ['phylum', 'class', 'order', 'family', 'genus', 'species']
kingdom=['Eukaryota']
number_of_lines = 25 if len(sys.argv) < 4 else int(sys.argv[4])


def get_batch(line):
    """
    Get a set of lines and get the json file
    :param line: bold accessions to process
    :return: tsv string with accession and lineage
    """
    def loop(r):
        try:
            j = r.json()
            records = j['bold_records']['records']
            tsv=''
            for acc in records:
                tax = records[acc]['taxonomy']
                lineage = ';'.join(kingdom + [
                    tax[level]['taxon']['name'] if level in tax else ''
                    for level in SIX])
                tsv += '%s\t%s\n' % (acc, lineage)
            return tsv
        except json.decoder.JSONDecodeError as err:
            print(err)
            print('Decode error with status', r.status_code)
            with open('failed.dump', 'a') as f:
                f.write('\n'.join(line.split('|')))
    #print('Processing from', line[0], 'to', line[-1])
    line = '|'.join([x.strip() for x in line if x is not None])
    url = 'http://www.boldsystems.org/index.php/API_Public/specimen'
    payload = dict(ids=line, format='json')
    r = requests.get(url, params=payload)
    status = r.status_code
    if status == 200:
        tsv = loop(r)
    else:
        print('Failed with status', status, 'showing', r.text )
        print('Trying up to 4 times or dumping')
        count = 0
        while status != 200 or count <= 4:
            r = requests.get(url, params=payload)
            status = r.status_code
            time.sleep(0.1)
            count += 1
        if count <= 4:
            tsv = loop(r)
        else:
            print('Request failed, dumping failed accessions to failed.dump')
            with open('failed.dump', 'a') as f:
                f.write('%s\n' % '\n'.join(line.split('|')))
            tsv = None
    return tsv


def grouper(iterable, n):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=None)


def wc(fn):
    with open(fn) as f:
        return sum(1 for _ in f)


ifn = sys.argv[1]
ofn = sys.argv[2]
cpus = int(sys.argv[3])
n_lines = wc(ifn)
iterations = int(n_lines / number_of_lines)
with open(ifn, 'r') as input_file, open(ofn, 'w') as output_file :
    tsv = Parallel(n_jobs=cpus)(delayed(get_batch)(
        current_line) for current_line in tqdm(grouper(
        input_file, number_of_lines), total=iterations))
    output_file.write(''.join([x for x in tsv if x is not None]))
