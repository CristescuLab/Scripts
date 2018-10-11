import re
from itertools import product
import sys

iupac = {
    "R": "(A|G)",
    "Y": "(C|T)",
    "S": "(G|C)",
    "W": "(A|T)",
    "K": "(G|T)",
    "M": "(A|C)",
    "B": "(C|G|T)",
    "D": "(A|G|T)",
    "H": "(A|C|T)",
    "V": "(A|C|G)",
    "N": "(A|C|G|T)"
}

def expander(s):
    pat = r"\(([^)]*)\)"
    pieces = re.split(pat, s)
    pieces = [piece.split("|") for piece in pieces]
    for p in product(*pieces):
        yield "".join(p)

primers=sys.argv[1:]

done=[]
for primer in primers:
    for key, val in iupac:
        primer = primer.replace(key, val)
    done += [t for t in expander(primer)]

print('\n'.join(done))
