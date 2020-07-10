"""
"""

from pickle import load, dump
from collections import defaultdict
import sys

in_file = sys.argv[1]
try:
    ou_file = sys.argv[2]
except IndexError:
    ou_file = in_file

ori = load(open(in_file, 'rb'))

dico = dict()
size = resolution = None

key = next(iter(ori))
size       = ori[key]['size']
resolution = ori[key]['resolution']

for key, val in ori.items():
    dico[key] = defaultdict(int)
    dico[key].update((pos, (val['passage'][pos],
                            int(round(val['sum_nrm'][pos] * 1000, 0)),
                            int(round(val['sqr_nrm'][pos] * 1000, 0))))
                     for pos in val['sum_raw'])
    dico[key]['counter'] = val['counter']
    if size != val['size'] or resolution != val['resolution']:
        print (size, val['size'], resolution, val['resolution'])
        raise Exception('differing size or resoluton')

dico['size']       = size
dico['resolution'] = resolution

out = open(ou_file, 'wb')
dump(dico, out)
out.close()
