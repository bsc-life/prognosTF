import sys
try:
    from cPickle import load
except ImportError:
    from pickle import load
from matplotlib import pyplot as plt
import numpy as np


pickle_f = sys.argv[1]
size = int(sys.argv[2])

groups = load(open(pickle_f, 'rb'))

group = list(groups.keys())[0]


matrix = [[groups[group]['sum_raw'][(i, j)] for i in range(size)] for j in range(size)]
plt.subplot(221)
plt.imshow(np.log2(matrix), origin='lower')
for i in range(len(matrix)):
    for j in range(len(matrix)):
        if matrix[i][j]:
            plt.text(j, i, matrix[i][j], ha='center', va='center', size=6)
plt.colorbar()

matrix = [[groups[group]['sum_nrm'][(i, j)] / groups[group]['counter'] for i in range(size)] for j in range(size)]
plt.subplot(222)
plt.imshow((matrix), origin='lower')
plt.colorbar()

matrix = [[groups[group]['sum_nrm'][(i, j)] / groups[group]['passage'][(i, j)] for i in range(size)] for j in range(size)]
plt.subplot(223)
plt.imshow((matrix), origin='lower')
plt.colorbar()

matrix = [[groups[group]['passage'][(i, j)] for i in range(size)] for j in range(size)]
plt.subplot(224)
plt.imshow(np.log2(matrix), origin='lower')
for i in range(len(matrix)):
    for j in range(len(matrix)):
        if matrix[i][j]:
            plt.text(j, i, matrix[i][j], ha='center', va='center', size=6)
plt.colorbar()

plt.show()
