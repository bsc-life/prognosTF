import sys
from cPickle import load
from matplotlib import pyplot as plt
import numpy as np


pickle_f = sys.argv[1]
size = int(sys.argv[2])

coord, (sum_raw, sqr_raw, sum_nrm, sqr_nrm, passage, tot) = load(open(pickle_f))

matrix = [[sum_raw[(i, j)] for i in range(size)] for j in range(size)]

plt.subplot(121)
plt.imshow(np.log2(matrix), origin='lower')
for i in range(len(matrix)):
    for j in range(len(matrix)):
        if matrix[i][j]:
            plt.text(j, i, matrix[i][j])

plt.colorbar()

matrix = [[passage[(i, j)] for i in range(size)] for j in range(size)]

plt.subplot(122)
plt.imshow(np.log2(matrix), origin='lower')
for i in range(len(matrix)):
    for j in range(len(matrix)):
        if matrix[i][j]:
            plt.text(j, i, matrix[i][j])

plt.colorbar()

plt.show()
