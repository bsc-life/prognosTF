import sys
from cPickle import load
from matplotlib import pyplot as plt
import numpy as np
from collections import OrderedDict




pickle_f = sys.argv[1]
size = int(sys.argv[2])

full_matrix  = load(open('data/matrix.pickle'))

resolution   = 10000
in_feature   = False
peak_files   = ['data/peaks_protA.bed', 'data/peaks_protB.bed']
chrom_sizes  = OrderedDict([('1', 500 * resolution + 1000),
                            ('2', 300 * resolution + 1000),
                            ('3', 200 * resolution + 1000)])
windows_span = size // 2


def read_line_feature(line):
    '''
    Get information per peak of a feature +/-
    '''
    c, p1, p2, f = line.split()[:4]
    return c, (int(p1) + int(p2)) / 2 / resolution, f

def read_line_no_feature(line):
    '''
    Get information per peak
    '''
    c, p1, p2 = line.split()[:3]
    return c, (int(p1) + int(p2)) / 2 / resolution, ''

def read_line_no_feature_but(line):
    '''
    Get information per peak
    '''
    c, p1, p2 = line.split()[:3]
    return c, (int(p1) + int(p2)) / 2 / resolution, '{}:{}-{}'.format(c, p1, p2)

peaks1 = open(peak_files[0], "r")

try:
    peaks2 = open(peak_files[1], "r")
    same = False
except IndexError:
    same = True
    peaks2 = peaks1

# findout if bed file contain features, or only coordinates
line = next(peaks1)
try:
    read_line_feature(line)
    read_line1 = read_line_feature
except ValueError:
    if in_feature:
        read_line1 = read_line_no_feature_but
    else:
        read_line1 = read_line_no_feature
peaks1.seek(0)

line = next(peaks2)
try:
    read_line_feature(line)
    read_line2 = read_line_feature
except ValueError:
    read_line2 = read_line_no_feature

max_chrom = dict((c, chrom_sizes[c] // resolution - windows_span)
                 for c in chrom_sizes)

peaks1.seek(0)  # needs to be here in case peaks1 and peak2 are the same
bin_coordinate1 = set((c, p, f) for c, p, f in map(read_line1, peaks1)
                      if c in max_chrom and windows_span <= p <= max_chrom[c])

peaks2.seek(0)  # needs to be here in case peaks1 and peak2 are the same
bin_coordinate2 = set((c, p, f) for c, p, f in map(read_line2, peaks2)
                      if c in max_chrom and windows_span <= p <= max_chrom[c])

peaks1.seek(0)
npeaks1 = sum(1 for _ in peaks1)
peaks2.seek(0)
npeaks2 = sum(1 for _ in peaks2)

# sort peaks
bin_coordinate1 = sorted(bin_coordinate1)
bin_coordinate2 = sorted(bin_coordinate2)

print(bin_coordinate1)
print(bin_coordinate2)

sections = {}
total = 0
for c in chrom_sizes:
    sections[c] = total
    total += chrom_sizes[c] // resolution + 1
print (sections)


real_matrix = [[0 for i in range(size)] for j in range(size)]
count = 0
for c1, b1, _ in bin_coordinate1:
    i = sections[c1] + b1
    for c2, b2, _ in bin_coordinate2:
        if c1 != c2:
            continue
        j = sections[c2] + b2
        if abs(j - i) >= windows_span:
            count += 1
            for k in range(size):
                for l in range(size):
                    real_matrix[k][l] += full_matrix[i + k - windows_span][j + l - windows_span]
print (count)


coord, (sum_raw, sqr_raw, sum_nrm, sqr_nrm, passage, tot) = load(open(pickle_f))

matrix = [[sum_raw[(i, j)] for i in range(size)] for j in range(size)]

plt.subplot(221)
plt.imshow(np.log2(matrix), origin='lower')
for i in range(len(matrix)):
    for j in range(len(matrix)):
        if matrix[i][j]:
            plt.text(j, i, matrix[i][j], ha='center', va='center', size=6)

plt.colorbar()

matrix = [[passage[(i, j)] for i in range(size)] for j in range(size)]

plt.subplot(222)
plt.imshow(np.log2(matrix), origin='lower')
for i in range(len(matrix)):
    for j in range(len(matrix)):
        if matrix[i][j]:
            plt.text(j, i, matrix[i][j], ha='center', va='center', size=6)

plt.colorbar()

plt.subplot(223)
plt.imshow(np.log2(real_matrix), origin='lower')
for i in range(len(real_matrix)):
    for j in range(len(real_matrix)):
        if real_matrix[i][j]:
            plt.text(j, i, real_matrix[i][j], ha='center', va='center', size=6)

plt.colorbar()

plt.show()
