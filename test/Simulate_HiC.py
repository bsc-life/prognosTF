#!/usr/bin/env python


from random       import random, seed
from collections  import OrderedDict
from fractions    import gcd
from subprocess   import Popen
from functools    import reduce

import numpy as np
from numpy.random import negative_binomial

from matplotlib   import pyplot as plt


def load_genome(chroms):
    sections = {}
    total = 0
    for c in chroms:
        sections[c] = total
        total += chroms[c]
    bins = {}
    for c in chroms:
        for i in range(chroms[c]):
            bins[sections[c] + i] = c, i
    d = reduce(lambda a, b: gcd(a, b), chroms.values())
    weighted_chroms = [c for c in chroms for _ in range(int(chroms[c]**1.5 / d))]
    return sections, bins, weighted_chroms

def main():
    seed_num = 1
    reso = 10000
    chroms = OrderedDict([('1', 500), ('2', 300), ('3', 200)])
    genome_size = sum(chroms.values())

    sections, bins, weighted_chroms = load_genome(chroms)

    seed(seed_num)
    np.random.seed(seed_num)

    npeaks = 40
    cmprts_pos = {}
    bad_cols = {}
    prob = 0.2
    step = 10
    for c in chroms:
        bad_cols[c] = set()
        for _ in range(chroms[c] // 10):
            bad_cols[c].add(int(random() * chroms[c]))
        cmprts_pos[c] = []
        end = 0
        beg = 0
        while end < chroms[c]:
            if random() < prob:
                cmprts_pos[c].append((beg, end))
                beg = end
            end += step
        cmprts_pos[c].append((beg, end))
    cmprts = {}
    for c in cmprts_pos:
        cmprts[c] = {'A': set(p for i, (beg, end) in enumerate(cmprts_pos[c])
                              for p in range(beg, end) if i %2),
                     'B': set(p for i, (beg, end) in enumerate(cmprts_pos[c])
                              for p in range(beg, end) if not i %2)}

    peaks = set()
    peaks1 = set()
    peaks2 = set()
    for c in range(npeaks):
        bin1 = int(random() * (genome_size - 2))
        if random() < 0.4:
            peaks1.add(bin1)
        else:
            peaks2.add(bin1)
        peaks.add(bin1)

    loops = set()
    for bin1 in peaks:
        for bin2 in peaks:
            if random() < 0.1:
                continue
            if bin1 in peaks1:
                range1 = 3
            else:
                range1 = 2
            if bin2 in peaks1:
                range2 = 3
            else:
                range2 = 2
            for i in range(range1):
                for j in range(range2):
                    loops.add((bin1 + i, bin2 + j))
                    loops.add((bin2 + j, bin1 + i))

    print('generating SAM')
    out = open('data/fake.sam', 'w')
    out.write('@HD\tVN:1.5\tSO:coordinate\n')
    for c in chroms:
        out.write('@SQ\tSN:%s\tLN:%d\n' % (c, chroms[c] * reso - 1))
    matrix = [[0 for _ in range(sum(chroms.values()))]
              for _ in range(sum(chroms.values()))]
    nrnd = 1000000
    bin_prob = 0.005
    nbs = iter(negative_binomial(1, bin_prob, size=nrnd))
    count = 0
    while count < nrnd:
        c1 = weighted_chroms[int(random() * len(weighted_chroms))]
        pos1 = int(random() * chroms[c1])
        if random() > (float(chroms[c1]) / genome_size)**0.8:
            c2 = weighted_chroms[int(random() * len(weighted_chroms))]
            pos2 = int(random() * chroms[c2])
        else:
            c2 = c1
            pos2 = -1
            while pos2 < 0 or pos2 >= chroms[c2]:
                try:
                    if random() < 0.15:
                        wanted_cmprt = 'A' if pos1 in cmprts[c1]['A'] else 'B'
                        while pos2 not in cmprts[c2][wanted_cmprt]:
                            pos2 = pos1 + (next(nbs) * (-1 if random() > 0.5 else 1))
                    else:
                        pos2 = pos1 + (next(nbs) * (-1 if random() > 0.5 else 1))
                except StopIteration:
                    nbs = iter(negative_binomial(1, bin_prob, size=nrnd))
        if pos1 in bad_cols[c1] or pos2 in bad_cols[c2]:
            if random() < 0.5:
                continue
        bin1 = sections[c1] + pos1
        bin2 = sections[c2] + pos2
        if random() < 0.5:
            if (bin1, bin2) not in loops:
                continue
        out.write('SRR.{0}\t1024\t{1}\t{2}\t1\t75P\t{3}\t{4}\t75\t*\t*\n'.format(
            i, c1, int(reso / 2 + pos1 * reso), c2, int(reso / 2 + pos2 * reso)))
        out.write('SRR.{0}\t1024\t{1}\t{2}\t1\t75P\t{3}\t{4}\t75\t*\t*\n'.format(
            i, c2, int(reso / 2 + pos2 * reso), c1, int(reso / 2 + pos1 * reso)))
        matrix[bin1][bin2] += 1
        matrix[bin2][bin1] += 1
        count += 1
    out.close()

    print('generating BAM')
    Popen('samtools sort -@ 8 -O BAM {} > {}'.format('data/fake.sam',
                                                     'data/fake.bam'),
          shell=True).communicate()
    Popen('rm -f {}'.format('data/fake.sam'), shell=True).communicate()
    Popen('samtools index -@ 8 {}'.format('data/fake.bam'),
          shell=True).communicate()

    # print('tadbit normalize -w tmp --bam {} -r {}'.format('data/fake.bam', reso)).communicate()
    Popen('tadbit normalize -w tmp --bam {} -r {}'.format('data/fake.bam', reso), shell=True).communicate()
    Popen('mv tmp/04_normalization/biases* data/biases.pickle', shell=True).communicate()
    Popen('rm -rf tmp', shell=True).communicate()

    plt.figure(figsize=(41, 30))
    plt.imshow(np.log2(matrix), interpolation='None', origin='lower')
    plt.colorbar()
    plt.savefig('data/matrix.png')

    print('Saving BEDs')
    out = open('data/peaks_protA.bed', 'w')
    out.write(''.join('{0}\t{1}\t{2}\n'.format(
        bins[p][0],
        bins[p][1] * reso + int(random() * 1000),
        bins[p][1] * reso + int(random() * 1000)) for p in sorted(peaks1)))
    out.close()

    out = open('data/peaks_protB.bed', 'w')
    out.write(''.join('{0}\t{1}\t{2}\n'.format(
        bins[p][0],
        bins[p][1] * reso + int(random() * 1000),
        bins[p][1] * reso + int(random() * 1000)) for p in sorted(peaks2)))
    out.close()

    out = open('data/compartments.bed', 'w')
    out.write(''.join('{}\t{}\t{}\t{}\n'.format(
        c, p * reso, p * reso + reso,
        (1 if p in cmprts[c]['A'] else -1) * (0.2 + 0.8 * random()))
                      for c in chroms
                      for p in range(chroms[c])))
    out.close()


if __name__ == "__main__":
    exit(main())
