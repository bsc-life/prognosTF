#! /usr/bin/env python

"""
"""

from collections import defaultdict


def parse_peaks(peak_files, resolution, in_feature, chrom_sizes, windows_span):

    def read_line_feature(line):
        '''
        Get information per peak of a feature +/-
        '''
        c, p1, p2, f = line.split()[:4]
        return c, (int(p1) + int(p2)) // 2 // resolution, f

    def read_line_no_feature(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        return c, (int(p1) + int(p2)) // 2 // resolution, ''

    def read_line_no_feature_but(line):
        '''
        Get information per peak
        '''
        c, p1, p2 = line.split()[:3]
        return c, (int(p1) + int(p2)) // 2 // resolution, '{}:{}-{}'.format(c, p1, p2)

    peaks1 = open(peak_files[0], "r")
    try:
        peaks2 = open(peak_files[1], "r")
    except IndexError:
        peaks2 = peaks1

    # find out if bed file contain features, or only coordinates
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
    if len(peak_files) == 1:
        bin_coordinate2 = bin_coordinate1
    else:
        bin_coordinate2 = sorted(bin_coordinate2)

    peaks1.close()
    if len(peak_files) > 1:
        peaks2.close()

    return bin_coordinate1, bin_coordinate2, npeaks1, npeaks2


def generate_pairs(bin_coordinate1, bin_coordinate2, resolution, windows_span,
                   max_dist, window, section_pos):

    wsp = (windows_span * 2) + 1
    mdr = max_dist / resolution

    # put pairs in intervals
    if window == 'inter':
        test = lambda a, b: (a[0] != b[0]
                             and a != b)
    elif window == 'intra':
        test = lambda a, b: (a[0] == b[0]
                             and wsp <= abs(b[1] - a[1]) <= mdr
                             and a != b)
    elif window == 'all':
        test = lambda a, b: ((a[0] == b[0]
                              and wsp <= abs(b[1] - a[1]) <= mdr
                              and a != b) or (a[0] != b[0] and a != b))
    else:
        lower, upper = window
        test = lambda a, b: (a[0] == b[0]
                             and wsp <= abs(b[1] - a[1]) <= mdr
                             and a != b
                             and lower < abs(a[1] - b[1]) <= upper)

    if bin_coordinate1 is bin_coordinate2:  # we want only one side
        pairs = ((a, b) for i, a in enumerate(bin_coordinate1, 1)
                 for b in bin_coordinate2[i:]
                 if test(a, b))
    else:
        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if test(a, b))

    # Sort pairs of coordinates according to genomic position of the
    # smallest of each pair, and store it into a new list
    final_pairs = []
    for (chr1, bs1, f1), (chr2, bs2, f2) in pairs:
        pos1 = section_pos[chr1][0] + bs1
        pos2 = section_pos[chr2][0] + bs2

        beg1 = pos1 - windows_span
        end1 = pos1 + windows_span + 1
        beg2 = pos2 - windows_span
        end2 = pos2 + windows_span + 1

        if beg1 > beg2:
            beg1, end1, beg2, end2 = beg2, end2, beg1, end1

        what = f1 + f2
        final_pairs.append((beg1, end1, beg2, end2, what))

    final_pairs.sort()

    return final_pairs


def submatrix_coordinates(final_pairs, badcols, wsp, counter):
    '''
    Input BED file(s) of ChIP peaks and bin into desired resolution of Hi-C
    '''

    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features


    # in buf we store a list of coordinates to be yielded
    # when buf spans for twice the window span we sort it and empty it
    buf = []
    buf_beg = 0
    for beg1, end1, beg2, end2, what in final_pairs:
        range1 = [(x, p1) for x, p1 in enumerate(range(beg1, end1))
                  if p1 not in badcols]
        range2 = [(y, p2) for y, p2 in enumerate(range(beg2, end2))
                  if p2 not in badcols]

        if not range1 or not range2:
            continue
        counter[what] +=1
        for x, p1 in range1:
            for y, p2 in range2:
                buf.append(((p1, p2), x, y, what))

        if end1 - buf_beg > wsp:
            buf.sort()
            top = end1 - wsp
            p1 = min(buf)[0][0]
            while p1 < top:
                (p1, p2), x, y, what = buf.pop(0)
                yield (p1, p2), x, y, what
            buf_beg = p1
    buf.sort()
    while buf:
        yield buf.pop(0)


def readfiles(genomic_file, iter_pairs):
    # create empty meta-waffles
    fh1 = open(genomic_file)
    a, b, raw, nrm = next(fh1).split('\t')
    pos1 = (int(a), int(b))
    pos2, x, y, group = next(iter_pairs)

    try:
        while True:
            if pos2 > pos1:
                a, b, raw, nrm = next(fh1).split('\t')
                pos1 = (int(a), int(b))
            elif pos1 == pos2:
                raw = int(raw)
                nrm = float(nrm)
                yield pos1, x, y, raw, nrm, group
                pos2, x, y, group = next(iter_pairs)
                if pos1 != pos2:  # some cells in the peak file are repeated
                    a, b, raw, nrm = next(fh1).split('\t')
                    pos1 = (int(a), int(b))
            else:
                pos2, x, y, group = next(iter_pairs)
    except StopIteration:
        fh1.close()


def interactions_at_intersection(genomic_mat, iter_pairs, submatrices, bins):
    def _write_submatrices(X, Y, x, y, raw, nrm, group):
        c1, b1 = bins[X]
        c2, b2 = bins[Y]
        out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            c1, b1, c2, b2, x, y, raw, nrm, group))

    if bins:
        out = open(submatrices, 'w')
        do_the_thing = _write_submatrices
    else:
        do_the_thing = lambda a, b, c, d, e, f, g: None

    groups = {}
    readfiles_iterator = readfiles(genomic_mat, iter_pairs)
    for (X, Y), x, y, raw, nrm, group in readfiles_iterator:
        try:
            groups[group]['sum_raw'][x, y] += raw
            groups[group]['sqr_raw'][x, y] += raw**2
            groups[group]['sum_nrm'][x, y] += nrm
            groups[group]['sqr_nrm'][x, y] += nrm**2
            groups[group]['passage'][x, y] += 1
            do_the_thing(X, Y, x, y, raw, nrm, group)
        except KeyError:
            groups[group] = {
                'sum_raw' : defaultdict(int),
                'sqr_raw' : defaultdict(int),
                'sum_nrm' : defaultdict(float),
                'sqr_nrm' : defaultdict(float),
                'passage' : defaultdict(int)}
            groups[group]['sum_raw'][x, y] += raw
            groups[group]['sqr_raw'][x, y] += raw**2
            groups[group]['sum_nrm'][x, y] += nrm
            groups[group]['sqr_nrm'][x, y] += nrm**2
            groups[group]['passage'][x, y] += 1
            do_the_thing(X, Y, x, y, raw, nrm, group)

    if bins:
        out.close()
    return groups
