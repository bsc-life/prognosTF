#! /usr/bin/env python

"""
"""

from argparse    import ArgumentParser
from collections import OrderedDict, defaultdict
from datetime    import datetime
from time        import time
from pickle      import load, dump, HIGHEST_PROTOCOL, _Unpickler
from subprocess  import Popen, PIPE

import os
import errno

from pysam       import AlignmentFile


def mkdir(dnam):
    dnam = os.path.abspath(dnam)
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(dnam):
            raise


def printime(msg, silent=True):
    if silent:
        return
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')


def chromosome_from_bam(inbam, resolution, get_bins=False):
    ## peaks file sorted per chromosome
    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references, [x // resolution + 1
                                                    for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    bins = {}
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        if get_bins:
            for n, i in enumerate(range(*section_pos[crm])):
                bins[i] = (crm, n)
        total += sections[crm]

    chrom_sizes = OrderedDict(zip(bamfile.references,
                                  [x for x in bamfile.lengths]))

    return section_pos, chrom_sizes, bins


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
        same = False
    except IndexError:
        same = True
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
    if same:
        bin_coordinate2 = bin_coordinate1
    else:
        bin_coordinate2 = sorted(bin_coordinate2)

    return bin_coordinate1, bin_coordinate2, npeaks1, npeaks2, same


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


def submatrix_coordinates(final_pairs, badcols, wsp):
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
    def split_line1(line):
        a, b, c, d = line.split()
        return (int(a), int(b)), c, d

    # create empty meta-waffles
    fh1 = open(genomic_file)
    pos1, raw, nrm = split_line1(next(fh1))
    pos2, x, y, e = next(iter_pairs)

    try:
        while True:
            if pos2 > pos1:
                pos1, raw, nrm = split_line1(next(fh1))
            elif pos1 == pos2:
                raw = int(raw)
                nrm = float(nrm)
                yield pos1[0], pos1[1], x, y, raw, nrm, e
                pos2, x, y, e = next(iter_pairs)
                if pos1 != pos2:  # some cells in the peak file are repeated
                    pos1, raw, nrm = split_line1(next(fh1))
            else:
                pos2, x, y, e = next(iter_pairs)
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
    for X, Y, x, y, raw, nrm, group in readfiles(genomic_mat, iter_pairs):
        try:
            groups[group]['counter'] += 1
        except KeyError:
            groups[group] = {
                'counter' : 1,
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


def main():
    opts = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_files   = opts.peak_files
    outfile      = opts.outfile
    windows_span = opts.windows_span
    max_dist     = opts.max_dist
    window       = opts.window
    genomic_mat  = opts.genomic_mat
    ncpus        = opts.ncpus
    in_feature   = opts.first_is_feature
    biases       = opts.biases
    submatrices  = opts.submatrices
    silent       = opts.silent
    badcols      = _Unpickler(open(biases, "rb"), encoding='latin1').load()['badcol']

    window_name = window
    if window not in  ['inter', 'intra', 'all']:
        window = [int(x) / resolution for x in window.split('-')]
        if window[0] >= window[1]:
            raise Exception('ERROR: beginning of window should be smaller '
                            'than end')

    mkdir(os.path.split(outfile)[0])

    # get chromosome coordinates and conversor genomic coordinate to bins
    section_pos, chrom_sizes, bins = chromosome_from_bam(
        inbam, resolution, get_bins=True if submatrices else False)

    # define pairs of peaks
    peak_coord1, peak_coord2, npeaks1, npeaks2, same = parse_peaks(
        peak_files, resolution, in_feature, chrom_sizes, windows_span)
    printime('Total of different/usable peak bin coordinates in {}:'.format(
        peak_files[0]), silent)
    printime(('   - {} (out of {})').format(
        len(peak_coord1), npeaks1), silent)
    if not same:
        printime('Total of different/usable peak bin coordinates in {}:'.format(
            peak_files[1]), silent)
        printime(('   - {} (out of {})').format(
            len(peak_coord2), npeaks2), silent)

    printime('- Generating pairs of coordinates...', silent)
    pair_peaks = generate_pairs(peak_coord1, peak_coord2, resolution,
                                windows_span, max_dist, window, section_pos)

    iter_pairs = submatrix_coordinates(pair_peaks, badcols, (windows_span * 2) + 1)

    # retrieve interactions at peak pairs using genomic matrix
    # sum them by feature and store them in dictionary
    printime(' - Reading genomic matrix and peaks', silent)
    groups = interactions_at_intersection(genomic_mat, iter_pairs, submatrices, bins)
    printime(' - Finished extracting', silent)

    out = open(outfile, 'wb')
    dump(groups, out, protocol=HIGHEST_PROTOCOL)
    out.close()

    printime('all done!', silent)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('--peaks', dest='peak_files', required=True,
                        nargs="+", metavar='PATH',
                        help='''one or two pairwise peaks files to
                        compute average submatrix (norm and raw). These files
                        should contain at least two columns, chromosome and
                        position, and may contain an extra column of feature. If
                        present, the result will be returned according to the
                        possible combination of this feature''')
    parser.add_argument('--bam', dest='inbam', required=True,
                        metavar='PATH', help='Input HiC-BAM file')
    parser.add_argument('--biases', dest='biases', default=True, help='Biases',
                        required=True)
    parser.add_argument('--genomic_matrix', dest='genomic_mat', default=True,
                        metavar='GENOMIC_MATRIX', required=True,
                        help='''Path to genomic matrix in 3 columns format
                        (should be sorted with `sort -k1,2n`)''')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True,
                        metavar='INT', default=False, type=int,
                        help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        metavar='PATH', help='path to output file (pickle format)')
    parser.add_argument('--all_submatrices', dest='submatrices', default='',
                        metavar='PATH', help='''if PATH is provided here, stores
                        all the individual submatrices generated''')
    parser.add_argument('-s', dest='windows_span', required=True, type=int,
                        metavar='INT',
                        help='''Windows span around center of the peak (in bins; the
                        total windows size is 2 times windows-span + 1)''')
    parser.add_argument('-m', '--max_dist', dest='max_dist', metavar='INT',
                        default=float('inf'), type=int,
                        help='''[%(default)s] Max dist between center peaks''')
    parser.add_argument('-w', '--window', dest='window', required=False,
                        default='intra', metavar='INT-INT', type=str,
                        help='''[%(default)s] If only interested in some
                        intervals to check: "-w 1000000-2000000"
                        correspond to the window interval, from 1Mb to 2Mb.
                        Use "-w inter" for inter-chromosomal regions, "-w intra" for
                        intra-chromosomal, "-w all" for all combinations
                        (without distance restriction)''')
    parser.add_argument('--first_is_feature', dest='first_is_feature', default=False,
                        action='store_true', help='''When 2 BED files are input,
                        the peaks in the first BED should be also considered as
                        feature. This is to create average sub-matrices for
                        each peak in the first BED file.''')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')
    parser.add_argument('-C', '--cpus', dest='ncpus', type=int, default=8,
                        help='''[%(default)s] number of cpus to be used for parsing
                        the HiC-BAM file''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
