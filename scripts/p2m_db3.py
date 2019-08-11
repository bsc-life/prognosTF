#! /usr/bin/env python

"""
"""

from argparse    import ArgumentParser
from collections import OrderedDict, defaultdict
from datetime    import datetime
from time        import time
from pickle      import load, dump, HIGHEST_PROTOCOL
from subprocess  import Popen, PIPE

import os
import errno

from pysam       import AlignmentFile


def mkdir(dnam):
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(dnam):
            raise


def printime(msg):
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')


def binning_bed(peak_files, resolution, windows_span, max_dist,
                chrom_sizes, windows, in_feature, section_pos, badcols):
    '''
    Input BED file(s) of ChIP peaks and bin into desired resolution of Hi-C
    '''

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
        peaks2 = open(peak_files[1])
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

    bin_coordinate1 = set((c, p, f) for c, p, f in map(read_line1, peaks1)
                          if windows_span < p < max_chrom[c])

    peaks2.seek(0)  # needs to be here in case peaks1 and peak2 are the same
    bin_coordinate2 = set((c, p, f) for c, p, f in map(read_line2, peaks2)
                          if windows_span < p < max_chrom[c])

    printime('Total of different bin coordinates: {} and {}'.format(
        len(bin_coordinate1), len(bin_coordinate2)))

    # sort peaks
    bin_coordinate1 = sorted(bin_coordinate1)
    if same:
        bin_coordinate2 = bin_coordinate1[:]
    else:
        bin_coordinate2 = sorted(bin_coordinate2)

    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features

    wsp = (windows_span * 2) + 1
    mdr = max_dist / resolution

    printime('- Generating pairs of coordinates...')
    # put pairs in intervals
    for window in windows:
        if window == 'inter':
            test = lambda a, b: (a[0] != b[0]
                                 and a != b)
        elif window == 'intra':
            test = lambda a, b: (a[0] == b[0]
                                 and wsp <= abs(b[1] - a[1]) <= mdr
                                 and a != b)
        else:
            lower, upper = window
            test = lambda a, b: (a[0] == b[0]
                                 and wsp <= abs(b[1] - a[1]) <= mdr
                                 and a != b
                                 and lower < abs(a[1] - b[1]) <= upper)

        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if test(a, b))

        buf = []
        for a in bin_coordinate1:
            for b in bin_coordinate2:
                if test(a, b):
                    buf_beg = a
                    break
            else:
                continue
            break
        for (chr1, bs1, f1), (chr2, bs2, f2) in pairs:
            beg1, end1 = bs1 - windows_span, bs1 + windows_span
            beg2, end2 = bs2 - windows_span, bs2 + windows_span
            what = f1 + f2

            if beg1 > beg2 and chr1 == chr2:
                beg1, end1, beg2, end2 = beg2, end2, beg1, end1

            pos1 = section_pos[chr1][0]
            pos2 = section_pos[chr2][0]
            start_bin1 = pos1 + (beg1 / resolution)
            end_bin1   = pos1 + (end1 / resolution) + 1
            start_bin2 = pos2 + (beg2 / resolution)
            end_bin2   = pos2 + (end2 / resolution) + 1

            range1 = [(x, p1) for x, p1 in enumerate(range(start_bin1, end_bin1))
                      if p1 not in badcols]
            range2 = [(y, p2) for y, p2 in enumerate(range(start_bin2, end_bin2))
                      if p2 not in badcols]

            if not range1 or not range2:
                continue
            for x, p1 in range1:
                for y, p2 in range2:
                    buf.append(p1, p2, x, y, what)
            if end_bin1 - buf_beg > windows_span * 2:
                buf.sort()
                top = end_bin1 - windows_span
                p1 = 0
                while p1 < top:
                    p1, p2, x, y, what = buf.pop(0)
                    yield p1, p2, x, y, what
                buf_beg = p1
        while buf:
            yield buf.pop(0)


def readfiles(genomic_file, iter_pairs):
    def split_line1(line):
        a, b, c, d = line.split()
        return (int(a), int(b)), c, d

    # create empty meta-waffles
    printime(' - Reading genomic matrix and peaks')
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
    printime(' - Finished extracting')


def main():
    opts = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_files   = opts.peak_files
    outdir       = opts.outdir
    windows_span = opts.windows_span
    max_dist     = opts.max_dist
    windows      = opts.windows
    genomic_mat  = opts.genomic_mat
    ncpus        = opts.ncpus
    in_feature   = opts.first_is_feature
    biases       = opts.biases
    tmpdir       = opts.tmpdir
    if not tmpdir:
        tmpdir = os.path.join(outdir, 'tmp')

    mkdir(outdir)
    mkdir(tmpdir)

    windows = [[int(x) / resolution for x in win.split('-')]
               if win not in  ['inter', 'intra'] else win for win in windows]
    if 'intra' not in windows and 'inter' not in windows:
        for b, e in windows:
            if b >= e:
                raise Exception('ERROR: begining of windows should be smaller '
                                'than end')
    ## peaks file sorted per chromosome
    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references, [x / resolution + 1
                                                    for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    bins = {}
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        for n, i in enumerate(range(*section_pos[crm])):
            bins[i] = (crm, n)
        total += sections[crm]

    chrom_sizes = OrderedDict(zip(bamfile.references,
                                  [x for x in bamfile.lengths]))

    badcols = load(open(biases))['badcol']

    # get all lists... we assume that that will be ordered as wanted
    iter_pairs = binning_bed(peak_files, resolution, windows_span, max_dist,
                             chrom_sizes, windows, in_feature, section_pos, badcols)

    p = Popen('less', stdin=PIPE)
    proc = Popen(("sort -k9 -S 10% --parallel={0} "
                  "--temporary-directory={1} -o {2}").format(
                      ncpus, tmpdir, os.path.join(outdir, 'final_per_cell_sorted.tsv')),
                 shell=True)
    for X, Y, x, y, raw, nrm, group in readfiles(genomic_mat, iter_pairs):
        c1, b1 = bins[X]
        c2, b2 = bins[Y]
        try:
            proc.stdin.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                c1, b1, c2, b2, x, y, raw, nrm, group))
        except IOError as e:
            if e.errno == errno.EPIPE or e.errno == errno.EINVAL:
                # normal stop
                break
            else:
                raise
    p.stdin.close()
    p.wait()

    printime('Summing peaks')
    groups = defaultdict(int)
    for chunk in peaks:
        for l in open(peaks[chunk]['tmp_path'] + '_groups'):
            g, v = l.split()
            groups[g] += int(v)

    # initialize defaultdict and prev variable
    sum_raw = defaultdict(int)
    sqr_raw = defaultdict(int)
    sum_nrm = defaultdict(float)
    sqr_nrm = defaultdict(float)
    passage = defaultdict(int)
    fhandler = open(os.path.join(outdir, 'final_per_cell_sorted.tsv'))
    line = next(fhandler)

    c1, b1, c2, b2, x, y, raw, nrm, group = line.split()
    prev = group
    fhandler.seek(0)
    out = open(os.path.join(outdir, 'final_sum.pickle'), 'wb')
    for line in fhandler:
        c1, b1, c2, b2, x, y, raw, nrm, group = line.split()
        if group != prev:
            dump([prev, (sum_raw, sqr_raw, sum_nrm, sqr_nrm, passage, groups[prev])],
                 out, protocol=HIGHEST_PROTOCOL)
            sum_raw = defaultdict(int)
            sqr_raw = defaultdict(int)
            sum_nrm = defaultdict(float)
            sqr_nrm = defaultdict(float)
            passage = defaultdict(int)
            prev = group
        raw = int(raw)
        nrm = float(nrm)
        y = int(y)
        x = int(x)
        sum_raw[x, y] += raw
        sqr_raw[x, y] += raw**2
        sum_nrm[x, y] += nrm
        sqr_nrm[x, y] += nrm**2
        passage[x, y] += 1
    dump([group, (sum_raw, sqr_raw, sum_nrm, sqr_nrm, passage, groups[group])],
         out, protocol=HIGHEST_PROTOCOL)
    out.close()

    # clean
    os.system('rm -rf ' + tmpdir)

    printime('all done!')


def get_options():
    parser = ArgumentParser(usage="-i Peaks -r INT [options]")

    windows = ('255000-1000000',
               '1000000-2500000',
               '2500000-5000000',
               '5000000-10000000',
               '10000000-15000000',
               '15000000-20000000',
               '20000000-30000000')

    parser.add_argument('-i', '--peaks', dest='peak_files', required=True,
                        nargs="+", metavar='PATH',
                        help='''one or two pairwise peaks files to
                        compute average submatrix (norm and raw). These files
                        should contain at least two columns, chromosome and
                        position, and may contain an extra column of feature. If
                        present, the resuilt will be retrurned according to the
                        possible combination of this feature''')
    parser.add_argument('-b', '--bam', dest='inbam', required=True,
                        metavar='PATH', help='Input HiC-BAM file')
    parser.add_argument('-b', '--biases', dest='biases', default=True, help='Biases',
                        required=True)
    parser.add_argument('-r', '--resolution', dest='resolution', required=True,
                        metavar='INT', default=False, type=int,
                        help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--outdir', dest='outdir', default='',
                        metavar='PATH', help='output directory')
    parser.add_argument('-s', dest='windows_span', required=True, type=int,
                        metavar='INT',
                        help='''Windows span around center of the peak (in bins; the
                        total windows size is 2 times windows-span + 1)''')
    parser.add_argument('-m', '--max_dist', dest='max_dist', metavar='INT',
                        default=float('inf'), type=int,
                        help='''[%(default)s] Max dist between center peaks''')
    parser.add_argument('-w', '--windows', dest='windows', required=False,
                        default=windows, metavar='INT-INT', type=str, nargs="+",
                        help='''[%(default)s] If only interested in some
                        intervals to check: "-w 1000000-2000000 2000000-5000000"
                        correspond to 2 window intervals, one from 1Mb to 2Mb
                        and one from 2Mb to 5Mb. Use "-w inter" for
                        inter-chromosomal regions, or "-w intra" for
                        intra-chromosomal (whithout distance restriction)''')
    parser.add_argument('--first_is_feature', dest='first_is_feature', default=False,
                        action='store_true', help='''When 2 BED files are input,
                        the peaks in the first BED should be also considered as
                        feature. This is to create average sub-matrices for
                        each peak in the first BED file.''')
    parser.add_argument('-C', '--cpus', dest='ncpus', type=int, default=8,
                        help='''[%(default)s] number of cpus to be used for parsing
                        the HiC-BAM file''')
    parser.add_argument('-t', '--tmp', dest='tmpdir', default=None,
                        help='Tmpdir to store coordinates files')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
