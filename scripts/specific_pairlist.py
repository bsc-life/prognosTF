#! /usr/bin/env python

"""
"""

from argparse    import ArgumentParser
from collections import OrderedDict
from datetime    import datetime
from time        import time
from subprocess  import Popen

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


def binning_bed(peak_files, resolution, windows_span, max_dist, outdir,
                name, chrom_sizes, windows, ncpus, tmpdir, **kwargs):
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

    peaks1 = open(peak_files[0], "r")
    try:
        peaks2 = open(peak_files[1])
    except IndexError:
        peaks2 = peaks1

    # findout if bed file contain features, or only coordinates
    line = next(peaks1)
    try:
        read_line_feature(line)
        read_line1 = read_line_feature
    except ValueError:
        read_line1 = read_line_no_feature
    peaks1.seek(0)

    line = next(peaks2)
    try:
        read_line_feature(line)
        read_line2 = read_line_feature
    except ValueError:
        read_line2 = read_line_no_feature

    bin_coordinate1 = set((c, p, f) for c, p, f in map(read_line1, peaks1)
                          if p > windows_span)  # TODO: also remove position too close to the end

    peaks2.seek(0)  # needs to be here in case peaks1 and peak2 are the same
    bin_coordinate2 = set((c, p, f) for c, p, f in map(read_line2, peaks2)
                          if p > windows_span)  # TODO:  also remove position too close to the end

    printime('Total of different bin coordinates: {} and {}'.format(
        len(bin_coordinate1), len(bin_coordinate2)))

    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features

    wsp = (windows_span * 2) + 1
    mdr = max_dist / resolution

    printime('- Generating pairs of coordinates...')

    intervals = {}

    num_win = sorted([w for w in windows if w not in ['inter', 'intra']])
    # we can go faster if windows are not ovelapping
    for i, (beg, end) in enumerate(num_win[:-1]):
        if max((beg, end)) > min(num_win[i + 1]):
            overlapping = True
            break
    else:
        overlapping = False

    # put pairs in intervals
    if 'inter' in windows:
        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if a[0] != b[0] and a != b)
        intervals['inter'] = {}
        for (chrom1, bs1, f1), (chrom2, bs2, f2) in pairs:
            intervals['inter'].setdefault(
                (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
    elif 'intra' in windows:
        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if a[0] == b[0] and wsp <= abs(b[1] - a[1]) <= mdr and a != b)
        intervals['intra'] = {}
        for (chrom1, bs1, f1), (chrom2, bs2, f2) in pairs:
            intervals['intra'].setdefault(
                (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
    else:
        pairs = ((a, b) for a in bin_coordinate1 for b in bin_coordinate2
                 if a[0] == b[0] and wsp <= abs(b[1] - a[1]) <= mdr and a != b)
        if overlapping:
            for (chrom1, bs1, f1), (chrom2, bs2, f2) in pairs:
                distance = abs(bs1 - bs2)
                for lower, upper in (w for w in windows):
                    if lower < distance <= upper:
                        intervals.setdefault((lower, upper), {})
                        intervals[(lower, upper)].setdefault(
                            (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
        else:
            for (chrom1, bs1, f1), (chrom2, bs2, f2) in pairs:
                distance = abs(bs1 - bs2)
                for lower, upper in (w for w in windows):
                    if lower < distance <= upper:
                        intervals.setdefault((lower, upper), {})
                        intervals[(lower, upper)].setdefault(
                            (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
                        break

    # define categories and write pairs
    for window in intervals:
        printime('Writing interval: {}'.format(window))
        for f1, f2 in intervals[window]:
            extra = ('both' if f1 == f2 else f1) + '_' + f2
            if window not in ['intra', 'inter']:
                window_name = '{}_{}'.format(window[0] * resolution,
                                             window[1] * resolution)
            else:
                window_name = window
            fname = os.path.join(tmpdir, '%s_%s_%s.tsv' % (
                name, window_name, extra))
            out = open(fname, 'w')
            for c1, s1, c2, s2 in intervals[window][(f1, f2)]:
                start1, end1 = s1 - windows_span, s1 + windows_span
                start2, end2 = s2 - windows_span, s2 + windows_span
                # check chromosome length
                new_start1, new_end1 = start1 * resolution, end1 * resolution
                new_start2, new_end2 = start2 * resolution, end2 * resolution
                if new_end1 <= chrom_sizes[c1] and new_end2 <= chrom_sizes[c1]:
                    out.write('%s\t%d\t%d\t%s\t%d\t%d\n' % (
                        c1, new_start1, new_end1, c2, new_start2, new_end2))
            out.close()

            # sort it
            Popen(("sort -k4,6 -S 10% --parallel={0} {1} "
                   "--temporary-directory={2} > {1}_sorted").format(
                       ncpus, fname, tmpdir), shell=True).communicate()

            out = open(os.path.join(outdir, '%s_%s_%s.tsv' % (
                name, window_name, extra)), 'w')
            out.write('\n'.join('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{3}:{4}'.format(
                *line.split()) for line in open(fname + '_sorted')) + '\n')
            out.close()


def main():
    opts = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_files   = opts.peak_files
    outdir       = opts.outdir
    name         = opts.name
    windows_span = opts.windows_span
    max_dist     = opts.max_dist
    windows      = opts.windows
    ncpus        = opts.ncpus
    tmpdir       = opts.tmpdir
    if not tmpdir:
        tmpdir = os.path.join(outdir, 'tmp')

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
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    chrom_sizes = OrderedDict(zip(bamfile.references,
                                  [x for x in bamfile.lengths]))

    binning_bed(peak_files, resolution, windows_span, max_dist, outdir, name,
                chrom_sizes, windows, ncpus, tmpdir)

    # clean
    os.system('rm -rf ' + tmpdir)

    printime('Sublists written!')


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
    parser.add_argument('-r', '--resolution', dest='resolution', required=True,
                        metavar='INT', default=False, type=int,
                        help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--outdir', dest='outdir', default='',
                        metavar='PATH', help='output directory')
    parser.add_argument('-n', '--name', dest='name', default='pairs',
                        help='[%(default)s] Output name prefix')
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
    parser.add_argument('-C', '--cpus', dest='ncpus', type=int, default=8,
                        help='''[%(default)s] number of cpus to be used for parsing
                        the HiC-BAM file''')
    parser.add_argument('-t', '--tmp', dest='tmpdir', default=None,
                        help='Tmpdir to store coordinates files')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
