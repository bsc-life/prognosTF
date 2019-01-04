from argparse    import ArgumentParser
from collections import OrderedDict
from datetime    import datetime
from itertools   import combinations
from os          import path

from pysam       import AlignmentFile


def binning_bed(peak_file, resolution, windows_span, max_dist, outdir,
                name, chrom_sizes, windows, **kwargs):
    '''Input bed file of ChIP peaks and bin into desired resolution of Hi-C'''

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

    peaks = open(peak_file, "r")

    # findout if bed file contain features, or only coordinates
    line = peaks.next()
    try:
        read_line_feature(line)
        read_line = read_line_feature
    except ValueError:
        read_line = read_line_no_feature
    peaks.seek(0)

    windows_span /= resolution

    bin_coordinate = set((c, p, f) for c, p, f in map(read_line, peaks)
                          if p > windows_span)  # take into account to add windows span both sides

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Total of different bin coordinates: ', len(bin_coordinate)

    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features

    wsp = (windows_span * 2) + 1
    mdr = max_dist / resolution

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), '- Writing list of coordinates...'

    intervals = {}
    pairs = ((a, b) for a, b in combinations(bin_coordinate, 2)
             if a[0] == b[0] and wsp <= abs(b[1] - a[1]) <= mdr)
    num_win = sorted([w for w in windows if w not in ['inter', 'intra']])
    # we can go faster if windows are not ovelapping
    for i, (b, e) in enumerate(num_win[:-1]):
        if max((b, e)) > min(num_win[i + 1]):
            overlapping = True
            break
    else:
        overlapping = False
    # put pairs in intervals
    if overlapping:
        for (chrom1, bs1, f1), (chrom2,  bs2, f2) in pairs:
            distance = abs(bs1 - bs2)
            for lower, upper in (w for w in windows if w not in ['inter', 'intra']):
                if lower < distance <= upper:
                    intervals.setdefault((lower, upper), {})
                    intervals[(lower, upper)].setdefault(
                    (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
    else:
        for (chrom1, bs1, f1), (chrom2,  bs2, f2) in pairs:
            distance = abs(bs1 - bs2)
            for lower, upper in (w for w in windows if w not in ['inter', 'intra']):
                if lower < distance <= upper:
                    intervals.setdefault((lower, upper), {})
                    intervals[(lower, upper)].setdefault(
                    (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
                    break
    if 'inter' in windows:
        pairs = ((a, b) for a, b in combinations(bin_coordinate, 2)
                if a[0] != b[0])
        intervals['inter'] = {}
        for (chrom1, bs1, f1), (chrom2,  bs2, f2) in pairs:
            intervals['inter'].setdefault(
                (f1, f2), []).append((chrom1, bs1, chrom2, bs2))
    if 'intra' in windows:
        pairs = ((a, b) for a, b in combinations(bin_coordinate, 2)
                if a[0] == b[0] and wsp <= abs(b[1] - a[1]) <= mdr)
        intervals['intra'] = {}
        for (chrom1, bs1, f1), (chrom2,  bs2, f2) in pairs:
            intervals['intra'].setdefault(
                (f1, f2), []).append((chrom1, bs1, chrom2, bs2))

    # define categories and write pairs
    for window in intervals:
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'Writing interval: ', window
        for f1, f2 in intervals[window]:
            extra = ('both' if f1 == f2 else f1) + '_' + f2
            if window != 'inter':
                window_name = '%d_%d' % (beg * resolution, end * resolution)
            else:
                window_name = window
            w = open(path.join(outdir, '%s_%s_%s.tsv' % (
                name, window_name, extra)), 'w')
            for c1, s1, c2, s2 in intervals[window][(f1, f2)]:
                start1, end1 = s1 - windows_span, s1 + windows_span
                start2, end2 = s2 - windows_span, s2 + windows_span
                # check chromosome length
                new_start1, new_end1 = start1 * resolution, end1 * resolution
                new_start2, new_end2 = start2 * resolution, end2 * resolution
                if new_end1 <= chrom_sizes[c1] and new_end2 <= chrom_sizes[c1]:
                    w.write('%s\t%d\t%d\t%s\t%d\t%d\n' % (
                       c1, new_start1, new_end1, c2, new_start2, new_end2))
            w.close()


def main():
    opts = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_file    = opts.peak_file
    outdir       = opts.outdir
    name         = opts.name
    windows_span = opts.windows_span
    max_dist     = opts.max_dist
    windows      = opts.windows

    windows = [[int(x) / resolution for x in win.split('-')]
               if win not in  ['inter', 'intra'] else win for win in windows]
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

    binning_bed(peak_file, resolution, windows_span, max_dist, outdir, name,
                chrom_sizes, windows)

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Sublists written!'


def get_options():
    parser = ArgumentParser(usage="-i Peaks -r INT [options]")

    windows = ('255000-1000000',
               '1000000-2500000',
               '2500000-5000000',
               '5000000-10000000',
               '10000000-15000000',
               '15000000-20000000',
               '20000000-30000000')

    parser.add_argument('-i','--peak', dest='peak_file', required=True,
                        metavar='PATH', help='''Pairwise peaks to compute average submatrix
                        (norm and raw)''')
    parser.add_argument('-b','--bam', dest='inbam', required=True,
                        metavar='PATH', help= 'Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True,
                        metavar='INT', default=False, type=int,
                        help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--outdir', dest='outdir', default='',
                        metavar='PATH', help='output directory')
    parser.add_argument('-n', '--name', dest='name', default='pairs',
                        help='[%(default)s] Output name prefix')
    parser.add_argument('-s', dest='windows_span',required=True, type=int,
                        metavar='INT',
                        help='''Windows span around center of the peak (Total
                        windows size is 2 times windows-span + 1)''')
    parser.add_argument('-m','--max_dist', dest='max_dist', metavar='INT',
                        default=float('inf'), type=int,
                        help='''[%(default)s] Max dist between center peaks''')
    parser.add_argument('-w','--windows', dest='windows', required=False,
                        default=windows, metavar='INT-INT', type=str,nargs="+",
                        help='''[%(default)s] If only interested in some
                        intervals to check: "-w 1000000-2000000 2000000-5000000"
                        correspond to 2 window intervals, one from 1Mb to 2Mb
                        and one from 2Mb to 5Mb. Use "-w inter" for
                        inter-chromosomal regions, or "-w intra" for
                        intra-chromosomal (whithout distance restriction)''')

    opts = parser.parse_args()
    return opts

if __name__=='__main__':
    exit(main())
