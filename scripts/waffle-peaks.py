#! /usr/bin/env python

"""
"""

from os.path     import split as os_split

from argparse    import ArgumentParser
try:  # python 3
    from pickle        import dump, HIGHEST_PROTOCOL, _Unpickler as Unpickler
except ImportError:  # python 2
    from pickle        import dump, HIGHEST_PROTOCOL, Unpickler

try:
    from meta_waffle       import chromosome_from_bam, parse_peaks, generate_pairs
    from meta_waffle       import submatrix_coordinates, interactions_at_intersection
    from meta_waffle.utils import printime, mkdir
except ImportError:  # meta-waffle is not installed.. but it's still ok!!!
    from os.path  import join as os_join
    from sys.path import insert

    insert(0, os_join(os_split(os_split(__file__)[0])[0], 'meta_waffle'))


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
    in_feature   = opts.first_is_feature
    biases       = opts.biases
    submatrices  = opts.submatrices
    silent       = opts.silent
    badcols      = Unpickler(open(biases, "rb"), encoding='latin1').load()['badcol']

    if window not in  ['inter', 'intra', 'all']:
        window = [int(x) / resolution for x in window.split('-')]
        if window[0] >= window[1]:
            raise Exception('ERROR: beginning of window should be smaller '
                            'than end')

    mkdir(os_split(outfile)[0])

    # get chromosome coordinates and conversor genomic coordinate to bins
    section_pos, chrom_sizes, bins = chromosome_from_bam(
        inbam, resolution, get_bins=submatrices!='')

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

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
