__author__ = 'sgalan'

import os
import matplotlib
matplotlib.use('Agg')

from subprocess                      import Popen
from multiprocessing                 import cpu_count
from argparse                        import ArgumentParser
from collections                     import OrderedDict
from datetime                        import datetime
from random                          import getrandbits

from pytadbit.parsers.hic_bam_parser import get_biases_region, _iter_matrix_frags
from pytadbit.parsers.hic_bam_parser import read_bam, filters_to_bin, printime
from pytadbit.utils.file_handling    import mkdir
from pytadbit.utils.extraviews       import nicer

from pysam                           import AlignmentFile


def write_matrix(inbam, resolution, biases, outfile,
                 filter_exclude=(1, 2, 3, 4, 6, 7, 8, 9, 10),
                 region1=None, start1=None, end1=None, clean=True,
                 region2=None, start2=None, end2=None,
                 tmpdir='.', ncpus=8, verbose=True, window=None):

    if not isinstance(filter_exclude, int):
        filter_exclude = filters_to_bin(filter_exclude)

    _, rand_hash, bin_coords, chunks = read_bam(
        inbam, filter_exclude, resolution, ncpus=ncpus,
        region1=region1, start1=start1, end1=end1,
        region2=region2, start2=start2, end2=end2,
        tmpdir=tmpdir, verbose=verbose)

    bamfile = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references,
                               [x / resolution + 1 for x in bamfile.lengths]))

    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    if biases:
        bias1, bias2, decay, bads1, bads2 = get_biases_region(biases, bin_coords)
    else:
        bads1 = bads2 = {}

    if bads1 is bads2:
        badcols = bads1
    else:  # should never happen
        badcols = bads1
        badcols.update(bads2)

    if verbose:
        printime('  - Writing matrices')

    mkdir(os.path.split(os.path.abspath(outfile))[0])
    # write the rest of the file to be sorted
    out = open(outfile, 'w')
    nheader = 0
    for i, c in enumerate(bamfile.references):
        out.write('# CHROM\t{}\t{}\n'.format(c, bamfile.lengths[i]))
        nheader += 1
    out.write('# RESOLUTION\t{}\n'.format(resolution))
    nheader += 1
    out.write('# BADCOLS\t{}\n'.format(','.join(map(str, badcols.keys()))))
    nheader += 1

    if window == 'all':
        outside = lambda c_, j_, k_: False
    elif window == 'intra':
        outside = lambda c_, j_, k_: c_ == ''
    elif window == 'inter':
        outside = lambda c_, j_, k_: c_ != ''
    else:
        min_, max_ = window
        outside = lambda c_, j_, k_: (k_ - j_) < min_ or (k_ - j_) > max_

    # pull all sub-matrices and write full matrix
    for c, j, k, v in _iter_matrix_frags(chunks, tmpdir, rand_hash,
                                         verbose=verbose, clean=clean):
        if k < j or j in badcols or k in badcols:  # we keep only half matrix
            continue
        if outside(c, j, k):
            continue
        try:
            n = v / bias1[j] / bias2[k] / decay[c][k - j]  # normalize
        except KeyError:
            n = v / bias1[j] / bias2[k]  # normalize
        out.write('{}\t{}\t{}\t{}\n'.format(j, k, v, n))
    out.close()

    # this is the last thing we do in case something goes wrong
    if clean:
        os.system('rm -rf %s' % (os.path.join(tmpdir, '_tmp_%s' % (rand_hash))))
    return nheader


def sort_BAMtsv(nheader, outfile, tmp):
    tsv = outfile
    printime('Sorting BAM matrix: {}'.format(tsv))
    # sort file first and second column and write to same file
    print(("(head -n {0} {1} && tail -n +{0} {1} | "
               "sort -k1n -k2n -S 10% -T {2}) > {1}").format(
                   nheader, tsv, tmp))
    _ = Popen(("(head -n {0} {2} && tail -n +{1} {2} | "
               "sort -k1n -k2n -S 10% -T {3}) > {2}_").format(
                   nheader, nheader + 1, tsv, tmp), shell=True).communicate()
    os.system("mv {0}_ {0}".format(tsv))


def main():
    opts        = get_options()
    inbam       = opts.inbam
    resolution  = opts.resolution
    outfile     = opts.outfile
    biases_file = opts.biases_file
    window      = opts.window

    if window not in  ['inter', 'intra', 'all']:
        window = [int(x) / resolution for x in window.split('-')]
        if window[0] >= window[1]:
            raise Exception('ERROR: beginning of window should be smaller '
                            'than end')

    nheader = write_matrix(inbam, resolution, biases_file, outfile,
                           ncpus=opts.ncpus, clean=opts.clean, window=window)

    rand_hash = "%016x" % getrandbits(64)
    tmpdir = os.path.join('.', '_tmp_%s' % (rand_hash))
    mkdir(tmpdir)

    #sort all files for only read once per pair of peaks to extract
    sort_BAMtsv(nheader, outfile, tmpdir)

    os.system('rm -rf {}'.format(tmpdir))

    printime('Done.')


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-bam', '--bam', dest='inbam', required=True, default=False,
                        help='Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--out', dest='outfile', required=True, default=False,
                        help='Output file to store counts')
    parser.add_argument('-b', '--biases', dest='biases_file', required=True, default=False,
                        help='Pickle file with biases')
    parser.add_argument('--keep_tmp', dest='clean', default=True, action='store_false',
                        help='Keep temporary files for debugging')
    parser.add_argument('-C', dest='ncpus', default=cpu_count(),
                        type=int, help='Number of CPUs used to read BAM')
    parser.add_argument('-w', '--window', dest='window', required=False,
                        default='all', metavar='INT-INT', type=str,
                        help='''[%(default)s] If only interested in some
                        intervals to check: "-w 1000000-2000000"
                        correspond to the window interval, from 1Mb to 2Mb.
                        Use "-w inter" for inter-chromosomal regions, "-w intra" for
                        intra-chromosomal, "-w all" for all combinations
                        (without distance restriction)''')
    opts = parser.parse_args()

    return opts


if __name__ == '__main__':
    exit(main())
