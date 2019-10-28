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


def write_matrix(inbam, resolution, biases, outdir,
                 filter_exclude=(1, 2, 3, 4, 6, 7, 8, 9, 10),
                 region1=None, start1=None, end1=None, clean=True,
                 region2=None, start2=None, end2=None,
                 tmpdir='.', ncpus=8, verbose=True):

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

    if verbose:
        printime('  - Writing matrices')

    fnam = os.path.join(outdir, '{}_bam_{}.tsv'.format(os.path.split(outdir)[-1], nicer(resolution, sep='')))
    mkdir(outdir)
    out = open(fnam, 'w')

    # pull all sub-matrices and write full matrix
    for c, j, k, v in _iter_matrix_frags(chunks, tmpdir, rand_hash,
                                         verbose=verbose, clean=clean):
        if k < j or j in bads1 or k in bads2:  # we keep only half matrix
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


def sort_BAMtsv(outdir, tmp, resolution):
    tsv = os.path.join(outdir, "{}_bam_{}.tsv".format(
        os.path.split(outdir)[-1], nicer(resolution, sep='')))
    printime('Sorting BAM matrix: {}'.format(tsv))
    # sort file first and second column and write to same file
    _ = Popen("sort -k1n -k2n -S 10% {0} -T {1} -o {0}".format(tsv, tmp),
              shell=True).communicate()


def main():
    opts = get_options()
    inbam = opts.inbam
    resolution = opts.resolution
    outdir = opts.outdir
    biases_file = opts.biases_file

    write_matrix(inbam, resolution, biases_file, outdir,
                 ncpus=opts.ncpus, clean=opts.clean)

    rand_hash = "%016x" % getrandbits(64)
    tmpdir = os.path.join('.', '_tmp_%s' % (rand_hash))
    mkdir(tmpdir)

    #sort all files for only read once per pair of peaks to extract
    sort_BAMtsv(outdir, tmpdir, resolution)

    os.system('rm -rf {}'.format(tmpdir))

    printime('Done.')


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-bam', '--bam', dest='inbam', required=True, default=False,
                        help='Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--out', dest='outdir', required=True, default=False,
                        help='Outdir to store counts')
    parser.add_argument('-b', '--biases', dest='biases_file', required=True, default=False,
                        help='Pickle file with biases')
    parser.add_argument('--keep_tmp', dest='clean', default=True, action='store_false',
                        help='Keep temporary files for debugging')
    parser.add_argument('-C', dest='ncpus', default=cpu_count(),
                        type=int, help='Number of CPUs used to read BAM')
    opts = parser.parse_args()

    return opts


if __name__ == '__main__':
    exit(main())
