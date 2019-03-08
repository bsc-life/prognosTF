__author__ = 'sgalan'

import os

from glob                            import glob
from subprocess                      import Popen
from multiprocessing                 import cpu_count
from argparse                        import ArgumentParser
from collections                     import OrderedDict
from cPickle                         import load
from math                            import isnan
from StringIO                        import StringIO
from datetime                        import datetime

from pytadbit.parsers.hic_bam_parser import get_biases_region, _iter_matrix_frags
from pytadbit.parsers.hic_bam_parser import printime, get_matrix
from pytadbit.parsers.hic_bam_parser import read_bam, filters_to_bin
from pytadbit.utils.extraviews       import nicer
from pytadbit.utils.file_handling    import mkdir

import pysam


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

    bamfile = pysam.AlignmentFile(inbam, 'rb')
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

    fnam = os.path.join(outdir, '{}_bam_{}kb.tsv'.format(region1, resolution / 1000))

    out = open(os.path.join(fnam), 'w')

    # pull all sub-matrices and write full matrix
    for c, j, k, v in _iter_matrix_frags(chunks, tmpdir, rand_hash,
                                         verbose=verbose, clean=clean):
        if k < j: # we are only going to keep half of the matrix
            continue
        if j not in bads1 and k not in bads2:
            n = v / bias1[j] / bias2[k] / decay[c][abs(j-k)]
            pos1 = j + section_pos[region1][0]
            pos2 = k + section_pos[region1][0]
            out.write('{}\t{}\t{}\t{}\n'.format(pos1, pos2, v, n))

    out.close()

    # this is the last thing we do in case something goes wrong
    os.system('rm -rf %s' % (os.path.join(tmpdir, '_tmp_%s' % (rand_hash))))

    if  verbose:
        printime('\nDone.')


def sort_BAMtsv(outdir, resolution):
    all_tsv = glob(outdir + "*_bam_%ikb.tsv"%(resolution / 1000))
    for tsv in all_tsv:
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Sorting BAM matrix: ', tsv)
        # sort file first and second column and write to same file
        fh = Popen("sort -k1n -k2n -S 20% "+ tsv+ " -o "+ tsv, shell=True)
        fh.communicate()


def main():
    opts = get_options()
    inbam = opts.inbam
    resolution = opts.resolution
    outdir = opts.outdir
    biases_file = opts.biases_file

    bamfile = pysam.AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references,
                               [x / resolution + 1 for x in bamfile.lengths]))

    # remove nan from biases ONE-D
    for chromosome in sections.keys():
        print(datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
              ' - Splitting peak pairs per chromosome...')
        write_matrix(inbam, resolution, biases_file, outdir, region1=chromosome,
                     ncpus=opts.ncpus)

    #sort all files for only read once per pair of peaks to extract
    sort_BAMtsv(outdir, resolution)


def get_options():
    parser = ArgumentParser(usage="-bam Bam -r INT [options]")
    parser.add_argument('-bam', '--bam', dest='inbam', required=True, default=False,
                        help='Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--out', dest='outdir', required=True, default=False,
                        help='Outdir to store counts')
    parser.add_argument('-b', '--biases', dest='biases_file', required=True, default=False,
                        help='Pickle file with biases')
    parser.add_argument('-C', dest='ncpus', required=True, default=cpu_count(),
                        type=int, help='Number of CPUs used to read BAM')
    opts = parser.parse_args()

    return opts


if __name__ == '__main__':
    exit(main())
