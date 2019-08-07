
from collections     import OrderedDict, defaultdict
from argparse        import ArgumentParser
from datetime        import datetime
from time            import time, sleep
from subprocess      import Popen
from glob            import glob
from multiprocessing import Pool
from cPickle         import dump, load

import os
import errno

import numpy as np
from pysam import AlignmentFile


def printime(msg):
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')


def mkdir(dnam):
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(dnam):
            raise


def extract_coordinates(peak_fname, tmp_fname, resolution, section_pos,
                        badcols):
    '''
    Chunk file into multiple, and write them in parallel per file write coord
    of 10,000 peak pairs Write a dictionary depending of pairs of peak per
    target
    '''
    out = open(tmp_fname, 'w')
    for line in open(peak_fname):
        chr1, beg1, end1, chr2, beg2, end2 = line.split()
        beg1, end1, beg2, end2 = int(beg1), int(end1), int(beg2), int(end2)
        if beg1 > beg2:
            beg1, end1, beg2, end2 = beg2, end2, beg1, end1
        pos1 = section_pos[chr1][0]
        pos2 = section_pos[chr2][0]
        start_bin1 = pos1 + (beg1 / resolution)
        end_bin1   = pos1 + (end1 / resolution) + 1
        start_bin2 = pos2 + (beg2 / resolution)
        end_bin2   = pos2 + (end2 / resolution) + 1

        range1 = [(x, p1) for x, p1 in enumerate(xrange(start_bin1, end_bin1))
                  if p1 not in badcols]
        range2 = [(y, p2) for y, p2 in enumerate(xrange(start_bin2, end_bin2))
                  if p2 not in badcols]
        if range1 and range2:
            out.write('\n'.join('{}\t{}\t{}\t{}'.format(p1, p2, x, y)
                                for x, p1 in range1 for y, p2 in range2) + '\n')
    out.close()


def readfiles(file1, file2):
    def split_line1(l):
        a, b, c, d = l.split()
        return (int(a), int(b)), c, d
    def split_line2(l):
        a, b, c, d = map(int, l.split())
        return (a, b), c, d

    # create empty meta-waffles
    sum_raw = defaultdict(int)
    sqr_raw = defaultdict(int)
    sum_nrm = defaultdict(float)
    sqr_nrm = defaultdict(float)
    passage = defaultdict(int)
    printime(' - Reading BAM and peaks from %s...' % file2)
    fh1 = open(file1)
    try:
        fh2 = open(file2)
    except IOError:
        return
    pos1, raw, nrm = split_line1(next(fh1))
    pos2, x, y = split_line2(next(fh2))
    try:
        while True:
            if pos2 > pos1:
                pos1, raw, nrm = split_line1(fh1.next())
            elif pos1 == pos2:
                raw = int(raw)
                nrm = float(nrm)
                sum_raw[x, y] += raw
                sum_nrm[x, y] += nrm
                sqr_raw[x, y] += raw**2
                sqr_nrm[x, y] += nrm**2
                passage[x, y] += 1
                pos2_ = pos2
                pos2, x, y = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated
                    pos1, raw, nrm = split_line1(fh1.next())
            else:
                passage[x, y] += 1
                pos2, x, y = split_line2(fh2.next())
    except StopIteration:
        fh1.close()
        fh2.close()
    printime(' - Finished extracting %s' % file2)
    return sum_raw, sum_nrm, sqr_raw, sqr_nrm, passage


def mean_metamatrix(sum_raw, sum_nrm, sqr_raw, sqr_nrm, passage, outdir, label):
    '''To obtain the mean matrix, divide raw and norm per pasages, plot if wanted'''
    size = max(passage.keys())[0] + 1

    array_avg_raw = np.zeros((size, size))
    array_avg_nrm = np.zeros((size, size))
    array_std_raw = np.zeros((size, size))
    array_std_nrm = np.zeros((size, size))
    array_passage = np.zeros((size, size))

    for x in sum_raw:
        N = float(passage[x])
        array_passage[x] = N
        array_avg_raw[x] = sum_raw[x] / N
        array_avg_nrm[x] = sum_nrm[x] / N
        array_std_raw[x] = (sqr_raw[x] / N - (sum_raw[x] / N)**2)**0.5
        array_std_nrm[x] = (sqr_nrm[x] / N - (sum_nrm[x] / N)**2)**0.5

    np.savetxt(os.path.join(outdir, 'mean_raw_' + label + '.txt'), array_avg_raw)
    np.savetxt(os.path.join(outdir, 'mean_nrm_' + label + '.txt'), array_avg_nrm)
    np.savetxt(os.path.join(outdir, 'stdv_raw_' + label + '.txt'), array_std_raw)
    np.savetxt(os.path.join(outdir, 'stdv_nrm_' + label + '.txt'), array_std_nrm)
    np.savetxt(os.path.join(outdir, 'passages_' + label + '.txt'), array_passage)


def main():
    opts = get_options()

    resolution   = opts.resolution
    peak_file    = opts.peak_file
    tmpdir       = opts.tmpdir
    outdir       = opts.outdir
    ncpus        = opts.ncpus
    biases       = opts.biases
    genomic_mat  = opts.genomic_mat
    inbam        = opts.inbam

    badcols = load(open(biases))['badcol']

    if not tmpdir:
        tmpdir = os.path.join(outdir, 'tmp')

    mkdir(outdir)
    mkdir(tmpdir)

    # get chromosome lengths
    bamfile  = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references,
                               [x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    # given sublist as input
    sublists = glob(os.path.join(peak_file, '*'))
    peaks = {}
    for sublist in sublists:
        peak_fname = os.path.split(sublist)[-1]
        peaks[os.path.split(sublist)[-1]] = {
            'ori_path': sublist,
            'tmp_path': os.path.join(tmpdir, peak_fname + '_tmp')}

    # extract all sublists of peaks
    pool = Pool(ncpus)
    printime('Extracting coordinates')
    for peak in peaks:
        pool.apply_async(extract_coordinates,
                         args=(peaks[peak]['ori_path'], peaks[peak]['tmp_path'],
                               resolution, section_pos, badcols))
    pool.close()
    pool.join()

    # sort all sublists of peaks
    printime('Sorting coordinates')
    for peak in peaks:
        peaks[peak]['sorted'] = peaks[peak]['tmp_path'] + '_sorted'
        Popen(("sort -k1,2n -S 20% --parallel={0} {1} "
               "--temporary-directory={2} > {3}").format(
                   opts.ncpus, peaks[peak]['tmp_path'], tmpdir,
                   peaks[peak]['sorted']), shell=True).communicate()

    # extract submatrices
    pool = Pool(ncpus)
    printime('Getting interactions')
    procs = {}
    for peak in peaks:
        readfiles(genomic_mat, peaks[peak]['sorted'])
        procs[peak] = pool.apply_async(readfiles, (genomic_mat, peaks[peak]['sorted']))
    pool.close()
    pool.join()

    # # clean
    # for peak in peaks:
    #     os.system('rm -f '  + peaks[peak]['tmp_path'])
    # os.system('rm -rf ' + tmpdir)

    # save meta-waffles
    for peak in peaks:
        try:
            sum_raw, sum_nrm, sqr_raw, sqr_nrm, passage = procs[peak].get()
        except TypeError:
            continue
        if passage:
            printime(' - Generating meta-waffle %s...' % peak)
            # get mean matrix, raw and norm divided pasages
            mean_metamatrix(sum_raw, sum_nrm, sqr_raw, sqr_nrm, passage, outdir, peak)
            # save dicts
            out = open(os.path.join(outdir, 'submats_%s.pickle' % (peak)), 'wb')
            dump({'passage': passage,
                  'sum_raw': sum_raw, 'sum_nrm': sum_nrm,
                  'sqr_raw': sqr_raw, 'sqr_nrm': sqr_nrm}, out)
            out.close()
        else:
            print('No information at this interval: ', peak)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', '--peak', dest='peak_file', required=True, default=False,
                        help='''Pairwise peaks to compute average submatrix (norm and raw)''')
    parser.add_argument('-o', '--outdir', dest='outdir', required=True,
                        metavar='PATH', help='output directory')
    parser.add_argument('-bam', '--bam', dest='inbam', required=True, default=False,
                        help='Input HiC-BAM file')
    parser.add_argument('-b', '--biases', dest='biases', default=True, help='Biases',
                        required=True)
    parser.add_argument('-m', '--genomic_mat', dest='genomic_mat', default=True,
                        metavar='GENOMIC_MATRIX', required=True,
                        help='''Path to genomic matrix in 3 columns format
                        (should be sorted with `sort -k1,2n`)''')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-t', '--tmp', dest='tmpdir', default=None,
                        help='Tmpdir to store coordinates files')
    parser.add_argument('-C', '--cpus', dest='ncpus', type=int, default=8,
                        help='''[%(default)s] number of cpus to be used for parsing
                        the HiC-BAM file''')

    opts = parser.parse_args()

    return opts

if __name__ == '__main__':
    exit(main())
