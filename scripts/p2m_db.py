from collections     import OrderedDict, defaultdict
from argparse        import ArgumentParser
from datetime        import datetime
from time            import time
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
                        badcols, ncpus, chunk):
    '''
    Chunk file into multiple, and write them in parallel per file write coordinates
    of 10,000 peak pairs Write a dictionary depending of pairs of peak per
    target
    '''
    max_line_size = 200  # space to search for the next "\n"
    fhandler = open(peak_fname)
    fhandler.seek(0, 2)
    fsize = fhandler.tell()
    fhandler.seek(0)

    beg = int(fsize / ncpus) * chunk
    if beg:
        for pos in range(max_line_size):
            fhandler.seek(beg + pos)
            if fhandler.read(1) == '\n':
                beg += pos + 1
                break

    end = int(fsize / ncpus) * (chunk + 1)
    for pos in range(max_line_size):
        fhandler.seek(end + pos)
        if fhandler.read(1) == '\n':
            end += pos + 1
            break
    if beg == end:
        raise Exception('ERROR: more CPUs than lines... perhaps no need to'
                        ' do it in parallel...')
    fhandler.seek(beg)

    out = open(tmp_fname, 'w')
    read = beg
    groups = {}  # to store groups to be made for the "meta" part
    for line in fhandler:
        chr1, beg1, end1, chr2, beg2, end2 = line.split()
        beg1, end1, beg2, end2 = int(beg1), int(end1), int(beg2), int(end2)
        if beg1 > beg2:
            beg1, end1, beg2, end2 = beg2, end2, beg1, end1
        # convert to genomic matrix coordinates (with all chromosomes concatenated)
        pos1 = section_pos[chr1][0]
        pos2 = section_pos[chr2][0]
        start_bin1 = pos1 + (beg1 / resolution)
        end_bin1   = pos1 + (end1 / resolution) + 1
        start_bin2 = pos2 + (beg2 / resolution)
        end_bin2   = pos2 + (end2 / resolution) + 1

        range1 = ((x, p1) for x, p1 in enumerate(range(start_bin1, end_bin1))
                  if p1 not in badcols)
        range2 = [(y, p2) for y, p2 in enumerate(range(start_bin2, end_bin2))
                  if p2 not in badcols]
        groups.setdefault(read, []).append((chr2, beg2))
        if range1 and range2:
            out.write('\n'.join('{}\t{}\t{}\t{}\t{}'.format(
                p1, p2, x, y, read)
                                for x, p1 in range1
                                for y, p2 in range2) + '\n')
        read += len(line)
        if read >= end:
            break
    out.close()

    out = open(tmp_fname + '_conv', 'w')
    dump(groups, out)
    out.close()


def readfiles(file1, file2, tmp_fname):
    def split_line1(line):
        a, b, c, d = line.split()
        return (int(a), int(b)), c, d
    def split_line2(line):
        a, b, c, d, e = map(int, line.split())
        return (a, b), c, d, e  # e is the position in the input file

    # create empty meta-waffles
    printime(' - Reading genomic matrix and peaks from:\n     -> %s' % file2)
    fh1 = open(file1)
    try:
        fh2 = open(file2)
    except IOError:
        return
    pos1, raw, nrm = split_line1(next(fh1))
    pos2, x, y, e = split_line2(next(fh2))

    out = open(tmp_fname, 'w')
    try:
        while True:
            if pos2 > pos1:
                pos1, raw, nrm = split_line1(fh1.next())
            elif pos1 == pos2:
                raw = int(raw)
                nrm = float(nrm)
                out.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                    pos1[0], pos1[1], x, y, raw, nrm, e))
                pos2_ = pos2
                pos2, x, y, e = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated
                    pos1, raw, nrm = split_line1(fh1.next())
            else:
                pos2, x, y, e = split_line2(fh2.next())
    except StopIteration:
        fh1.close()
        fh2.close()
    out.close()
    printime(' - Finished extracting %s' % file2)


def mean_metamatrix(sum_raw, sum_nrm, sqr_raw, sqr_nrm, passage, outdir, label):
    '''To obtain the mean matrix, divide raw and norm per passages, plot if wanted'''
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

    resolution  = opts.resolution
    peak_file   = opts.peak_file
    tmpdir      = opts.tmpdir
    outdir      = opts.outdir
    ncpus       = opts.ncpus
    biases      = opts.biases
    genomic_mat = opts.genomic_mat
    inbam       = opts.inbam

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
    bins = {}
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        for n, i in enumerate(range(*section_pos[crm])):
            bins[i] = (crm, n)
        total += sections[crm]

    # given path to file as input
    peaks = {}
    for chunk in range(ncpus):
        peak_fname = os.path.split(peak_file)[-1]
        peaks[chunk] = {
            'ori_path': peak_file,
            'tmp_path': os.path.join(tmpdir, '{}_{}tmp'.format(peak_fname, chunk))}

    # extract all sub-lists of peaks and translate into genomic matrix coordinates
    pool = Pool(ncpus)
    printime('Extracting coordinates')
    for chunk in range(ncpus):
        # extract_coordinates(peak_file, peaks[chunk]['tmp_path'],
        #                     resolution, section_pos, badcols, ncpus, chunk)
        # break
        if opts.force or not os.path.exists(peaks[chunk]['tmp_path']):
             pool.apply_async(extract_coordinates,
                              args=(peak_file, peaks[chunk]['tmp_path'],
                                    resolution, section_pos, badcols, ncpus, chunk))
    pool.close()
    pool.join()

    # sort all sub-lists of peaks (sort is already multiprocessing)
    printime('Sorting coordinates')
    skipped = True
    for chunk in peaks:
        peaks[chunk]['sorted'] = peaks[chunk]['tmp_path'] + '_sorted'
        if opts.force or not os.path.exists(peaks[chunk]['tmp_path'] + '_sorted'):
            Popen(("sort -k1,2n -S 10% --parallel={0} {1} "
                   "--temporary-directory={2} > {3}").format(
                       ncpus, peaks[chunk]['tmp_path'], tmpdir,
                       peaks[chunk]['sorted']), shell=True).communicate()
            skipped = False

    # extract sub-matrices of interactions
    pool = Pool(ncpus)
    printime('Getting interactions')
    for chunk in peaks:
        pool.apply_async(
            readfiles, (genomic_mat, peaks[chunk]['sorted'],
                        peaks[chunk]['tmp_path']))  # overwrite tmp_file
    pool.close()
    pool.join()

    # clean
    # for chunk in peaks:
    #     os.system('rm -f '  + peaks[chunk]['tmp_path'])
    # os.system('rm -rf ' + tmpdir)

    # sort by sub-matrix (7th column is the position of the sub-matrix in input)
    if opts.force or not skipped:
        printime('Sorting sub-matrices')
        for chunk in peaks:
            Popen(("sort -k7n -S 10% --parallel={0} {1} "
                   "--temporary-directory={2} > {3}").format(
                       ncpus, peaks[chunk]['tmp_path'], tmpdir,
                       peaks[chunk]['sorted']), shell=True).communicate()

    printime('Extracting individual information')
    out = open(os.path.join(outdir, 'final_per_cell.tsv'), 'w')
    for chunk in peaks:
        for line in open(peaks[chunk]['tmp_path'] + '_sorted'):
            X, Y, x, y, raw, nrm, read = line.split()
            c1, b1 = bins[int(X)]
            c2, b2 = bins[int(Y)]
            out.write(''.join('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
                c1, b1, c2, b2, x, y, raw, nrm, read)))
    out.close()

    # sort by sub-matrix
    printime('Sorting again sub-matrices')
    Popen(("sort -k9n -S 10% --parallel={0} {1} "
           "--temporary-directory={2} > {3}").format(
               ncpus, os.path.join(outdir, 'final_per_cell.tsv'), tmpdir,
               os.path.join(outdir, 'final_per_cell_sorted.tsv')),
          shell=True).communicate()

    printime('Summing peaks')
    # initialize defaultdict and prev variable
    sum_raw = defaultdict(int)
    sqr_raw = defaultdict(int)
    sum_nrm = defaultdict(float)
    sqr_nrm = defaultdict(float)
    passage = defaultdict(int)
    fhandler = open(os.path.join(outdir, 'final_per_cell_sorted.tsv'))
    line = next(fhandler)
    c1, b1, c2, b2, x, y, raw, nrm, read = line.split()
    prev = read

    fhandler.seek(0)
    out = open(os.path.join(outdir, 'final_sum.pickle'), 'w')
    for line in fhandler:
        c1, b1, c2, b2, x, y, raw, nrm, read = line.split()
        if read != prev:
            c, b, e = __get_original_coord(int(prev), peak_file)
            dump([(c, b, e), (sum_raw, sqr_raw, sum_nrm, sqr_nrm, passage)],
                 out)
            sum_raw = defaultdict(int)
            sqr_raw = defaultdict(int)
            sum_nrm = defaultdict(float)
            sqr_nrm = defaultdict(float)
            passage = defaultdict(int)
            prev = read
        raw = int(raw)
        nrm = float(nrm)
        x = int(x)
        y = int(y)
        sum_raw[x, y] += raw
        sum_nrm[x, y] += nrm
        sqr_raw[x, y] += raw**2
        sqr_nrm[x, y] += nrm**2
        passage[x, y] += 1
    out.close()


def __get_original_coord(pos, fname):
    fh = open(fname)
    fh.seek(pos)
    c1, b1, e1, c2, b2, e2 = next(fh).split()
    return c2, b2, e2


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
    parser.add_argument('--force', dest='force', default=False, action='store_true',
                        help='overwrites existing files')
    parser.add_argument('-t', '--tmp', dest='tmpdir', default=None,
                        help='Tmpdir to store coordinates files')
    parser.add_argument('-C', '--cpus', dest='ncpus', type=int, default=8,
                        help='''[%(default)s] number of cpus to be used for parsing
                        the HiC-BAM file''')

    opts = parser.parse_args()

    return opts

if __name__ == '__main__':
    exit(main())
