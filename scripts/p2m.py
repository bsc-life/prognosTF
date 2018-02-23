from collections import OrderedDict
from argparse    import ArgumentParser
from datetime    import datetime
from time        import sleep
import subprocess
import glob
import multiprocessing as mu
from functools import partial
import os
import errno
from collections import defaultdict
from cPickle import dump, load

import numpy as np
from pysam import AlignmentFile


def mkdir(dnam):
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(dnam):
            raise

def extract_coordinates(peak_fname, resolution, section_pos, tmpdir, badcols):
    '''Chunk file into multiple, and write them in parallel per file write coord of 10,000 peak pairs
    Write a dictionary depending of pairs of peak per target'''
    peak_tmp_fname = os.path.split(peak_fname)[-1]
    out = open(os.path.join(tmpdir, peak_tmp_fname + '_tmp'), 'w')
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

        for x, p1 in enumerate(xrange(start_bin1, end_bin1)):
            if p1 in badcols:
                continue
            for y, p2 in enumerate(xrange(start_bin2, end_bin2)):
                if p2 in badcols:
                    continue
                out.write('{}\t{}\t{}\t{}\n'.format(p1, p2, x, y))
    out.close()


def readfiles(file1, file2):
    def split_line1(l):
        a, b, c, d = l.split()
        return (int(a), int(b)), c, d
    def split_line2(l):
        a, b, c, d = map(int, l.split())
        return (a, b), c, d
    avg_raw  = defaultdict(int)
    avg_nrm  = defaultdict(float)
    avg_pass = defaultdict(int)
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), ' - Reading BAM and peaks from %s...' % file2
    fh1 = open(file1)
    fh2 = open(file2)
    pos1, raw, nrm = split_line1(fh1.next())
    pos2, x, y = split_line2(fh2.next())
    try:
        while True:
            if pos2 > pos1:
                pos1, raw, nrm = split_line1(fh1.next())
            elif pos1 == pos2:
                avg_raw[x,y] += int(raw)
                avg_nrm[x,y] += float(nrm)
                avg_pass[x,y] += 1
                pos2_ = pos2
                pos2, x, y = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated
                    pos1, raw, nrm = split_line1(fh1.next())
            else:
                avg_pass[x,y] += 1
                pos2, x, y = split_line2(fh2.next())
    except StopIteration:
        fh1.close()
        fh2.close()
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Finished extracting %s' % file2
    return avg_raw, avg_nrm, avg_pass


def mean_metamatrix(avg_raw, avg_nrm, avg_pass, outdir, label):
    '''To obtain the mean matrix, divide raw and norm per passages, plot if wanted'''
    size = max(avg_pass.keys())[0] + 1

    array_raw = np.zeros((size, size))
    array_nrm = np.zeros((size, size))
    for x in avg_raw.keys():
        array_raw[x] = float(avg_raw[x]) / float(avg_pass[x])
        array_nrm[x] = float(avg_nrm[x]) / float(avg_pass[x])

    np.savetxt(os.path.join(outdir, 'mean_raw_' + label + '.txt'), array_raw)
    np.savetxt(os.path.join(outdir, 'mean_nrm_' + label + '.txt'), array_nrm)
    # if want to plot wait until all finishes (same colorbarscale)


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

    ## peaks file sorted per chromosome
    bamfile  = AlignmentFile(inbam, 'rb')
    sections = OrderedDict(zip(bamfile.references,[x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    # given sublist as input
    sublists = glob.glob(os.path.join(peak_file, '*'))
    peaks = {}
    for sublist in sublists:
        peaks[os.path.split(sublist)[-1]] = {'ori_path': sublist}

    # extract all sublists of peaks
    pool = mu.Pool(ncpus)
    for peak in peaks:
        pool.apply_async(extract_coordinates, 
                         args=(peaks[peak]['ori_path'], resolution, 
                         section_pos, tmpdir, badcols))
    pool.close()
    pool.join()

    # sort all sublists of peaks
    procs = []
    for peak in peaks:
        peaks[peak]['tmp_path'] = os.path.join(tmpdir, peak + '_tmp')
        peaks[peak]['sorted_path'] = peaks[peak]['tmp_path'] + '_sorted'
        procs.append(subprocess.Popen(("sort -k1,2n -S 20% --parallel={0} {1} "
                                       "--temporary-directory={2} > {3}").format(
                opts.ncpus, peaks[peak]['tmp_path'], tmpdir, 
                peaks[peak]['sorted_path']), shell=True))
    while procs:
        for p in procs:
            if p.poll() is None:
                continue
            if p.poll() == 0:
                procs.remove(p)
            else:
                print p.communicate()
                raise Exception('ERROR: problem with sorting\n')
        sleep(0.1)
    for peak in peaks:
        os.system('rm -f ' + peaks[peak]['tmp_path'])

    # extract submatrices
    avg_raw  = {}
    avg_nrm  = {}
    avg_pass = {}
    for peak in peaks:
        avg_raw [peak] = defaultdict(int)
        avg_nrm [peak] = defaultdict(float)
        avg_pass[peak] = defaultdict(int)
    
    pool = mu.Pool(ncpus)
    procs = {}
    for peak in peaks:
        procs[peak] = pool.apply_async(readfiles, (genomic_mat, peaks[peak]['sorted_path']))
    pool.close()
    pool.join()

    # save meta-waffles
    for peak in peaks:
        avg_raw, avg_nrm, avg_pass = procs[peak].get()

        if avg_pass:
            mean_metamatrix(avg_raw, avg_nrm, avg_pass, outdir, peak) # get mean matrix, raw and norm divided passages
            out_raw = open(os.path.join(outdir, 'raw_%s.pickle' % (peak)),'wb')
            out_nrm = open(os.path.join(outdir, 'nrm_%s.pickle' % (peak)),'wb')
            out_pas = open(os.path.join(outdir, 'pass_%s.pickle' % (peak)),'wb')
            dump(avg_raw,out_raw)
            dump(avg_nrm,out_nrm)
            dump(avg_pass,out_pas)
            out_raw.close()
            out_nrm.close()
            out_pas.close()
        else:
            print 'No information at this interval: ', peak

    # clean
    for peak in peaks:
        os.system('rm -f '  + peaks[peak]['tmp_path'])
        os.system('rm -rf ' + tmpdir)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', '--peak', dest='peak_file',required=True, default=False,
                        help='''Pairwise peaks to compute average submatrix (norm and raw)''')
    parser.add_argument('-o', '--outdir',dest='outdir',required=True, 
                        metavar='PATH', help='output directory')
    parser.add_argument('-bam', '--bam',dest='inbam',required=True, default=False,
                        help= 'Input HiC-BAM file')
    parser.add_argument('-b', '--biases',dest='biases',default=True, help = 'Biases', 
                        required=True)
    parser.add_argument('-m', '--genomic_mat', dest='genomic_mat',default=True, 
                        metavar='GENOMIC_MATRIX', required=True,
                        help = '''Path to genomic matrix in 3 columns format 
                        (should be sorted with `sort -k1,2n`)''')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-t', '--tmp',dest='tmpdir', default=None,
                        help='Tmpdir to store coordinates files')
    parser.add_argument('-C', '--cpus',dest='ncpus',type=int, default=8,
                        help='''[%(default)s] number of cpus to be used for parsing the HiC-BAM file''')

    opts = parser.parse_args()

    return opts

if __name__=='__main__':
    exit(main())
