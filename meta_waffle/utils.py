"""
"""

import os
import errno
from datetime    import datetime
from time        import time
from collections import OrderedDict, defaultdict
from pickle      import dump, load
from __future__  import print_function

try:
    from pysam       import AlignmentFile
except ImportError:
    print('WARNING: pysam not found.')


def sum_2_waffles(waffle1, waffle2):
    size = waffle1['size']
    for k in ['sum_raw', 'sqr_raw', 'passage', 'sqr_nrm', 'sum_nrm']:
        pik = waffle1[k]
        pik2 = waffle2[k]
        waffle1[k] = defaultdict(int, (((p1, p2), pik[p1, p2] + pik2[p1, p2])
                                       for p1 in range(size)
                                       for p2 in range(size)))
    waffle1['counter'] += waffle1['counter']


def sum_groups(waffle_list, waffle_out, split_features=False, clean=False, verbose=True):
    """
    Merges several outputs from waffle_peaks into a single pickle file

    :param waffle_list: list of pahts
    :param waffle_out: path to output file
    :param False clean: to remove input files
    """
    pi = load(open(waffle_list[0], 'rb'))
    missed  = set()
    counter = dict((k, pi[k]['counter']) for k in pi)
    size    = dict((k, pi[k]['size']) for k in pi)
    total_files = len(waffle_list)
    if verbose:
        print('Summing sub-waffles...')
        print ('   - {:6,d} / {:6,d}'.format(1, total_files), end='\r', flush=True)
    for nfile, fnam in enumerate(waffle_list[1:], 2):
        try:
            pi2 = load(open(fnam, 'rb'))
        except:
            if verbose:
                print('missing', fnam)
            missed.add(fnam)
            continue
        if verbose:
            print ('   - {:6,d} / {:6,d}'.format(nfile, total_files), end='\r', flush=True)
        for key in pi2:
            if not key in pi:
                pi[key] = pi2[key]
                counter[key] = pi[key]['counter']
                size[key] = pi[key]['size']
                continue
            sum_2_waffles(pi[key], pi2[key])
            counter[key] += pi2[key]['counter']
    if verbose:
        print('\ndone.')
        print('Saving meta-waffle')
    if split_features:
        if not os.path.exists(waffle_out):
            os.makedirs(waffle_out)
        for feature in pi:
            out_file = os.path.join(waffle_out, feature + '.pickle')
            if os.path.exists(out_file):
                sum_2_waffles(pi[feature], load(open(out_file, 'rb')))
            pickle_out = open(out_file, 'wb')
            dump(pi[feature], pickle_out)
            pickle_out.close()
    else:
        pickle_out = open(waffle_out, 'wb')
        dump(pi, pickle_out)
        pickle_out.close()
    if clean:
        if verbose:
            print('Removing sub-waffles...')
        for fnam in waffle_list:
            if fnam in missed:
                continue
            os.remove(fnam)
    if verbose:
        print('\ndone.')
    return missed, counter


def parse_fasta(infasta):
    fh = open(infasta)
    genome = OrderedDict()
    lseq = 0
    chrom = None
    for line in fh:
        if line.startswith('>'):
            if chrom is not None:
                genome[chrom] = lseq
            chrom = line[1:].split()[0]
            lseq = 0
        else:
            lseq += len(line) - 1
    genome[chrom] = lseq
    return genome


def chromosome_from_bam(inbam, resolution, get_bins=False):
    bamfile = AlignmentFile(inbam, 'rb')
    chrom_sizes = OrderedDict()
    for i, c in enumerate(bamfile.references):
        chrom_sizes[c] = bamfile.lengths[i]
    return chromosome_from_header(chrom_sizes, resolution, get_bins=get_bins)


def chromosome_from_fasta(infasta, resolution, get_bins=False):
    chrom_sizes = parse_fasta(infasta)
    return chromosome_from_header(chrom_sizes, resolution, get_bins=get_bins)


def chromosome_from_header(chrom_sizes, resolution, get_bins=False):
    sections = OrderedDict((c, v // resolution + 1) for c, v in chrom_sizes.items())

    total = 0
    section_pos = dict()
    bins = {}
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        if get_bins:
            for n, i in enumerate(range(*section_pos[crm])):
                bins[i] = (crm, n)
        total += sections[crm]

    return section_pos, chrom_sizes, bins


def mkdir(dnam):
    dnam = os.path.abspath(dnam)
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(dnam):
            raise


def printime(msg, silent=True):
    if silent:
        return
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')
