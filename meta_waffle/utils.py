"""
"""

import os
import errno
from datetime    import datetime
from time        import time
from collections import OrderedDict

try:
    from pysam       import AlignmentFile
except ImportError:
    print('WARNING: pysam not found.')


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
    try:  # python 2
        bamfile = AlignmentFile(inbam, 'rb')
    except:   # python 3
        bamfile = AlignmentFile(inbam, 'r')
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
