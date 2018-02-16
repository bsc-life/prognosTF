from argparse                     import ArgumentParser
import pysam
from collections import OrderedDict
from datetime import datetime
import itertools
import glob


def read_line(line):   ##Get information per peak of a feature +/-
   c, p1, p2, f = line.split()[:4]
   return c, int(p1), int(p2), f

def intervals_to_windows(intervals,windows):
    return [windows[k-1] for k in intervals]

def binning_bed(peak_file, resolution, windows_span, max_dist, outdir,
                name, chrom_sizes, **kwargs):
    '''Input bed file of ChIP peaks and bin into desired resolution of Hi-C'''

    peaks = open(peak_file,"r")

    bin_coordinate = set((c, (p1 + p2) / 2 / resolution, f)
                         for c, p1, p2, f in map(read_line, peaks)
                            if (p1 + p2) / 2 > windows_span)    ## take into account to add windows span both sides


    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Total of different bin coordinates: ',len(bin_coordinate)

    wsp = ((windows_span * 2) / resolution) + 1
    mdr = max_dist / resolution

    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist
    # - Different combination of the features

    pairs = ([a, b] for a,b in itertools.combinations(bin_coordinate, 2)
             if a[0] == b[0] and wsp <= abs(b[1] - a[1]) <= mdr)

    windows = ((255000 , 1000000),
              (1000000 , 2500000),
              (2500000 , 5000000),
              (5000000 , 10000000),
              (10000000, 15000000),
              (15000000, 20000000),
              (20000000, 30000000))

    windows_intervals = kwargs.get('windows_intervals', None)

    if windows_intervals != None:
        windows = intervals_to_windows(windows_intervals,windows)

    personalized_intervals = kwargs.get('personalized_intervals', None)
    if personalized_intervals != None: # give the opportunity to define desired intervals of distance to study
        file_intervals = open(personalized_intervals,'r')
        list_windows = []
        for line in file_intervals:
            a, b = map(int, line.split())
            list_windows.append((a,b))
        windows = tuple(list_windows)


    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'- Writing list of coordinates...'

    intervals = dict(((b, e), []) for b, e in (windows))

    for (chromosome1, bs1, f1),(chromosome2,  bs2, f2) in pairs:
       distance = abs(bs1 - bs2)
       for (lower,upper) in windows:
           if lower / resolution < distance <= upper / resolution:
               intervals[(lower, upper)].append((chromosome1, bs1, f1, chromosome2, bs2, f2))

    adding = windows_span / resolution

    for beg, end in intervals:
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'Writing interval: ', beg, end
        w = open('%s%s_%d_%d.tsv'%(outdir, name, beg, end), 'wa')
        for line in intervals[(beg,end)]:
            c1, s1, f1, c2, s2, f2  = line
            start1, end1 = s1 - adding, s1 + adding
            start2, end2 = s2 - adding, s2 + adding
            #  check chromosome length
            new_start1, new_end1 = start1 * resolution, end1 * resolution
            new_start2, new_end2 = start2 * resolution, end2 * resolution
            if new_end1 > chrom_sizes[c1] or new_end2 > chrom_sizes[c1]:
                continue
            else:
                w.write('%s:%s-%s\t%s\t%s:%s-%s\t%s\n' % (c1, str(new_start1), str(new_end1),f1, c2, str(new_start2), str(new_end2),f2))
        w.close()

    ## open files splitted by distance and split again but depending different combinations
    intervals_files = glob.glob(outdir+'*.tsv')

    for file in intervals_files:
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'Splitting file: ', file
        name_id = file.split('/')[-1].split('.')[0]
        same_pos = open('%s%s_both_+.tsv'%(outdir,name_id),'wa')
        same_neg = open('%s%s_both_-.tsv'%(outdir, name_id),'wa')
        same_posneg = open('%s%s_+_-.tsv'%(outdir, name_id),'wa')
        same_negpos = open('%s%s_-_+.tsv'%(outdir, name_id),'wa')
        r = open(file,'r')
        for line in r:
            pair1, f1, pair2, f2 = line.split()
            if f1 == f2:
                if f1 == '+':
                    same_pos.write('%s\t%s\n' % (pair1, pair2))
                else:
                    same_neg.write('%s\t%s\n' % (pair1, pair2))
            else:
                if f1 == '+':
                    same_posneg.write('%s\t%s\n' % (pair1, pair2))

                else:
                    same_negpos.write('%s\t%s\n' % (pair1, pair2))
        same_pos.close()
        same_neg.close()
        same_negpos.close()
        same_posneg.close()


def main():
    opts = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_file    = opts.peak_file
    outdir       = opts.outdir
    name         = opts.name
    windows_span = opts.windows_span
    max_dist     = opts.max_dist

    if opts.windows_intervals:
        windows_intervals = opts.windows_intervals
    else:
        windows_intervals = False

    if opts.personalized_intervals:
        personalized_intervals = opts.personalized_intervals
    else:
        personalized_intervals = False


    ## peaks file sorted per chromosome
    bamfile = pysam.AlignmentFile(inbam,'rb')
    sections = OrderedDict(zip(bamfile.references,[x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    chrom_sizes = OrderedDict(zip(bamfile.references, [x for x in bamfile.lengths]))

    if opts.windows_intervals != False:
        print 'Only extracting selected intervals'
        binning_bed(peak_file, resolution,windows_span,max_dist,outdir,name,chrom_sizes,windows_intervals=windows_intervals)
    else:
        if opts.personalized_intervals != False :
            print 'Analyzing desired intervals of distances...'
            binning_bed(peak_file,resolution,windows_span,max_dist,outdir,name,chrom_sizes,personalized_intervals=personalized_intervals)
        else:
            binning_bed(peak_file,resolution,windows_span,max_dist,outdir,name,chrom_sizes)

    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'Sublists written!'


def get_options():
    parser = ArgumentParser(usage="-i Peaks -r INT [options]")

    parser.add_argument('-i','--peak', dest='peak_file',required=True, default=False,
                        help='''Pairwise peaks to compute average submatrix (norm and raw)''')
    parser.add_argument('-bam','--bam',dest='inbam',required=True, default=False,
                        help= 'Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-o','--outdir',dest='outdir',default=True,help='output directory')
    parser.add_argument('-n','--name',dest='name',default=True, help = 'Output name')
    parser.add_argument('-w','-w', dest='windows_span',required=True, default=False,type=int,
                        help='''Windows span around center of the peak''')
    parser.add_argument('-m','-max', dest='max_dist',required=True, default=False,type=int,
                        help='''Max dist between center peaks''')
    parser.add_argument('-windows_intervals','-windows_intervals', dest='windows_intervals',required=False,
                        default=False,type=int,nargs="+",
                        help='''If only interested in some intervals to check:
                        1,2,3,4,5,6,7" correspond to different windows''')
    parser.add_argument('-personalized_intervals','-personalized_intervals', dest='personalized_intervals',required=False,
                        default=False,help='''If want to study specific intervals of distances, tsv file with them''')

    opts = parser.parse_args()

    return opts

if __name__=='__main__':
    exit(main())