from argparse                     import ArgumentParser
import pysam
from collections import OrderedDict
from datetime import datetime
import subprocess
import uuid
import glob
import multiprocessing as mu
import functools
import os
from collections import defaultdict
from cPickle import dump,load
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np


def extract_coordinates(peak_list,resolution,section_pos,tmpdir,chromosome,badcols):
    '''Chunk file into multiple, and write them in parallel per file write coord of 10,000 peak pairs
    Write a dictionary depending of pairs of peak per target'''

    unique_filename = str(uuid.uuid4())
    w = open(tmpdir+chromosome+'_{}.tmp'.format(unique_filename),'wa')
    for line in peak_list:
        peak1, peak2 = line.split()
        chr1, beg1, end1 = peak1.replace(':','-').split('-')
        chr2, beg2, end2 = peak2.replace(':','-').split('-')
        pos = section_pos[chr1][0]
        if beg1 < beg2:
            start_bin1 = pos + (int(beg1) / resolution)
            end_bin1 = pos + (int(end1) / resolution)

            start_bin2 = pos + (int(beg2) / resolution)
            end_bin2 = pos + (int(end2) / resolution)
        else:
            start_bin1 = pos + (int(beg2) / resolution)
            end_bin1 = pos + (int(end2) / resolution)

            start_bin2 = pos + (int(beg1) / resolution)
            end_bin2 = pos + (int(end1) / resolution)
        for x, p1 in enumerate(xrange(start_bin1, end_bin1+1)):
            for y, p2 in enumerate(xrange(start_bin2, end_bin2+1)):
                if p1 in badcols or p2 in badcols:
                    continue
                w.write('{}\t{}\t{}\t{}\n'.format(p1, p2, x, y))
    w.close()


def eq_pos(pos1, pos2):
    return pos1 == pos2


def greater_pos(pos1, pos2):
    return pos1 > pos2


def readfiles(file1,file2,chromosome):
    def split_line1(l):
        a, b, c, d = l.split()
        return (int(a), int(b)), int(c), float(d)
    def split_line2(l):
        a, b, c, d = map(int, l.split())
        return (a, b), c, d
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Reading BAM and peaks ...'
    fh1 = open(file1)
    fh2 = open(file2)
    pos1, raw, nrm = split_line1(fh1.next())
    pos2, x, y = split_line2(fh2.next())
    try:
        while True:
            if eq_pos(pos1, pos2):
                avg_raw[(x,y)] += raw
                avg_nrm[(x,y)] += nrm
                avg_pass[(x,y)] += 1
                pos2_ = pos2
                pos2, x, y = split_line2(fh2.next())
                if pos2_ != pos2:  # some cells in the peak file are repeated
                    pos1, raw, nrm = split_line1(fh1.next())
            elif greater_pos(pos1, pos2):
                avg_pass[(x,y)] += 1
                pos2, x, y = split_line2(fh2.next())
            else:
                pos1, raw, nrm = split_line1(fh1.next())

    except StopIteration:
        fh1.close()
        fh2.close()
    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Finished'


def read_line(line):
   c, p1, p2 = line.split()[:3]
   return c, int(p1), int(p2)


def intervals_to_windows(intervals,windows):
    return [windows[k-1] for k in intervals]


def binning_bed(peak_file, resolution, windows_span, max_dist, outdir,
                name, chrom_sizes, **kwargs):
    '''Input bed file of ChIP peaks and bin into desired resolution of Hi-C'''

    peaks = open(peak_file,"r")

    bin_coordinate = set((c, (p1 + p2) / 2 / resolution)
                         for c, p1, p2 in map(read_line, peaks)
                            if (p1 + p2) / 2 > windows_span)    ## take into account to add windows span both sides


    print datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Total of different bin coordinates: ',len(bin_coordinate)

    wsp = ((windows_span * 2) / resolution) + 1
    mdr = max_dist / resolution

    # get all combinations of bin peaks:
    # - same chromosomes
    # - Not overlapping peaks, windows span added to bin*2 (wsp)
    # - below max_dist

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

    for (chromosome1, bs1),(chromosome2,  bs2) in pairs:
       distance = abs(bs1 - bs2)
       for (lower,upper) in windows:
           if lower / resolution < distance <= upper / resolution:
               intervals[(lower, upper)].append((chromosome1, bs1, chromosome2, bs2))

    adding = windows_span / resolution

    for beg, end in intervals:
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),'Writing interval: ', beg, end
        w = open('%s%s_%d_%d.tsv'%(outdir, name, beg, end), 'wa')
        for line in intervals[(beg,end)]:
            c1, s1, c2, s2 = line
            start1, end1 = s1 - adding, s1 + adding
            start2, end2 = s2 - adding, s2 + adding
            #  check chromosome length
            new_start1, new_end1 = start1 * resolution, end1 * resolution
            new_start2, new_end2 = start2 * resolution, end2 * resolution
            if new_end1 > chrom_sizes[c1] or new_end2 > chrom_sizes[c1]:
                continue
            else:
                w.write('%s:%s-%s\t%s:%s-%s\n' % (c1, str(new_start1), str(new_end1), c2, str(new_start2), str(new_end2)))

        w.close()
def colorbar_vals(matrices):
    vmax = []
    vmin = []
    for mat in matrices:
        nrm = load(open(mat,'rb'))
        p = mat.replace('nrm','pass')
        pas = load(open(p,'rb'))
        arr = np.zeros((51,51))
        for x in nrm.keys():
            arr[x]= float(nrm[x])/float(pas[x])
        max_diff = max(abs(1-np.min(arr)), abs(1-np.max(arr)))
        mmax = 1 + max_diff
        mmin = 1 - max_diff
        vmax.append(mmax), vmin.append(mmin)
    return max(vmax),min(vmin)

def target_plot(matrix,target,vmax,vmin,distance,out):
    fig, ax = plt.subplots(figsize=(5,5))

    plt.imshow(matrix,cmap='RdBu_r',interpolation='none',vmin = vmin, vmax = vmax)
    plt.title('Normalized '+target+': '+distance,fontsize=15, position=(0.5, 1.05))

    ax.set_xticks([0,5,10,15,20,25,30,35,40,45,50])
    labels = ['-125','-100','-75','-50','-25',target,'25','50','75','100','125']
    ax.set_xticklabels(labels,rotation=90,fontsize=12)
    ax.set_yticks([0,5,10,15,20,25,30,35,40,45,50])
    ax.set_yticklabels(labels,fontsize=12)

    plt.xlabel('Bins from ChIP peak',fontsize=14)
    plt.ylabel('Bins from ChIP peak',fontsize=14)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="3%", pad=0.3)
    plt.colorbar(cax=cax)
    plt.tight_layout()
    #plt.savefig(out+target+'_'+distance+'.pdf')
    plt.show()


def mean_metamatrix(avg_raw, avg_nrm, avg_pass, outdir, label):
    '''To obtain the mean matrix, divide raw and norm per passages, plot if wanted'''
    size = max(avg_pass.keys())[0] + 1

    array_raw = np.zeros((size, size))
    array_nrm = np.zeros((size, size))
    for x in avg_raw.keys():
        array_raw[x] = float(avg_raw[x]) / float(avg_pass[x])
        array_nrm[x] = float(avg_nrm[x]) / float(avg_pass[x])

    np.savetxt(outdir+'mean_raw_'+label+'.txt',array_raw)
    np.savetxt(outdir+'mean_nrm_'+label+'.txt', array_nrm)
    # if want to plot wait until all finishes (same colorbarscale)


def main():
    opts = get_options()

    inbam        = opts.inbam
    resolution   = opts.resolution
    peak_file    = opts.peak_file
    tmpdir       = opts.tmpdir
    outdir       = opts.outdir
    ncpus        = opts.ncpus
    name         = opts.name
    biases       = opts.biases
    mats         = opts.mats
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

    if opts.specific:
        specific = opts.specific
    else:
        specific = False


    bamfile = pysam.AlignmentFile(inbam,'rb')
    sections = OrderedDict(zip(bamfile.references,[x / resolution + 1 for x in bamfile.lengths]))
    total = 0
    section_pos = dict()
    for crm in sections:
        section_pos[crm] = (total, total + sections[crm])
        total += sections[crm]

    chrom_sizes = OrderedDict(zip(bamfile.references, [x for x in bamfile.lengths]))


    if specific == False: # write sublists
        ## peaks file sorted per chromosome

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

        #get sublists in directory
        sublists = glob.glob(outdir+'*.tsv')

    else:
        # given sublist as input
        sublists = glob.glob(peak_file+'*')


    for l in sublists:
        label = l.split('/')[-1].split('.')[0]
        print  datetime.now().strftime('%Y-%m-%d %H:%M:%S'), 'Getting finger: ', l
        #split file  peaks per chromosome
        print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Splitting peak pairs per chromosome...'
        fh = subprocess.Popen("awk -F: '{print >> " + '"' +  tmpdir + '"' +  " $1; close($1)}' %s " % ((l)),
                                  shell=True)
        fh.communicate() # wait until the process finishes
        badcols = load(open(biases))['badcol']
        chromosomes_file = glob.glob(tmpdir+'*')
        global avg_raw
        avg_raw = defaultdict(int)
        global avg_nrm
        avg_nrm = defaultdict(float)
        global avg_pass
        avg_pass = defaultdict(int)

        for peak_list in chromosomes_file:
            chromosome = peak_list.split('/')[-1]
            print peak_list
            peakfile = open(peak_list,'r')
            total_lines = peakfile.readlines()
            starts = [x for x in range(0,len(total_lines),10000)] # chunk to do in paralel 10000 bins at time
            lines = []
            print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Generate: ',len(starts),' tmp files for ',chromosome
            for n,i in enumerate(starts):
                if i < starts[-1]:
                    lines.append(total_lines[i:starts[n+1]])
                else:
                    lines.append(total_lines[i:len(total_lines)])

            # parallel job to write coordinates, splitting peak file, write tmp files
            print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Writing coordinates files...'
            pool = mu.Pool(ncpus)
            pool.map(functools.partial(extract_coordinates, chromosome=chromosome,resolution=resolution,section_pos=section_pos,tmpdir=tmpdir,badcols=badcols), lines)
            pool.close()
            pool.join()
            tmp_chr = glob.glob(tmpdir+'%s_*'%chromosome)

            if len(tmp_chr) > 1:
                print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Merging and sorting...'
                for tmp in tmp_chr:
                    os.system("cat "+tmp+">> %s%s_m"%(tmpdir,chromosome))
                for tmp in tmp_chr:
                    os.system("rm "+tmp)
                fh = subprocess.Popen("sort -k1,2n --parallel=24 %s%s_m --temporary-directory=/scratch/tmp/ > %s%s_sorted"%(tmpdir,chromosome,tmpdir,chromosome),
                                      shell=True)
                fh.communicate()
                os.system("rm %s%s_m"%(tmpdir,chromosome))
            else:
                print datetime.now().strftime('%Y-%m-%d %H:%M:%S'),' - Sorting...'
                fh = subprocess.Popen("sort -k1,2n --parallel=24 %s --temporary-directory=/scratch/tmp/ > %s%s_sorted"%(tmp_chr[0],tmpdir,chromosome),
                                      shell=True)
                fh.communicate()
                os.system("rm "+tmp_chr[0])


            # read bam chromosome and peak file same time
            file1 = mats+'%s_bam_%ikb.tsv'%(chromosome,resolution/1000)
            file2 = '%s%s_sorted'%(tmpdir,chromosome)
            check_size = open(file2,'r')
            nsize = check_size.readlines()
            if len(nsize) > 0:
                readfiles(file1,file2,chromosome)
                os.system("rm %s%s_sorted"%(tmpdir,chromosome))
                os.system("rm %s%s"%(tmpdir,chromosome))
            else:
                print 'No information for ', chromosome, ' all badcolumns :('
                os.system("rm %s%s_sorted"%(tmpdir, chromosome))
                os.system("rm %s%s"%(tmpdir, chromosome))

        if len(avg_pass) != 0:
            mean_metamatrix(avg_raw, avg_nrm, avg_pass, outdir, label) # get mean matrix, raw and norm divided passages
            out_raw=open(outdir+'raw_%s.pickle'%(label),'wb')
            out_nrm = open(outdir+'nrm_%s.pickle'%(label),'wb')
            out_pas = open(outdir+'pass_%s.pickle'%(label),'wb')
            dump(avg_raw,out_raw)
            dump(avg_nrm,out_nrm)
            dump(avg_pass,out_pas)
            out_raw.close()
            out_nrm.close()
            out_pas.close()
        else:
            print 'No information at this interval: ', label


def get_options():
    parser = ArgumentParser(usage="-i Peaks -r INT [options]")

    parser.add_argument('-i','--peak', dest='peak_file',required=True, default=False,
                        help='''Pairwise peaks to compute average submatrix (norm and raw)''')
    parser.add_argument('-bam','--bam',dest='inbam',required=True, default=False,
                        help= 'Input HiC-BAM file')
    parser.add_argument('-r', '--resolution', dest='resolution', required=True, default=False,
                        type=int, help='wanted resolution from generated matrix')
    parser.add_argument('-t','--tmp',dest='tmpdir',required=True, default=False,
                        help='Tmpdir to store coordinates files')
    parser.add_argument('-o','--outdir',dest='outdir',default=True,help='output directory')
    parser.add_argument('-C','--cpus',dest='ncpus',type=int,default=8,
                        help='''[%(default)s] number of cpus to be used for parsing the HiC-BAM file''')
    parser.add_argument('-n','--name',dest='name',default=True, help = 'Output name')
    parser.add_argument('-b','--biases',dest='biases',default=True, help = 'Biases')
    parser.add_argument('-mats','--mats',dest='mats',default=True, help = 'Folder where matrices are located')
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
    parser.add_argument('-specific','-specific', dest='specific',required=False,
                        default=False,help='''If sublists already generate use directory with sublists as input in peak_file input''')


    opts = parser.parse_args()

    return opts

if __name__=='__main__':
    exit(main())
