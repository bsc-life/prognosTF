#! /usr/bin/env python
"""
"""

from argparse    import ArgumentParser
try:  # python 3
    from pickle        import load
except ImportError:  # python 2
    from cPickle        import load

from scipy import stats
import matplotlib.pyplot as plt
import numpy as np

from serutils.stats.correlate import fit_with_uncertainty, latex_formula


def nicer(res, sep=' ', comma='', allowed_decimals=0):
    """
    TODO: go to utils
    writes resolution number for human beings.

    :param ' ' sep: character between number and unit (e.g. default: '125 kb')
    :param '' comma: character to separate groups of thousands
    :param 0 allowed_decimals: if 1 '1900 kb' would be written as '1.9 Mb'
    """
    format = lambda x: '{:,g}'.format(x).replace(',', comma)

    if not res:
        return format(res) + sep + 'b'
    if not res % 10**(9 - allowed_decimals):
        return format(res / 10.**9) + sep + 'Gb'
    if not res % 10**(6 - allowed_decimals):
        return format(res / 10.**6) + sep + 'Mb'
    if not res % 10**(3 - allowed_decimals):
        return format(res / 10.**3) + sep + 'kb'
    return format(res) + sep + 'b'


def main():

    opts = get_options()

    waffle_files   = opts.peak_files
    waffle = load(open(waffle_files, 'rb'))

    group = ''

    resolution = waffle[group]['resolution']
    size       = waffle[group]['size']
    counter    = waffle[group]['counter']

    matrix = [[waffle[group]['sum_nrm'][i, j] / counter
               for i in range(size)]
              for j in range(size)]

    #get coordinates:
    Phi = []
    R = []
    data = []
    divs = 10000
    mid = size // 2


    for k in range(size):
        tmp = []
        for i in range(-mid, size + mid):
            di = abs(mid - i)
            for j in range(-mid, size + mid):
                dj = abs(mid - j)
                if di + dj != k:
                    continue
                if not len(tmp) % 2:
                    try:
                        if i * j < 0:
                            tmp.insert(0, float('nan'))
                        else:
                            tmp.insert(0, matrix[j][i])
                    except IndexError:
                        tmp.insert(0, float('nan'))
                else:
                    try:
                        if i * j < 0:
                            tmp.append(float('nan'))
                        else:
                            tmp.append(matrix[j][i])
                    except IndexError:
                        tmp.append(float('nan'))
        ntmp = len(tmp)
        Phi.append(np.linspace(0, np.pi * 2, divs + 1))
        R.append(np.linspace(k, k, divs + 1))
        aug = int(divs / ntmp)
        res = divs / ntmp - aug
        i = 0
        tmp2 = []
        for t in tmp:
            v = int(res * i)
            if v:
                i = res * i - v
            tmp2.append([t] * (aug + v))
            i += 1

        tmp2 = [v for l in tmp2 for v in l]
        tmp2 += [tmp2[-1]] * (divs - len(tmp2))
        tmp2 = rotate(tmp2, aug // 2)
        data.append(np.asarray(tmp2))

    Phi.append(np.linspace(0, np.pi * 2, divs + 1))
    R.append(np.linspace(k + 1, k + 1, divs + 1))
    Phi = np.asarray(Phi)
    R = np.asarray(R)
    data = np.asarray(data)

    mid = size // 2
    xvals = []
    yvals = []
    for i in range(size):
        di = abs(mid - i)
        for j in range(size):
            dj = abs(mid - j)
            xvals.append(di + dj)
            yvals.append(matrix[i][j])
    xvals = np.asarray(xvals)
    yvals = np.asarray(yvals)

    plt.figure(figsize=(20, 16))

    ## POLAR
    axr = plt.subplot2grid((3, 14), (0, 0), rowspan=2, colspan=13, polar=True)
    axr.set_title('Corrected by distance average submatrices between peaks', size=13)
    # axr = plt.subplot(211, polar=True)
    m = axr.pcolormesh(Phi, R, data, linewidth=0)
    plt.colorbar(m, shrink=0.4, label='Averaged normalized interactions')
    yticks = np.asarray([0] + axr.get_yticks())
    axr.plot([0, 0, -np.pi / 2], [size, 0, size], 'k--', alpha=0.4, lw=1)

    axr.text(0         , size, ' peak', va='center')
    axr.text(-np.pi / 2, size, '\npeak', ha='center', va='top')

    axr.set_yticks(yticks)
    axr.set_yticklabels(['{}'.format(nicer(v * resolution)) for v in yticks])
    axr.set_xticks([])
    axr.set_ylim(0, size)
    axr.grid(ls='--')
    axr.set_rlabel_position(45)

    ## SQUARE
    axs = plt.subplot2grid((3, 14), (2, 0), colspan=7)
    axs.set_title('Average submatrices between peaks', size=12)
    _ = axs.imshow(matrix, origin='lower')
    axs.plot([mid, mid, -0.5], [-0.5, mid, mid], ls='--', color='k', alpha=0.4, lw=1)
    xticks = axs.get_xticks()
    if not mid in xticks:
        xticks = list(xticks)
        xticks.insert(len(xticks) // 2, mid)
        xticks = np.asarray(xticks)
    axs.set_xticks(xticks)
    axs.set_yticks(xticks)
    axs.set_xticklabels(['{}'.format(nicer((v - mid) * resolution)) if v != mid else 'peak'
                         for v in xticks], rotation=90)
    axs.set_yticklabels(['{}'.format(nicer((v - mid) * resolution)) if v != mid else 'peak'
                         for v in xticks])
    axs.set_xlim(-0.5, size - 0.5)
    axs.set_ylim(-0.5, size - 0.5)

    ## CORRELATION
    axl = plt.subplot2grid((3, 14), (2, 7), colspan=4)
    x = size**0.5 - xvals**0.5
    y = yvals
    func_string = 'A*x + C'
    p_x, p_y, z, confs, preds, r2 = fit_with_uncertainty(x, y, func_string=func_string, use_odr=True)
    formula = latex_formula(func_string, z)
    dots = axl.plot(x, y, 'o', alpha=0.2, ms=3)
    _ = axl.plot(x, y, 'o', alpha=0.5, ms=3.5, mec='k', mfc='none', mew=0.5)
    fit_line = axl.plot(p_x, p_y,color= 'firebrick', lw=2, label='Regression line')
    axl.fill_between(p_x, p_y + confs, p_y - confs, color='firebrick', alpha=0.3)
    axl.fill_between(p_x, p_y - preds, p_y + preds, color='grey', alpha=0.2)
    p1 = plt.Rectangle((0, 0), 1, 1, fc="firebrick", alpha=.3)
    p2 = plt.Rectangle((0, 0), 1, 1, fc="grey", alpha=.2)
    spear, pval = stats.spearmanr(xvals, yvals)

    axl.legend(fit_line + [p1, p2],
               ['Fit ($R^2=%.3f$):\ny = $%s$' % (r2, formula),
                '95% Confidence band',
                '95% Prediction band\n'],
               title='Spearman: {:.2f}\n  p-val: {:.3e}'.format(spear, pval),
               loc='upper left', frameon=False, bbox_to_anchor=[1,1])
    axl.set_xlabel('Distance $d$ between peaks ($\sqrt{%d} - \sqrt{d}$)' % (size))
    axl.set_ylabel('Normalized interactions')
    axl.set_xticks(size**0.5 - yticks**0.5)
    axl.set_xticklabels(['{}'.format(nicer(t * resolution)) for t in yticks], rotation=90)
    axl.set_title('Correlation between interactions and distance', size=12)
    axl.set_xlim((min(x) - abs(max(x)-min(x))*0.01, max(x) + abs(max(x)-min(x))*0.01))
    axl.grid()


def rotate(li, x):
    return li[-x % len(li):] + li[:-x % len(li)]




def get_options():
    parser = ArgumentParser()

    parser.add_argument('--peaks', dest='peak_files', required=True,
                        nargs="+", metavar='PATH',
                        help='''one or two pairwise peaks files to
                        compute average submatrix (norm and raw). These files
                        should contain at least two columns, chromosome and
                        position, and may contain an extra column of feature. If
                        present, the result will be returned according to the
                        possible combination of this feature''')
    # parser.add_argument('--bam', dest='inbam', required=True,
    #                     metavar='PATH', help='Input HiC-BAM file')
    # parser.add_argument('--biases', dest='biases', default=True, help='Biases',
    #                     required=True)
    parser.add_argument('--genomic_matrix', dest='genomic_mat', default=True,
                        metavar='GENOMIC_MATRIX', required=True,
                        help='''Path to genomic matrix in 3 columns format
                        (should be sorted with `sort -k1,2n`)''')
    # parser.add_argument('-r', '--resolution', dest='resolution', required=True,
    #                     metavar='INT', default=False, type=int,
    #                     help='wanted resolution from generated matrix')
    parser.add_argument('-o', '--outfile', dest='outfile', default='',
                        metavar='PATH', help='path to output file (pickle format)')
    parser.add_argument('--all_submatrices', dest='submatrices', default='',
                        metavar='PATH', help='''if PATH is provided here, stores
                        all the individual submatrices generated''')
    parser.add_argument('-s', dest='windows_span', required=True, type=int,
                        metavar='INT',
                        help='''Windows span around center of the peak (in bins; the
                        total windows size is 2 times windows-span + 1)''')
    parser.add_argument('-m', '--max_dist', dest='max_dist', metavar='INT',
                        default=float('inf'), type=int,
                        help='''[%(default)s] Max dist between center peaks''')
    parser.add_argument('-w', '--window', dest='window', required=False,
                        default='intra', metavar='INT-INT', type=str,
                        help='''[%(default)s] If only interested in some
                        intervals to check: "-w 1000000-2000000"
                        correspond to the window interval, from 1Mb to 2Mb.
                        Use "-w inter" for inter-chromosomal regions, "-w intra" for
                        intra-chromosomal, "-w all" for all combinations
                        (without distance restriction)''')
    parser.add_argument('--first_is_feature', dest='first_is_feature', default=False,
                        action='store_true', help='''When 2 BED files are input,
                        the peaks in the first BED should be also considered as
                        feature. This is to create average sub-matrices for
                        each peak in the first BED file.''')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
