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

    waffle_file = opts.peak_file
    output      = opts.outfile
    title       = opts.title

    waffle = load(open(waffle_file, 'rb'))

    group = ''

    print(waffle.keys())
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

    plt.figure(figsize=(14, 16))
    plt.subplots_adjust(left=0, bottom=0.07, right=0.92, top=0.9, wspace=0, hspace=0.3)

    if title:
        plt.suptitle(title, size=15)
    ## POLAR
    axr = plt.subplot2grid((3, 14), (0, 3), rowspan=2, colspan=11, polar=True)
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
    axl = plt.subplot2grid((3, 14), (2, 7), colspan=7)
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
                '95% Prediction band'],
               title='Spearman: {:.2f}\n  p-val: {:.3e}'.format(spear, pval),
               loc='upper left', frameon=False, bbox_to_anchor=[0.71, 1.37])
    axl.set_xlabel('Distance $d_i$ ($i$ from 0 to $N$) between peaks ($\sqrt{max_{0 \leq j \leq N}(d_j)} - \sqrt{d_i}$)')
    axl.set_ylabel('Normalized interactions')
    axl.set_xticks(size**0.5 - yticks**0.5)
    axl.set_xticklabels(['{}'.format(nicer(t * resolution)) for t in yticks], rotation=90)
    axl.set_title('Interactions vs distances', size=12)
    axl.set_xlim((min(x) - abs(max(x)-min(x))*0.01, max(x) + abs(max(x)-min(x))*0.01))
    axl.grid()

    plt.savefig(output, format=output.split('.')[-1])


def rotate(li, x):
    return li[-x % len(li):] + li[:-x % len(li)]


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', '--input', dest='peak_file', required=True,
                        metavar='PATH', help='''path to input pickle file''')
    parser.add_argument('-o', '--output', dest='outfile', required=True,
                        metavar='PATH', help='''path to output image (any format
                        based on file extension)''')
    parser.add_argument('--title', dest='title', default=None,
                        metavar='STR', help='''some quoted text to be used as
                        title for the plot''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
