"""
"""
from scipy import stats, odr
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
try:  # pysal 1.14
    from pysal import weights as pysal_weights
    from pysal.esda import moran
except ImportError:
    try:  # pysal 2.0
        from pysal.lib import weights as pysal_weights
        from pysal.explore.esda import moran
    except ImportError:
        print('WARNING: PYSAL not installed')

# from serutils.stats.correlate import fit_with_uncertainty, latex_formula


def get_confidence(x, y, p_y, new_x, df, conf=0.95):
    """
    Computes prediction band and confidence band from distribution and it's fit

    :params x: data array or list
    :params y: data array or list
    :params p_y: predicted data array corresponding to x
    :params new_x: data array to be predicted
    :params None df: number of degrees of freedom needed to fit
    :params .95 conf: desired confidence level, by default 0.95 (2 sigma)

    :returns: an array with the confidence values to add/subtract, another
       with prediction values, and the R-square value of the fit.
    """
    if not isinstance(x, np.ndarray):
        x = np.array(x)
    if not isinstance(y, np.ndarray):
        y = np.array(y)
    if not isinstance(p_y, np.ndarray):
        p_y = np.array(p_y)
    if not isinstance(new_x, np.ndarray):
        new_x = np.array(new_x)
    # number of samples in original fit
    n = x.size
    # alpha 1 minus the wanted probability
    alpha = 1. - conf
    # t distribution with n-2 degrees of freedom
    t = stats.t.ppf(1. - alpha / 2., n - df - 1)
    # mean of x
    mean_x = np.mean(x)
    # Error sum of squares
    sse = sum((y - p_y)**2)
    # Error mean square (estimate of the variance)
    mse = sse / (n - df - 1)
    # Square individual deviation
    sdi = (new_x - mean_x)**2
    # standard deviation
    sd = sum((x - mean_x)**2)
    # relative individual deviation
    sdi_sd = sdi / sd

    confs = t * np.sqrt(mse * (      1.0 / n + sdi_sd))
    preds = t * np.sqrt(mse * (1.0 + 1.0 / n + sdi_sd))

    # calculate R-square
    sst = sum((y - np.mean(y))**2)
    r2 = 1 - sse / sst

    return confs, preds, r2


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


def func(x, *z):
    a, b = z
    return a * x + b


def func_for_odr(z, x):
    a, b = z
    return a * x + b


def rotate(li, x):
    return li[-x % len(li):] + li[:-x % len(li)]


def plot_polar_waffle(matrix, size, divs=20000, resolution=1, axe=None, vmin=None, vmax=None):
    #get coordinates:
    Phi = []
    R = []
    data = []
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

    axr = axe if axe else plt.subplot()

    axr.set_title('Corrected by distance average submatrices between peaks', size=13)
    # axr = plt.subplot(211, polar=True)
    m = axr.pcolormesh(Phi, R, data, linewidth=0, vmin=vmin, vmax=vmax)
    plt.colorbar(m, shrink=0.4, label='Averaged normalized interactions')
    yticks = np.asarray([0] + axr.get_yticks())
    axr.plot([-np.pi, 0, -np.pi / 2], [size, 0, size], 'k--', alpha=0.4, lw=1)

    axr.text(-np.pi    , size, 'peak ' , ha='right' , va='center')
    axr.text(-np.pi / 2, size, '\npeak', ha='center', va='top')

    axr.set_yticks(yticks)
    axr.set_yticklabels(['{}'.format(nicer(v * resolution)) for v in yticks])
    axr.set_xticks([])
    axr.set_ylim(0, size)
    axr.grid(ls='--')
    axr.set_rlabel_position(45)


def correlate_distances(matrix, size, metric='loop'):
    mid = size // 2
    xvals = []
    yvals = []
    averages =[[], []]
    for i in range(size):
        di = abs(mid - i)
        for j in range(size):
            dj = abs(mid - j)
            xvals.append(di + dj)
            yvals.append(matrix[i][j])
            if i <= mid and j <= mid:
                averages[0].append(((i, j), xvals[-1], yvals[-1]))
            else:
                averages[1].append(((i, j), xvals[-1], yvals[-1]))
    xvals = np.asarray(xvals)
    yvals = np.asarray(yvals)

    spear, pval = stats.spearmanr(xvals, yvals)

    if metric == 'normal':
        x = size**0.5 - xvals**0.5
        x1 = x
        x2 = []
        y = (yvals - np.mean(yvals)) / np.std(yvals)
        y1 = y
        y2 = []
    elif metric == 'loop':
        x1 = [size**0.5 - v**0.5 for _, v, _ in averages[0]]
        y1 = np.asarray([v for _, _, v in averages[0]])
        y1 = (y1 - np.mean(y1)) / np.std(y1)

        x2 = [size**0.5 - v**0.5 for _, v, _ in averages[1]]
        y2 = np.asarray([v for _, _, v in averages[1]])
        y2 = (y2 - np.mean(y2)) / np.std(y2)

        x = np.asarray(x1 + x2)
        y = np.asarray(list(y1) + list(y2))

    z, _ = curve_fit(func, x, y, [1., 1.])

    linear_model = odr.Model(func_for_odr)
    data = odr.RealData(x, y, sy=np.var(x), sx=np.var(y))  # rescale in case X and Y have different ranges
    odr_obj = odr.ODR(data, linear_model, beta0=z)  # use previous fit as starting point
    out = odr_obj.run()
    z = out.beta

    x_range = min(x) - abs(max(x)-min(x))*0.01, max(x) + abs(max(x)-min(x))*0.01

    new_x = np.linspace(min(x) if x_range[0] is None else x_range[0],
                        max(x) if x_range[1] is None else x_range[1], 1000)

    # now predict y based on test x-values
    new_y = func(new_x, *z)

    # estimate confidence and prediction bands
    p_y = func(x, *z)
    confs, preds, r2 = get_confidence(x, y, p_y, new_x, 2, conf=0.95)

    return (spear, pval), (x, y, new_x, new_y, z, confs, preds, r2)


def plot_square_waffle(matrix, size, resolution=1, axe=None):
    mid = size // 2

    axs = axe if axe else plt.subplot()

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


def plot_waffle(waffle, title, output=None, plot=True, axe=None):

    resolution = waffle['resolution']
    size       = waffle['size']
    counter    = waffle['counter']

    matrix = [[waffle['sum_nrm'][i, j] / counter
               for i in range(size)[::-1]]
              for j in range(size)]

    if not plot:
        (spear, pval), (x, y, p_x, p_y, z, confs, preds, r2) = correlate_distances(
            matrix, size)
        return (spear, pval), z

    ## PLOT
    plt.figure(figsize=(14, 16))
    plt.subplots_adjust(left=0, bottom=0.07, right=0.92, top=0.9, wspace=0, hspace=0.3)

    if title:
        plt.suptitle(title, size=15)
    ## POLAR
    axr = plt.subplot2grid((3, 14), (0, 3), rowspan=2, colspan=11, polar=True)

    plot_polar_waffle(matrix, size, resolution=resolution, axe=axr)

    ## SQUARE
    axs = plt.subplot2grid((3, 14), (2, 0), colspan=7)
    plot_square_waffle(matrix, size, resolution, axe=axs)

    ## CORRELATION
    axl = plt.subplot2grid((3, 14), (2, 7), colspan=7)

    (spear, pval), (x, y, p_x, p_y, z, confs, preds, r2) = correlate_distances(matrix, size)
    plot_correlation(spear, pval, x, y, p_x, p_y, z, confs, preds,
                     r2, size, resolution=resolution, axe=axl)

    ## save
    if output:
        plt.savefig(output, format=output.split('.')[-1])

    return (spear, pval), z


def plot_correlation(spear, pval, x, y, p_x, p_y, z, confs, preds,
                     r2, size, resolution=1, axe=None):
    axl = axe if axe else plt.subplot()

    _ = axl.plot(x, y, 'o', alpha=0.2, ms=3)
    _ = axl.plot(x, y, 'o', alpha=0.5, ms=3.5, mec='k', mfc='none', mew=0.5)


    fit_line = axl.plot(p_x, p_y,color= 'firebrick', lw=2, label='Regression line')
    axl.fill_between(p_x, p_y + confs, p_y - confs, color='firebrick', alpha=0.3)
    axl.fill_between(p_x, p_y - preds, p_y + preds, color='grey', alpha=0.2)
    p1 = plt.Rectangle((0, 0), 1, 1, fc="firebrick", alpha=.3)
    p2 = plt.Rectangle((0, 0), 1, 1, fc="grey", alpha=.2)

    axl.legend(fit_line + [p1, p2],
               ['Fit ($R^2=%.3f$):\ny = $%.3fx + %.3f$' % (r2, z[0], z[1]),
                '95% Confidence band',
                '95% Prediction band'],
               title='Spearman: {:.2f}\n  p-val: {:.3e}'.format(spear, pval),
               loc='upper left', frameon=False, bbox_to_anchor=[0.71, 1.37])
    axl.set_xlabel(('Distance $d_i$ ($i$ from 0 to $N$) between peaks '
                    '($\sqrt{max_{0 \leq j \leq N}(d_j)} - \sqrt{d_i}$)'))
    axl.set_ylabel('Normalized interactions')
    maxv = int(max((size**0.5-p_x)**2))
    for nticks in range(10, 0, -1):
        if maxv / nticks == maxv // nticks:
            break
    yticks = np.linspace(0, maxv, nticks + 1)
    axl.set_xticks(size**0.5 - yticks**0.5)
    axl.set_xticklabels(['{}'.format(nicer(t * resolution)) for t in yticks], rotation=90)
    axl.set_title('Interactions vs distances', size=12)
    axl.set_xlim((min(x) - abs(max(x)-min(x))*0.01, max(x) + abs(max(x)-min(x))*0.01))
    axl.grid()


def get_MI(matrix, width=2, loop=False, seed=1):
    """
    Computing the Moran Index for each matrix cell
    -
    Moran Index is a measure of multi-dimensional spatial correlation.
    - A positive I value indicates that the feature is surrounded by features with similar values
      -> part of a cluster.
    - A negative I value indicates that the feature is surrounded by features with dissimilar values
      -> an outlier.
    The Local Moran's index can only be interpreted within the context of the computed Z score or p-value.
    """
    mi_stats= {}
    size = len(matrix)
    mi_stats['size'] = size
    matrixlog2 = np.log2(matrix)

    # Using random_seed to block the Moran Index one, and get consistency with all experiments
    np.random.seed(seed)

    n, ws = get_weights(matrix=matrix, size=size, width=width, loop=loop)

    w = pysal_weights.W(n, ws)
    lm = moran.Moran_Local([[matrixlog2[i][j] for i in range(size)]
                                       for j in range(size)],
                                      w, permutations=9999)

    mi_stats["moranI locals"] = lm.p_sim, lm.q

    gm = moran.Moran([[matrixlog2[i][j] for i in range(size)]
                                 for j in range(size)],
                                w, permutations=9999)

    mi_stats["moranI global"] = gm.I, (gm.VI_rand, gm.seI_rand, gm.z_rand, gm.p_rand)

    return mi_stats


def get_weights(matrix, size, width=2, loop=False):
    """
    Computing the weight of each matrix position in relation to the central cell,
    which is the maximum.
    """
    n = {}
    ws = {}
    hwidth = width / 2.
    ihwidth = 1 + 1. / hwidth
    # Going to negative values; degree of correlation + non-correlation
    if loop:
        calculus = lambda x, y: (ihwidth - x / hwidth) * y
    # Using only positive values; degree of correlation
    else:
        calculus = lambda x, y: 1 / x * y
    for i in range(size):
        for j in range(size):
            b = i + j * size
            n[b] = []
            ws[b] = []
            val = matrix[i][j]
            for k in range(-width, width + 1):
                if i + k >= size or i + k < 0:
                    continue
                for l in range(-width, width+1):
                    if j + l >= size or j + l < 0 or k==l==0:
                        continue
                    c = i + k + (j + l) * size
                    n[b].append(c)
                    ws[b].append(calculus(max(abs(k), abs(l)), val))
    return n, ws


def plot_MoranI(waffle, width=2, plot=True, output=None, axe=None):
    """
    Visualizing the computed Moran Index with average submatrices between peaks
    -
    The COType field distinguishes between a statistically significant (0.05 level) values.
        HH = A cluster of high values
        LL = A cluster of low values
        HL = An outlier in which a high value is surround primarily by low values
        LH = An outlier in which a low value is surrounded primarily by high values
    """

    ## MORAN INDEX
    size       = waffle['size']
    counter    = waffle['counter']

    matrix = [[waffle['sum_nrm'][i, j] / counter
               for i in range(size)]
              for j in range(size)]

    mi_stats = get_MI(matrix=matrix, width=width, loop=False)
    Moran_I, (VI_rand, seI_rand, z_rand, p_rand) = mi_stats['moranI global']

    if not plot:
        return Moran_I, (VI_rand, seI_rand, z_rand, p_rand)

    ## PLOT
    plt.figure(figsize=(10, 8))
    axl = axe if axe else plt.subplot()

    axl.set_title('Moran index in average submatrices between peaks')
    colors = ['firebrick', 'mediumturquoise', 'royalblue', 'lightsalmon'] # 1 HH, 2 LH, 3 LL, 4 HL
    x = []
    y = []
    c = []
    size = mi_stats['size']

    for pvc in [0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]:
        for k, pv in enumerate(mi_stats['moranI locals'][0]):
            if pv > pvc:
                continue
            i, j = divmod(k, size)
            x.append(i)
            y.append(j)
            c.append(colors[mi_stats['moranI locals'][1][k]-1])

    im = axl.imshow(np.log2(matrix), interpolation='None', origin='lower', cmap='Greys')
    axl.scatter(x, y, alpha=0.15, color=c)

    red_patch = mpatches.Patch(color=colors[0], label='High values')
    cyan_patch = mpatches.Patch(color=colors[1], label='Low value is surrounded by high values')
    blue_patch = mpatches.Patch(color=colors[2], label='Low values')
    orange_patch = mpatches.Patch(color=colors[3], label='High value is surround by low values')

    axl.legend(handles=[red_patch, blue_patch, orange_patch, cyan_patch], ncol=2,
               loc='upper center', bbox_to_anchor=(0.5, -0.075), frameon=False,
               title='Global Moran index: {:.2f}   p-val (rand): {:.3e}'.format(Moran_I, p_rand))
    axl.set_xlim(-0.5, size - 0.5)
    axl.set_ylim(-0.5, size - 0.5)
    plt.colorbar(im, ax=axl)

    ## SAVE
    if output:
        plt.savefig(output, format=output.split('.')[-1])
