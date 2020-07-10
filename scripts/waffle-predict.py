#! /usr/bin/env python
"""
"""

from argparse          import ArgumentParser
try:  # python 3
    from pickle        import _Unpickler as Unpickler
except ImportError:  # python 2
    from pickle        import Unpickler

from scipy.stats       import spearmanr

from meta_waffle.stats import matrix_to_decay, get_center


def main():
    opts = get_options()

    waffle_file = opts.infile
    output      = opts.outfile
    do_loop     = opts.do_loop

    out        = open(output, 'w')
    waffle     = Unpickler(open(waffle_file, 'rb')).load()
    group      = list(waffle.keys())[0]
    size       = waffle[group]['size']
    for group, data in waffle.items():
        counter = data['counter']
        try:
            matrix  = [[data['sum_nrm'][i, j] / counter
                        for i in range(size)[::-1]]
                       for j in range(size)]
            x, y = matrix_to_decay(matrix, size, metric='loop' if do_loop else 'normal')
            spear, pval = spearmanr(x, y)
        except ZeroDivisionError:
            spear = pval = float('nan')
        center = get_center(matrix, size)
        test = pval < 0.05 and spear > 0 and center > 1
        out.write('{}\t{}\t{}\t{}\t{}\n'.format(group, spear, pval, center, test))
    out.close()


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='infile', required=True,
                        metavar='PATH', help='path to input file (pickle format)')
    parser.add_argument('-o', dest='outfile', required=True,
                        metavar='PATH', help='''path to output text file with
                        summary statistics by group''')
    parser.add_argument('--loop', dest='do_loop', action='store_true',
                        help='''correct for the formation of
                        tight loops between peaks''')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
