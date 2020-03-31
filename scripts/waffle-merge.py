#! /usr/bin/env python
"""
"""
from argparse    import ArgumentParser
try:  # python 3
    from pickle        import load, dump, HIGHEST_PROTOCOL, _Unpickler as Unpickler
except ImportError:  # python 2
    from pickle        import load, dump, HIGHEST_PROTOCOL, Unpickler

from meta_waffle.utils import sum_groups


def main():
    opts = get_options()

    waffle_files = opts.infiles
    output       = opts.outfile

    waffle_tot = load(open(waffle_files[0], 'rb'))
    seen_groups = set(k for k in waffle_tot)
    dupl_groups = set()
    for waffle_file in waffle_files[1:]:
        waffle     = load(open(waffle_file, 'rb'))
        for k in waffle:
            if k in seen_groups:
                dupl_groups.add(k)
            else:
                seen_groups.add(k)
        sum_groups(waffle_tot, waffle)


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='infiles', required=True, nargs = '+',
                        metavar='PATH', help='path to input files (pickle format)')
    parser.add_argument('-o', dest='outfile', required=True,
                        metavar='PATH', help='''path to output pickle file with
                        the sum of all inputs''')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
