#! /usr/bin/env python
"""
"""
from argparse    import ArgumentParser
from meta_waffle.utils import sum_groups


def main():
    opts = get_options()

    waffle_files = opts.infiles
    output       = opts.outfile
    verbose      = not opts.silent

    missed, counter = sum_groups(waffle_files, output,
                                 split_features=opts.split_features,
                                 clean=opts.clean, verbose=verbose)

    if verbose:
        good = len(waffle_files) - len(missed)
        print('Merged {:,} files out of {:,} (summing {:,} submatrices).'.format(
            good, len(waffle_files), sum(counter.values())))
        if missed:
            print('\n  - Missed:')
            for fnam in missed:
                print('    . {}'.format(fnam))


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='infiles', required=True, nargs = '+',
                        metavar='PATH', help='path to input files (pickle format)')
    parser.add_argument('-o', dest='outfile', required=True,
                        metavar='PATH', help='''path to output pickle file with
                        the sum of all inputs. WARNING: should be a directory if
                        "split-feature" is used''')
    parser.add_argument('--split-features', dest='split_features', default=False,
                        action='store_true', help='Generate one waffle per feature')
    parser.add_argument('--clean', dest='clean', default=False,
                        action='store_true', help='Remove input files')
    parser.add_argument('--silent', dest='silent', default=False,
                        action='store_true', help='''shhhhhhttt''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
