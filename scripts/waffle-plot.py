#! /usr/bin/env python
"""
"""

from argparse    import ArgumentParser
try:  # python 3
    from pickle        import load
except ImportError:  # python 2
    from cPickle        import load

try:
    from meta_waffle.plots import plot_waffle
except ImportError:  # meta-waffle is not installed.. but it's still ok!!!
    from os.path  import join as os_join
    import sys

    sys.path.insert(0, os_join(os_split(os_split(__file__)[0])[0], 'meta_waffle'))
    from plots      import plot_waffle


def main():

    opts = get_options()

    waffle_file = opts.peak_file
    output      = opts.outfile
    title       = opts.title

    waffle = load(open(waffle_file, 'rb'))

    group = ''

    plot_waffle(waffle[group], title, output)


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
