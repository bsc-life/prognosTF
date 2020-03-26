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
    from os.path import join as os_join, split as os_split
    import sys

    sys.path.insert(0, os_join(os_split(os_split(__file__)[0])[0], 'meta_waffle'))
    from plots      import plot_waffle


def main():

    opts = get_options()

    waffle_file = opts.infile
    output      = opts.outfile
    title       = opts.title
    do_loop     = opts.do_loop

    waffle = load(open(waffle_file, 'rb'))

    if opts.group:
        group = opts.group
        if not group in waffle:
            raise Exception('ERROR: {} not in waffle. Use one of:\n{}\n'.format(
                group, '\n'.join(waffle.keys())
            ))
    else:
        group = list(waffle.keys())[0]

    plot_waffle(waffle[group], title, output=output,
                metric='loop' if do_loop else 'normal')


def get_options():
    parser = ArgumentParser()

    parser.add_argument('-i', dest='infile', required=True,
                        metavar='PATH', help='''path to input pickle file''')
    parser.add_argument('-o', dest='outfile', required=True,
                        metavar='PATH', help='''path to output image (any format
                        based on file extension)''')
    parser.add_argument('--loop', dest='do_loop', action='store_true',
                        help='''correct for the formation of
                        tight loops between peaks''')
    parser.add_argument('--group', dest='group',
                        metavar='STR', help='''Plot a given sub-group (default
                        the first is shown)''')
    parser.add_argument('--title', dest='title', default='',
                        metavar='STR', help='''some quoted text to be used as
                        title for the plot''')

    opts = parser.parse_args()
    return opts

if __name__ == '__main__':
    exit(main())
