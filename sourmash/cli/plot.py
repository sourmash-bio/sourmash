"""plot distance matrix made by 'compare'"""

def subparser(subparsers):
    subparser = subparsers.add_parser('plot')
    subparser.add_argument(
        'distances', help='output from "sourmash compare"'
    )
    subparser.add_argument(
        '--pdf', action='store_true',
        help='output PDF; default is PNG'
    )
    subparser.add_argument(
        '--labels', action='store_true',
        help='show sample labels on dendrogram/matrix'
    )
    subparser.add_argument(
        '--labeltext',
        help='filename containing list of labels (overrides signature names)'
    )
    subparser.add_argument(
        '--indices', action='store_false',
        help='show sample indices but not labels'
    )
    subparser.add_argument(
        '--vmin', default=0.0, type=float,
        help='lower limit of heatmap scale; default=%(default)f'
    )
    subparser.add_argument(
        '--vmax', default=1.0, type=float,
        help='upper limit of heatmap scale; default=%(default)f'
    )
    subparser.add_argument(
        '--subsample', type=int, metavar='N',
        help='randomly downsample to this many samples, max'
    )
    subparser.add_argument(
        '--subsample-seed', type=int, default=1, metavar='S',
        help='random seed for --subsample; default=1'
    )
    subparser.add_argument(
        '-f', '--force', action='store_true',
        help='forcibly plot non-distance matrices'
    )
    subparser.add_argument(
        '--output-dir', metavar='DIR', help='directory for output plots'
    )
    subparser.add_argument(
        '--csv', metavar='F',
        help='write clustered matrix and labels out in CSV format (with column'
        ' headers) to this file'
    )


def main(args):
    import sourmash
    return sourmash.commands.plot(args)
