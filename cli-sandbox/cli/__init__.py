from argparse import ArgumentParser
import sys

# Commands
from . import compare
from . import compute
from . import gather
from . import info
from . import plot
from . import search

# Subcommand groups
from . import lca
from . import sbt
from . import signature


VERSION = '2.2.0'


class CustomParser(ArgumentParser):
    def __init__(self, citation=True, **kwargs):
        super(CustomParser, self).__init__(**kwargs)
        self.citation = citation
        self._citation_printed = False

    def print_citation(self):
        if self._citation_printed:
            return
        print('Cite me plz')
        self._citation_printed = True

    def parse_args(self, args=None, namespace=None):
        if (args is None and len(sys.argv) == 1) or (args is not None and len(args) == 0):
            self.print_help()
            raise SystemExit(0)
        args = super(CustomParser, self).parse_args(args=args, namespace=namespace)
        if ('quiet' not in args or not args.quiet) and self.citation:
            self.print_citation()
        return args


def add_moltype_args(parser):
    parser.add_argument('--protein', dest='protein', action='store_true',
                        help='choose a protein signature (default: False)')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false',
                        help='do not choose a protein signature')
    parser.set_defaults(protein=False)

    parser.add_argument('--dayhoff', dest='dayhoff', action='store_true',
                        help='build Dayhoff-encoded amino acid signatures (default: False)')
    parser.add_argument('--no-dayhoff', dest='dayhoff',
                        action='store_false',
                        help='do not build Dayhoff-encoded amino acid signatures')
    parser.set_defaults(dayhoff=False)

    parser.add_argument('--hp', '--hydrophobic-polar', dest='hp', action='store_true',
                        help='build hydrophobic-polar-encoded amino acid signatures (default: False)')
    parser.add_argument('--no-hp', '--no-hydrophobic-polar', dest='hp',
                        action='store_false',
                        help='do not build hydrophobic-polar-encoded amino acid signatures')
    parser.set_defaults(hp=False)

    parser.add_argument('--dna', '--rna', dest='dna', default=None,
                        action='store_true',
                        help='choose a nucleotide signature (default: True)')
    parser.add_argument('--no-dna', '--no-rna', dest='dna',
                        action='store_false',
                        help='do not choose a nucleotide signature')
    parser.set_defaults(dna=None)


def get_parser():
    commands = ['compute', 'compare', 'search', 'plot', 'gather', 'lca', 'sbt', 'info', 'signature']
    commandstr = ' -- '.join(sorted(commands))

    desc = 'Compute, compare, manipulate, and analyze MinHash sketches of DNA sequences.'
    parser = CustomParser(prog='sourmash', description=desc)
    parser._optionals.title = 'Options'
    parser.add_argument('-v', '--version', action='version', version='sourmash '+ VERSION)
    parser.add_argument('-q', '--quiet', action='store_true', help='don\'t print citation information')
    sub = parser.add_subparsers(
        title='Commands', dest='cmd', metavar='cmd', help=commandstr,
        description='Invoke "sourmash <cmd> --help" for more details on executing each command.'
    )
    for cmd in commands:
        getattr(sys.modules[__name__], cmd).subparser(sub)
    parser._action_groups.reverse()
    return parser
