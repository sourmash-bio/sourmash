def add_moltype_args(parser, default_dna=None):
    parser.add_argument('--protein', dest='protein', action='store_true')
    parser.add_argument('--no-protein', dest='protein',
                        action='store_false')
    parser.set_defaults(protein=False)

    parser.add_argument('--dna', dest='dna', default=None,
                        action='store_true')
    parser.add_argument('--no-dna', dest='dna', action='store_false')
    parser.set_defaults(dna=default_dna)

