
import sys

from .signature import load_signatures
from .sourmash_args import traverse_find_sigs, get_moltype
from .sbtmh import SigLeaf, create_sbt_index
from .logging import notify, error


def build_sbt(signatures, ksize, moltype, scaled=False,
              downsample_n=False,
              bf_size=1e5, n_children=2):
    tree = create_sbt_index(bf_size, n_children=n_children)

    inp_files = traverse_find_sigs(signatures)

    n = 0
    ksizes = set()
    moltypes = set()
    nums = set()
    scaleds = set()

    for filename in inp_files:
        notify('loading {}', filename, end='\r')
        siglist = load_signatures(filename,
                                      ksize=ksize,
                                      select_moltype=moltype)
        siglist = list(siglist)
        if not siglist:
            notify \
                ('\nwarning: no signatures loaded at given ksize/molecule type from {}', filename)

        # load all matching signatures in this file
        ss = None
        for ss in siglist:
            ksizes.add(ss.minhash.ksize)
            moltypes.add(get_moltype(ss))
            nums.add(ss.minhash.num)

            if scaled:
                ss.minhash = ss.minhash.downsample_scaled(scaled)
            if downsample_n:
                ss.minhash = ss.minhash.downsample_n(downsample_n)
            scaleds.add(ss.minhash.scaled)

            leaf = SigLeaf(ss.md5sum(), ss)
            tree.add_node(leaf)
            n += 1
            if not ss:
                continue

            # check to make sure we aren't loading incompatible signatures
            if len(ksizes) > 1 or len(moltypes) > 1:
                error('multiple k-mer sizes or molecule types present; fail.')
                error('specify --dna/--protein and --ksize as necessary')
                error('ksizes: {}; moltypes: {}',
                      ", ".join(map(str, ksizes)), ", ".join(moltypes))
                sys.exit(-1)

            if nums == { 0 } and len(scaleds) == 1:
                pass # good
            elif scaleds == { 0 } and len(nums) == 1:
                pass # also good
            else:
                error('trying to build an SBT with incompatible signatures.')
                error('nums = {}; scaleds = {}', repr(nums), repr(scaleds))
                sys.exit(-1)

    return tree