from sourmash import signature
import sourmash_tst_utils as utils

def test_load_textmode(track_abundance):
    # ijson required a file in binary mode or bytes,
    # but we had an API example in the docs using 'rt'.
    # I fixed the docs, but I'm keeping this test here
    # to make sure we still support it =/
    sigfile = utils.get_test_data('genome-s10+s11.sig')
    with open(sigfile, 'rt') as sigfp:
        siglist = list(signature.load_signatures(sigfp))
    loaded_sig = siglist[0]
    assert loaded_sig.name == 'genome-s10+s11'
