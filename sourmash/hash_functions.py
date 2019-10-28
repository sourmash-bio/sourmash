from ._minhash import _HashFunctions as HashFunctions

def hashfunction_from_string(hash_str):
    if hash_str == "0.murmur64_DNA":
        return HashFunctions.murmur64_DNA
    elif hash_str == "0.murmur64_protein":
        return HashFunctions.murmur64_protein
    elif hash_str == "0.murmur64_dayhoff":
        return HashFunctions.murmur64_dayhoff
    else:
        raise Exception("unknown molecule type: {}".format(hash_str))
