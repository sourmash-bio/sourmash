from enum import Enum

class HashFunctions(Enum):
    murmur64_DNA = 1
    murmur64_protein = 2
    murmur64_dayhoff = 3

    @classmethod
    def from_string(cls, hash_str):
        if hash_str == "0.murmur64_DNA":
            return cls.murmur64_DNA
        elif hash_str == "0.murmur64_protein":
            return cls.murmur64_protein
        elif hash_str == "0.murmur64_dayhoff":
            return cls.murmur64_dayhoff
        else:
            raise Exception("unknown molecule type: {}".format(molecule))
