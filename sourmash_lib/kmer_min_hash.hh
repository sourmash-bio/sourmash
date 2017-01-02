#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <exception>
#include <string>
#include <unordered_map>

// #include "kmer_hash.hh"

////

typedef std::set<HashIntoType> CMinHashType;

typedef std::map<HashIntoType, uint64_t> CMinAbundanceType;

class minhash_exception : public std::exception
{
public:
    explicit minhash_exception(const std::string& msg = "Generic minhash exception")
        : _msg(msg) { }

    virtual ~minhash_exception() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};


class KmerMinHash
{
public:
    const uint32_t seed;
    const unsigned int num;
    const unsigned int ksize;
    const bool is_protein;
    CMinHashType mins;

    KmerMinHash(unsigned int n, unsigned int k, bool prot, uint32_t seed) :
      num(n), ksize(k), is_protein(prot), seed(seed) { };

    virtual void _shrink() {
        while (mins.size() > num) {
            CMinHashType::iterator mi = mins.end();
            mi--;
            mins.erase(mi);
        }
    }
    virtual void add_hash(HashIntoType h) {
        mins.insert(h);
        _shrink();
    }
  void add_word(std::string word) {
    HashIntoType hash = _hash_murmur(word, seed);
        add_hash(hash);
    }
    void add_sequence(const char * sequence, bool force=false) {
        if (strlen(sequence) < ksize) {
            return;
        }
        const std::string seq = _checkdna(sequence, force);
        if (!is_protein) {
            for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
                std::string kmer = seq.substr(i, ksize);
                std::string rc = _revcomp(kmer);
                if (kmer < rc) {
                    add_word(kmer);
                } else {
                    add_word(rc);
                }
            }
        } else {                      // protein
            std::string rc = _revcomp(seq);
            for (unsigned int i = 0; i < 3; i++) {
                std::string aa = _dna_to_aa(seq.substr(i, seq.length() - i));
                unsigned int aa_ksize = int(ksize / 3);
                std::string kmer;

                for (unsigned int j = 0; j < aa.length() - aa_ksize + 1; j++) {
                    kmer = aa.substr(j, aa_ksize);
                    add_word(kmer);
                }

                aa = _dna_to_aa(rc.substr(i, rc.length() - i));
                aa_ksize = int(ksize / 3);

                for (unsigned int j = 0; j < aa.length() - aa_ksize + 1; j++) {
                    kmer = aa.substr(j, aa_ksize);
                    add_word(kmer);
                }
            }
        }
    }

    std::string _dna_to_aa(const std::string& dna) {
        std::string aa;
        unsigned int dna_size = (dna.size() / 3) * 3; // floor it
        for (unsigned int j = 0; j < dna_size; j += 3) {
            std::string codon = dna.substr(j, 3);
            aa += (_codon_table)[codon];
        }
        return aa;
    }

    std::string _checkdna(const char * s, bool force=false) const {
        size_t seqsize = strlen(s);

        std::string seq = s;

        for (size_t i=0; i < seqsize; ++i) {
            switch(seq[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                break;
            default:
                if (force) {
                    seq[i] = 'N';
                } else {
                    std::string msg = "invalid DNA character in sequence: ";
                    msg += seq[i];
                    throw minhash_exception(msg);
                }
                break;
            }
        }
        return seq;
    }

    std::string _revcomp(const std::string& kmer) const {
        std::string out = kmer;
        size_t ksize = out.size();

        for (size_t i=0; i < ksize; ++i) {
            char complement;

            switch(kmer[i]) {
            case 'A':
                complement = 'T';
                break;
            case 'C':
                complement = 'G';
                break;
            case 'G':
                complement = 'C';
                break;
            case 'T':
                complement = 'A';
                break;
            case 'N':
                complement = 'N';
                break;
            default:
                std::string msg = "invalid DNA character in sequence: ";
                msg += kmer[i];

                throw minhash_exception(msg);
            }
            out[ksize - i - 1] = complement;
        }
        return out;
    }

    virtual void merge(const KmerMinHash& other) {
        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be merged");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be merged");
        }
        for (auto mi: other.mins) {
            mins.insert(mi);
        }
        _shrink();
    }
    virtual unsigned int count_common(const KmerMinHash& other) {
        CMinHashType combined;

        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be compared");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be compared");
        }

        CMinHashType::iterator mi;
        for (mi = mins.begin(); mi != mins.end(); ++mi) {
            combined.insert(*mi);
        }
        for (mi = other.mins.begin(); mi != other.mins.end(); ++mi) {
            combined.insert(*mi);
        }
        return mins.size() + other.mins.size() - combined.size();
    }
    virtual ~KmerMinHash() throw() { }

private:
    std::map<std::string, std::string> _codon_table = {
        {"TTT", "F"}, {"TTC", "F"},
        {"TTA", "L"}, {"TTG", "L"},

        {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},

        {"TAT", "Y"}, {"TAC", "Y"},
        {"TAA", "*"}, {"TAG", "*"},

        {"TGT", "C"}, {"TGC", "C"},
        {"TGA", "*"},
        {"TGG", "W"},

        {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"},

        {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},

        {"CAT", "H"}, {"CAC", "H"},
        {"CAA", "Q"}, {"CAG", "Q"},

        {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},

        {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"},
        {"ATG", "M"},

        {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},

        {"AAT", "N"}, {"AAC", "N"},
        {"AAA", "K"}, {"AAG", "K"},

        {"AGT", "S"}, {"AGC", "S"},
        {"AGA", "R"}, {"AGG", "R"},

        {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},

        {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},

        {"GAT", "D"}, {"GAC", "D"},
        {"GAA", "E"}, {"GAG", "E"},

        {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}
    };
};

class KmerMinAbundance: public KmerMinHash {
 public:
    CMinAbundanceType mins;
    HashIntoType max_mins;

  KmerMinAbundance(unsigned int n, unsigned int k, bool prot, uint32_t seed) :
      KmerMinHash(n, k, prot, seed) { };

    virtual void add_hash(HashIntoType h) {
        if (mins.size() < num) {
            mins[h] += 1;
            max_mins = std::max(max_mins, h);
            return;
        }

        if (h > max_mins) {
            return;
        } else {
            if (mins.find(h) != mins.end()) {
                mins[h] += 1;
            } else {
                mins.emplace(h, 1);
                mins.erase(max_mins);
                max_mins = (*std::max_element(mins.begin(), mins.end())).first;
            }
        }
        _shrink();
    }

    virtual void _shrink() {
        while (mins.size() > num) {
            mins.erase(max_mins);
            max_mins = (*std::max_element(mins.begin(), mins.end())).first;
        }
    }

    virtual void merge(const KmerMinAbundance& other) {
        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be merged");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be merged");
        }
        for (auto mi: other.mins) {
            mins[mi.first] += mi.second;
            max_mins = std::max(mi.first, max_mins);
        }
        _shrink();
    }

    virtual unsigned int count_common(const KmerMinAbundance& other) {
        CMinHashType combined;

        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be compared");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be compared");
        }

        for (auto mi: mins) {
            combined.insert(mi.first);
        }
        for (auto mi: other.mins) {
            combined.insert(mi.first);
        }
        return mins.size() + other.mins.size() - combined.size();
    }

    virtual unsigned int count_common(const KmerMinHash& other) {
        CMinHashType combined;

        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be compared");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be compared");
        }

        for (auto mi: mins) {
            combined.insert(mi.first);
        }
        for (auto mi: other.mins) {
            combined.insert(mi);
        }
        return mins.size() + other.mins.size() - combined.size();
    }

};

#endif // KMER_MIN_HASH_HH
