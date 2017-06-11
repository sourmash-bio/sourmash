#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <algorithm>
#include <set>
#include <map>
#include <queue>
#include <exception>
#include <string>

#include "../third-party/smhasher/MurmurHash3.h"

uint64_t _hash_murmur(const std::string& kmer,
                      const uint32_t seed) {
    uint64_t out[2];
    out[0] = 0; out[1] = 0;
    MurmurHash3_x64_128((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
}

typedef uint64_t HashIntoType;

typedef std::vector<HashIntoType> CMinHashType;

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

// Looks like a iterator but all it does is counts push_backs
struct Counter {
  struct value_type {
    template <typename T> value_type(const T &) {}
  };
  void push_back(const value_type &) { ++count; }
  size_t count = 0;
};

class KmerMinHash
{
public:
    const unsigned int num;
    const unsigned int ksize;
    const bool is_protein;
    const uint32_t seed;
    const HashIntoType max_hash;
    CMinHashType mins;

    KmerMinHash(unsigned int n, unsigned int k, bool prot, uint32_t s,
                HashIntoType mx)
        : // overflow num to represent "no maximum"
          num(n > 0 ? n : -1),
          ksize(k), is_protein(prot), seed(s),
          // overflow max_hash to represent "no maximum", this simplifies
          // the comparison in add_hash()
          max_hash(mx > 0 ? mx : -1) {
      if (n > 0) {
        mins.reserve(num + 1);
      }
      // only reserve a finite amount of space for unbounded MinHashes
      else {
        mins.reserve(1000);
      }
    };

    virtual void _shrink() {
        // pass
    }

    void check_compatible(const KmerMinHash& other) {
        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be compared");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be compared");
        }
        if (max_hash != other.max_hash) {
            throw minhash_exception("mismatch in max_hash; comparison fail");
        }
        if (seed != other.seed) {
            throw minhash_exception("mismatch in seed; comparison fail");
        }
	}

    virtual void add_hash(const HashIntoType h) {
      if (h <= max_hash) {
        if (mins.size() == 0) {
          mins.push_back(h);
          return;
        }
        else if (mins.back() > h or mins.size() < num) {
          auto pos = std::lower_bound(std::begin(mins), std::end(mins), h);

          // must still be growing, we know the list won't get too long
          if (pos == mins.cend()) {
            mins.push_back(h);
          }
          // inserting somewhere in the middle, if this value isn't already
          // in mins store it and shrink list if needed
          else if (*pos != h) {
            mins.insert(pos, h);
            if (mins.size() > num) {
              mins.pop_back();
            }
          }
        }
      }
    }
    void add_word(const std::string& word) {
        const HashIntoType hash = _hash_murmur(word, seed);
        add_hash(hash);
    }
    void add_sequence(const char * sequence, bool force=false) {
        if (strlen(sequence) < ksize) {
            return;
        }
        const std::string seq = sequence;
        if (!is_protein) {
            for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
                const std::string kmer = seq.substr(i, ksize);
                if (! _checkdna(kmer)) {
                    if (force) {
                        continue;
                    } else {
                        std::string msg = "invalid DNA character in input: ";
                        msg += seq[i];
                        throw minhash_exception(msg);
                    }
                }

                const std::string rc = _revcomp(kmer);

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

    bool _checkdna(const std::string seq) const {

        for (size_t i=0; i < seq.length(); ++i) {
            switch(seq[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                break;
            default:
                return false;
            }
        }
	return true;
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
        check_compatible(other);

        CMinHashType merged;
        merged.reserve(other.mins.size() + mins.size());
        std::set_union(other.mins.begin(), other.mins.end(),
                       mins.begin(), mins.end(),
                       std::back_inserter(merged));
        if (merged.size() < num) {
          mins = merged;
        }
        else {
          mins = CMinHashType(std::begin(merged), std::begin(merged) + num);
        }
    }
    virtual unsigned int count_common(const KmerMinHash& other) {
        check_compatible(other);

        Counter counter;
        std::set_intersection(mins.begin(), mins.end(),
                              other.mins.begin(), other.mins.end(),
                              std::back_inserter(counter));
        return counter.count;
    }

    virtual size_t size() {
        return mins.size();
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
    CMinAbundanceType mins_abund;
    HashIntoType max_mins;

    KmerMinAbundance(unsigned int n, unsigned int k, bool prot, uint32_t seed,
                     HashIntoType mx) :
        KmerMinHash(n, k, prot, seed, mx) { };

    virtual void add_hash(HashIntoType h) {
        if (max_hash && h > max_hash) {
            return;
        }

        if (!num || mins_abund.size() < num) {
            mins_abund[h] += 1;
            max_mins = std::max(max_mins, h);
            return;
        }

        if (num && h > max_mins) {
            return;
        } else {
            if (mins_abund.find(h) != mins_abund.end()) {
                mins_abund[h] += 1;
            } else {
                mins_abund.emplace(h, 1);
                mins_abund.erase(max_mins);
                max_mins = (*std::max_element(mins_abund.begin(),
                                              mins_abund.end())).first;
            }
        }
        _shrink();
    }

    virtual void _shrink() {
        if (num == 0) {
            return;
        }
        while (mins_abund.size() > num) {
            mins_abund.erase(max_mins);
            max_mins = (*std::max_element(mins_abund.begin(),
                                          mins_abund.end())).first;
        }
    }

    virtual void merge(const KmerMinAbundance& other) {
        check_compatible(other);

        for (auto mi: other.mins_abund) {
            mins_abund[mi.first] += mi.second;
            max_mins = std::max(mi.first, max_mins);
        }
        _shrink();
    }

    virtual unsigned int count_common(const KmerMinAbundance& other) {
        std::set<HashIntoType> combined;

		check_compatible(other);

        for (auto mi: mins_abund) {
            combined.insert(mi.first);
        }
        for (auto mi: other.mins_abund) {
            combined.insert(mi.first);
        }
        return mins_abund.size() + other.mins_abund.size() - combined.size();
    }

    virtual unsigned int count_common(const KmerMinHash& other) {
        std::set<HashIntoType> combined;

		check_compatible(other);

        for (auto mi: mins_abund) {
            combined.insert(mi.first);
        }
        for (auto mi: other.mins) {
            combined.insert(mi);
        }
        return mins_abund.size() + other.mins.size() - combined.size();
    }

    virtual size_t size() {
        return mins_abund.size();
    }

};

#endif // KMER_MIN_HASH_HH
