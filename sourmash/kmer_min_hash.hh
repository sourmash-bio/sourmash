#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <algorithm>
#include <iostream>
#include <set>
#include <map>
#include <queue>
#include <exception>
#include <string>

#include "../third-party/smhasher/MurmurHash3.h"

#define tbl \
  "                                                                "\
  /*ABCDEFGHIJKLMNOPQRSTUVWXYZ      abcdefghijklmnopqrstuvwxyz    */\
  " TVGH FCD  M KN   YSAABW R       TVGH FCD  M KN   YSAABW R"

uint64_t _hash_murmur(const std::string& kmer,
                      const uint32_t seed) {
    uint64_t out[2];
    out[0] = 0; out[1] = 0;
    MurmurHash3_x64_128((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
}

typedef uint64_t HashIntoType;

typedef std::vector<HashIntoType> CMinHashType;

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
    const bool dayhoff;
    const uint32_t seed;
    const HashIntoType max_hash;
    CMinHashType mins;

    KmerMinHash(unsigned int n, unsigned int k, bool prot, bool dyhoff, uint32_t s,
                HashIntoType mx)
        : num(n), ksize(k), is_protein(prot), dayhoff(dyhoff), seed(s), max_hash(mx) {
      if (n > 0) {
        mins.reserve(num + 1);
      }
      // only reserve a finite amount of space for unbounded MinHashes
      else {
        mins.reserve(1000);
      }
    };

    void check_compatible(const KmerMinHash& other) {
        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be compared");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be compared");
        }
        if (dayhoff != other.dayhoff) {
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
      if ((max_hash and h <= max_hash) or not max_hash) {
        if (mins.size() == 0) {
          mins.push_back(h);
          return;
        }
        else if (h <= max_hash or mins.back() > h or mins.size() < num) {
          auto pos = std::lower_bound(std::begin(mins), std::end(mins), h);

          // must still be growing, we know the list won't get too long
          if (pos == mins.cend()) {
            mins.push_back(h);
          }
          // inserting somewhere in the middle, if this value isn't already
          // in mins store it and shrink list if needed
          else if (*pos != h) {
            mins.insert(pos, h);
            if (num and mins.size() > num) {
              mins.pop_back();
            }
          }
        }
      }
    }

    virtual void remove_hash(const HashIntoType h) {
        auto pos = std::lower_bound(std::begin(mins), std::end(mins), h);
        if (pos != mins.cend() and *pos == h) {
          mins.erase(pos);
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
        std::string seq = sequence;
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        if (!is_protein) {
            for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
                const std::string kmer = seq.substr(i, ksize);
                if (! _checkdna(kmer)) {
                    if (force) {
                        continue;
                    } else {
                        std::string msg = "invalid DNA character in input k-mer: ";
                        msg += kmer;
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

    std::string translate_codon(std::string& codon) {
        std::string residue;

        if (codon.length() >= 2 && codon.length() <= 3){
            // If codon is length 2, pad with an N for ambiguous codon amino acids
            if (codon.length() == 2) {
                codon += "N";
            }
            auto translated = _codon_table.find(codon);

            if (translated != _codon_table.end()) {
                // "second" is the element mapped to by the codon
                // Because .find returns an iterator
                residue = translated -> second;
            } else {
                // Otherwise, assign the "X" or "unknown" amino acid
                residue = "X";
            }
        } else if (codon.length() == 1){
            // Then we only have one nucleotides and the amino acid is unknown
            residue = "X";
        } else {
            std::string msg = "Codon is invalid length: ";
            msg += codon;
            throw minhash_exception(msg);
        }
        return residue;
    }

    std::string _dna_to_aa(const std::string& dna) {
        std::string aa;
        std::string codon;
        std::string residue;
        unsigned int dna_size = (dna.size() / 3) * 3; // floor it
        for (unsigned int j = 0; j < dna_size; j += 3) {

            codon = dna.substr(j, 3);

            residue = translate_codon(codon);

            // Use dayhoff encoding of amino acids
            if (dayhoff) {
                std::string new_letter = aa_to_dayhoff(residue);
                aa += new_letter;
            } else {
                aa += residue;
            }

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

        auto from = out.begin();
        auto to = out.end();

        char c;
        for (to--; from <= to; from++, to--) {
            c = tbl[(int)*from];
            *from = tbl[(int)*to];
            *to = c;
        }

        return out;
    }

    std::string aa_to_dayhoff(const std::string& aa) const {
        // Convert an amino acid letter to dayhoff encoding
        std::string new_letter;

        auto dayhoff_encoded = _dayhoff_table.find(aa);
        if (dayhoff_encoded != _dayhoff_table.end()) {
            // "second" is the element mapped to by the codon
            // Because .find returns an iterator
            new_letter = dayhoff_encoded -> second;
        } else {
            // Otherwise, assign the "X" or "unknown" amino acid
            new_letter = "X";
        }
        return new_letter;
    }

    virtual void merge(const KmerMinHash& other) {
        check_compatible(other);

        CMinHashType merged;
        merged.reserve(other.mins.size() + mins.size());
        std::set_union(other.mins.begin(), other.mins.end(),
                       mins.begin(), mins.end(),
                       std::back_inserter(merged));
        if (merged.size() < num or !num) {
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

        {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"}, {"TCN", "S"},

        {"TAT", "Y"}, {"TAC", "Y"},
        {"TAA", "*"}, {"TAG", "*"},

        {"TGT", "C"}, {"TGC", "C"},
        {"TGA", "*"},
        {"TGG", "W"},

        {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"}, {"CTN", "L"},

        {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"}, {"CCN", "P"},

        {"CAT", "H"}, {"CAC", "H"},
        {"CAA", "Q"}, {"CAG", "Q"},

        {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"}, {"CGN", "R"},

        {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"},
        {"ATG", "M"},

        {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"}, {"ACN", "T"},

        {"AAT", "N"}, {"AAC", "N"},
        {"AAA", "K"}, {"AAG", "K"},

        {"AGT", "S"}, {"AGC", "S"},
        {"AGA", "R"}, {"AGG", "R"},

        {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"}, {"GTN", "V"},

        {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"}, {"GCN", "A"},

        {"GAT", "D"}, {"GAC", "D"},
        {"GAA", "E"}, {"GAG", "E"},

        {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}, {"GGN", "G"}
    };


// Dayhoff table from
// Peris, P., López, D., & Campos, M. (2008).
// IgTM: An algorithm to predict transmembrane domains and topology in
// proteins. BMC Bioinformatics, 9(1), 1029–11.
// http://doi.org/10.1186/1471-2105-9-367
//
// Original source:
// Dayhoff M. O., Schwartz R. M., Orcutt B. C. (1978).
// A model of evolutionary change in proteins,
// in Atlas of Protein Sequence and Structure,
// ed Dayhoff M. O., editor.
// (Washington, DC: National Biomedical Research Foundation; ), 345–352.
//
// | Amino acid    | Property              | Dayhoff |
// |---------------|-----------------------|---------|
// | C             | Sulfur polymerization | a       |
// | A, G, P, S, T | Small                 | b       |
// | D, E, N, Q    | Acid and amide        | c       |
// | H, K, R       | Basic                 | d       |
// | I, L, M, V    | Hydrophobic           | e       |
// | F, W, Y       | Aromatic              | f       |
    std::map<std::string, std::string> _dayhoff_table = {
        {"C", "a"},

        {"A", "b"}, {"G", "b"}, {"P", "b"}, {"S", "b"}, {"T", "b"},

        {"D", "c"}, {"E", "c"}, {"N", "c"}, {"Q", "c"},

        {"H", "d"}, {"K", "d"}, {"R", "d"},

        {"I", "e"}, {"L", "e"}, {"M", "e"}, {"V", "e"},

        {"F", "f"}, {"W", "f"}, {"Y", "f"}

    };

};

class KmerMinAbundance: public KmerMinHash {
 public:
    CMinHashType abunds;

    KmerMinAbundance(unsigned int n, unsigned int k, bool prot, bool dayhoff,
                     uint32_t seed, HashIntoType mx) :
        KmerMinHash(n, k, prot, dayhoff, seed, mx) { };

    virtual void add_hash(HashIntoType h) {
      if ((max_hash and h <= max_hash) or not max_hash) {
        // empty? add it, if within range / no range specified.
        if (mins.size() == 0) {
          mins.push_back(h);
          abunds.push_back(1);
          return;
        } else if (h <= max_hash or mins.back() > h or mins.size() < num) {
          // "good" hash - within range, smaller than current entry, or
          // still space.
          auto pos = std::lower_bound(std::begin(mins), std::end(mins), h);

          // at end -- must still be growing, we know the list won't get too
          // long
          if (pos == mins.cend()) {
            mins.push_back(h);
            abunds.push_back(1);
          } else if (*pos != h) {
          // didn't find hash already in mins, so
          // inserting somewhere in the middle; shrink list if needed.

            // calculate distance for use w/abunds *before* insert, as
            // 'mins.insert' may invalidate 'pos'.
            size_t dist = std::distance(begin(mins), pos);
            mins.insert(pos, h);
            abunds.insert(begin(abunds) + dist, 1);

            // now too big? if so, continue.
            if (mins.size() > num and not max_hash) {
              mins.pop_back();
              abunds.pop_back();
            }
          } else { // *pos == h - hash value already there, increment count.
            auto p = std::distance(begin(mins), pos);
            abunds[p] += 1;
          }
        }
      }
    }

    virtual void remove_hash(const HashIntoType h) {
        auto pos = std::lower_bound(std::begin(mins), std::end(mins), h);
        if (pos != mins.cend() and *pos == h) {
          mins.erase(pos);
          size_t dist = std::distance(begin(mins), pos);
          abunds.erase(begin(abunds) + dist);
        }
    }

    virtual void merge(const KmerMinAbundance& other) {
        check_compatible(other);

        CMinHashType merged_mins;
        CMinHashType merged_abunds;
        size_t max_size = other.mins.size() + mins.size();

        merged_mins.reserve(max_size);
        merged_abunds.reserve(max_size);

        auto it1_m = mins.begin();
        auto it2_m = other.mins.begin();
        auto out_m = std::back_inserter(merged_mins);

        auto it1_a = abunds.begin();
        auto it2_a = other.abunds.begin();
        auto out_a = std::back_inserter(merged_abunds);

        for (; it1_m != mins.end(); ++out_m, ++out_a) {
            if (it2_m == other.mins.end()) {
                /* we reached the end of other.mins,
                   so just copy the remainder of mins to the output */
                std::copy(it1_m, mins.end(), out_m);
                std::copy(it1_a, abunds.end(), out_a);
                break;
            }
            if (*it2_m < *it1_m) {
                /* other.mins is smaller than mins,
                   so copy it to output and advance other.mins iterators */
                *out_m = *it2_m;
                *out_a = *it2_a;
                ++it2_m;
                ++it2_a;
            } else if (*it2_m == *it1_m) {
                /* same value in both mins, so sums the abundances
                   on the output and advances all iterators */
                *out_m = *it1_m;
                *out_a = *it1_a + *it2_a;
                ++it1_m; ++it1_a;
                ++it2_m; ++it2_a;
            } else {
                /* mins is smaller than other.mins,
                   so copy it to output and advance the mins iterators */
                *out_m = *it1_m;
                *out_a = *it1_a;
                ++it1_m;
                ++it1_a;
            }
        }
        /* we reached the end of mins/abunds,
           so just copy the remainder of other to the output
           (other might already be at the end, in this case nothing happens) */
        std::copy(it2_m, other.mins.end(), out_m);
        std::copy(it2_a, other.abunds.end(), out_a);

        if (merged_mins.size() < num || !num) {
          mins = merged_mins;
          abunds = merged_abunds;
        } else {
          mins = CMinHashType(std::begin(merged_mins), std::begin(merged_mins) + num);
          abunds = CMinHashType(std::begin(merged_abunds), std::begin(merged_abunds) + num);
        }
    }

    virtual size_t size() {
        return mins.size();
    }

};

#endif // KMER_MIN_HASH_HH
