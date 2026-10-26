#ifndef GBWTGRAPH_MINIMAX_H
#define GBWTGRAPH_MINIMAX_H

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

/*
  minimax.h: k-mer counts for minimaxmer selection.

  Minimaxmers are from Yoshi Suzuki's himapper_rp (https://github.com/yoshihikosuzuki/himapper_rp):
  in each window, pick the k-mer with the fewest occurrences in the genome, then the most
  occurrences in the reads, then the smallest hash. Here the genome is the graph: genome_count is
  the number of distinct graph positions of the k-mer, the count that weighted minimizers threshold on.
  K-mers that are not in the graph are not candidates.

  The counts file is written by seed-eval's prototype/minimax_counts.cpp:
  "MMXCNT01", uint64 k, uint64 n, n uint64 canonical keys in sorted order, n uint32 values
  (genome_count << 16 | reads_count); little-endian. A canonical key is the smaller 2-bit key of
  the two strands.
*/

namespace gbwtgraph
{

//------------------------------------------------------------------------------

class MinimaxCounts
{
public:
  // genome_count is capped one below the 16-bit maximum, as in seed-eval's Python prototype, so that
  // the two order k-mers identically.
  constexpr static std::uint64_t GENOME_COUNT_MAX = 0xFFFE;
  constexpr static std::uint64_t READS_COUNT_MAX = 0xFFFF;
  // Bits of the key used for the bucket directory.
  constexpr static size_t BUCKET_BITS = 28;

  MinimaxCounts() = default;

  // Loads the counts file. Throws std::runtime_error if the file cannot be read or is invalid.
  explicit MinimaxCounts(const std::string& filename)
  {
    std::ifstream in(filename, std::ios::binary);
    if(!in) { throw std::runtime_error("MinimaxCounts: cannot open " + filename); }
    char magic[8];
    std::uint64_t k = 0, n = 0;
    in.read(magic, sizeof(magic));
    in.read(reinterpret_cast<char*>(&k), sizeof(k));
    in.read(reinterpret_cast<char*>(&n), sizeof(n));
    if(!in || std::string(magic, sizeof(magic)) != "MMXCNT01")
    {
      throw std::runtime_error("MinimaxCounts: " + filename + " is not a minimaxmer counts file");
    }
    if(k == 0 || k > 31) { throw std::runtime_error("MinimaxCounts: invalid k " + std::to_string(k) + " in " + filename); }
    this->kmer_length = k;
    this->keys.resize(n);
    this->values.resize(n);
    in.read(reinterpret_cast<char*>(this->keys.data()), n * sizeof(std::uint64_t));
    in.read(reinterpret_cast<char*>(this->values.data()), n * sizeof(std::uint32_t));
    if(!in || in.peek() != std::ifstream::traits_type::eof())
    {
      throw std::runtime_error("MinimaxCounts: " + filename + " has the wrong size for " + std::to_string(n) + " entries");
    }

    // Check the order, build the bucket directory, and compute a checksum for the fingerprint.
    size_t key_bits = 2 * k;
    this->shift = (key_bits > BUCKET_BITS ? key_bits - BUCKET_BITS : 0);
    this->bucket_start.assign((std::uint64_t(1) << (key_bits - this->shift)) + 1, 0);
    std::uint64_t checksum = 0xcbf29ce484222325ULL;
    for(size_t i = 0; i < n; i++)
    {
      if(i > 0 && this->keys[i] <= this->keys[i - 1])
      {
        throw std::runtime_error("MinimaxCounts: keys in " + filename + " are not sorted and unique");
      }
      this->bucket_start[(this->keys[i] >> this->shift) + 1]++;
      checksum = (checksum ^ this->keys[i]) * 0x100000001b3ULL;
      checksum = (checksum ^ this->values[i]) * 0x100000001b3ULL;
    }
    for(size_t b = 1; b < this->bucket_start.size(); b++) { this->bucket_start[b] += this->bucket_start[b - 1]; }

    std::ostringstream fingerprint;
    fingerprint << "k=" << k << ";n=" << n << ";fnv1a=" << std::hex << checksum;
    this->print = fingerprint.str();
  }

  size_t k() const { return this->kmer_length; }
  size_t size() const { return this->keys.size(); }

  // k, number of entries and a checksum of the contents.
  const std::string& fingerprint() const { return this->print; }

  // Finds the packed counts (genome_count << 16 | reads_count) of a canonical key.
  bool find(std::uint64_t canonical_key, std::uint32_t& packed) const
  {
    size_t bucket = canonical_key >> this->shift;
    if(bucket + 1 >= this->bucket_start.size()) { return false; }
    auto first = this->keys.begin() + this->bucket_start[bucket];
    auto last = this->keys.begin() + this->bucket_start[bucket + 1];
    auto iter = std::lower_bound(first, last, canonical_key);
    if(iter == last || *iter != canonical_key) { return false; }
    packed = this->values[iter - this->keys.begin()];
    return true;
  }

  /*
    Order value of a candidate: [genome_count:16][65535 - reads_count:16][hash >> 32], smaller is
    better. The high 32 bits of the hash break ties in the same order as the full hash, so with
    constant counts this is the plain minimizer order. Used for selection only: the low bits of a
    selected k-mer's order value are biased towards small values (it won its window), so they must
    not choose hash table cells; the index uses the key's own hash for that.
  */
  static size_t order(std::uint32_t packed, size_t hash)
  {
    std::uint64_t genome = std::min<std::uint64_t>(packed >> 16, GENOME_COUNT_MAX);
    std::uint64_t reads = std::min<std::uint64_t>(packed & 0xFFFF, READS_COUNT_MAX);
    return (genome << 48) | ((READS_COUNT_MAX - reads) << 32) | (hash >> 32);
  }

private:
  size_t kmer_length = 0;
  size_t shift = 0;
  std::vector<std::uint64_t> keys;
  std::vector<std::uint32_t> values;
  std::vector<std::uint64_t> bucket_start;
  std::string print;
};

//------------------------------------------------------------------------------

} // namespace gbwtgraph

#endif // GBWTGRAPH_MINIMAX_H
