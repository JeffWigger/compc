#ifndef COMPC_ELIAS_DELTA_H_
#define COMPC_ELIAS_DELTA_H_
#include <cstdint>
#include "compc/elias_base.hpp"
#include <tuple>
namespace compc
{

  template <typename T>
  class EliasDelta : public EliasBase<T>
  {
    public:
      EliasDelta() = default;
      EliasDelta(T zero_offset, bool map_negative_numbers_to_positive) : offset (zero_offset), map_negative_numbers (map_negative_numbers_to_positive){};
      ~EliasDelta() = default;
      uint8_t* compress(T*, std::size_t&) override;
      T* decompress(const uint8_t*, std::size_t, std::size_t) override;
      std::size_t get_compressed_length(const T*, std::size_t) override;
      ArrayPrefixSummary get_prefix_sum_array(const T*, std::size_t) override;
    private:
      uint32_t batch_size_small {100};
      uint32_t batch_size_large {100000};
  };
}

#endif  // COMPC_ELIAS_DELTA_H_
