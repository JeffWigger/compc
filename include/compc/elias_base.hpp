#ifndef COMPC_ELIAS_BASE_H_
#define COMPC_ELIAS_BASE_H_
#include <cmath>
#include <cstdint>
#include <iostream>
#include <tuple>
#include <vector>

#include "compc/compressor.hpp"
#include "compc/helpers.hpp"
namespace compc {

typedef struct ArrayPrefixSummary {
  int local_threads;
  std::vector<std::size_t> local_sums;
  uint32_t batch_size;
  std::size_t total_chunks;
  bool error;
} ArrayPrefixSummary;

template <typename T> class EliasBase : public Compressor<T> {
public:
  EliasBase() = default;
  EliasBase(T zero_offset) : offset(zero_offset){};
  EliasBase(T zero_offset, bool map_negative_numbers_to_positive)
      : offset(zero_offset), map_negative_numbers(map_negative_numbers_to_positive){};
  virtual ~EliasBase() = default;
  virtual ArrayPrefixSummary get_prefix_sum_array(const T*, std::size_t) = 0;
  T offset{0};
  bool map_negative_numbers{false};
};
} // namespace compc

#endif // COMPC_ELIAS_BASE_H_
